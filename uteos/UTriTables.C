/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SUTTG, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#include <stdexcept>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <limits>

//
// EOS library defines
//
#include "EOSParam.H"
#include "EOSModel.H"

//
// build table defines
//
#include "Tables.H"
#include "TriMesh.H"
#include "RMesh.H"
#include "BMesh.H"
#include "mpi.h"
#include "UTriTables.H"
#include "ErrorMap.H"

UTriTables::UTriTables( std::vector<Parameters> & model_params,
			const std::string model_name,
			const std::string mesh_type,
			const std::vector<double> & xvar_bounds,
			const std::vector<double> & yvar_bounds,
			const int logxy,
			const int totalthreads,
			const int boundary_samples,
			const int region_samples,
                        const std::string & node_algorithm,
			const int add_nodes,
                        const double step_multiplier,
                        const double step_error_level,
			const int mpirank,
			const int mpisize,
			const int debug_level)
  : Tables(model_params,model_name,mesh_type,xvar_bounds,yvar_bounds,
           logxy,totalthreads,mpirank,mpisize), bsamples(boundary_samples),
    rsamples(region_samples), nodealg(node_algorithm), addnodes(add_nodes),
    stepmult(step_multiplier), steperr(step_error_level), debug(debug_level)
{
  //
  // size boundary meshes appropriately
  //
  bmst.resize(maxtab-mintab);
  bmse.resize(maxtab-mintab);

  //
  // sanitize optimizer nubers
  //
  if (stepmult < 0.1) stepmult = 0.1;
  if (stepmult > 1.0) stepmult = 1.0;
}

UTriTables::~UTriTables()
{
  //
  // cleanup any meshes
  //
  for (std::size_t i=0;i<bmst.size();i++)
    for (std::size_t j=0;j<bmst[i].size();j++)
      if (bmst[i][j] != 0) delete bmst[i][j];
  for (std::size_t i=0;i<bmse.size();i++)
    for (std::size_t j=0;j<bmse[i].size();j++)
      if (bmse[i][j] != 0) delete bmse[i][j];

}

class TriErrSorter
{
public:
  TriErrSorter( const int triangle,
                const int node1,
                const int node2,
                const double coord,
                const double maxerr,
                const double bcoord1,
                const double bcoord2 )
    : t(triangle),n1(node1),n2(node2),c(coord),e(maxerr),b(2,0.)
  {
    b[0] = bcoord1;
    b[1] = bcoord2;
  }

  int t;
  int n1,n2;
  double c;
  double e;
  std::vector<double> b;

};

bool tesvduplicatefun( TriErrSorter a,
                       TriErrSorter b )
{
  if (a.n1 < b.n1) return true;
  else if (a.n1 == b.n1) {
    if (a.n2 < b.n2) return true;
    else if (a.n2 == b.n2) {
      if (a.c < b.c) return true;
      else if (a.c == b.c) {
	if (a.e > b.e) return true;
	else return false;
      }
      else return false;
    }
    else return false;
  }
  else return false;
}

bool tesvmaxerrfun( TriErrSorter a,
                    TriErrSorter b )
{
  if (a.e > b.e) return true;
  else return false;
}

bool mtesortfun( std::pair<int,double> a,
                 std::pair<int,double> b )
{
  if (a.second > b.second) return true;
  else if (a.second < b.second) return false;
  else if (a.first < b.first) return true;
  else return false;
}

double UTriTables::compute_region_error( std::vector<RMesh *> & rmeshes,
                                         std::vector<int> & tnodes,
                                         std::vector<std::pair<int,std::vector<double> > > & maxerrtris )
{
  maxerrtris.clear();

  //
  // create triangles on master mesh
  //
  if (mpimaster) {
    rmeshes[0]->find_tris();
  }

  //
  // get and share triangles to recompute
  //
  std::vector<int> comptris;
  if (mpimaster) {
    rmeshes[0]->find_unsaved_errors(tnodes,comptris);
  }
  int cnsize = comptris.size();
  MPI_Bcast(&cnsize,1,MPI_INT,0,MPI_COMM_WORLD);
  if (!mpimaster) comptris.resize(cnsize);
  MPI_Bcast(comptris.data(),cnsize,MPI_INT,0,MPI_COMM_WORLD);

  //
  // compute errors on all meshes
  //
  double merr = -1.;
  int merrtri = -1;
  std::vector<double> merrloc(2,-1.);
  std::vector<double> trierrs;
  std::vector<double> trierrlocs;
  for (std::size_t i=0;i<rmeshes.size();i++) {
    std::vector<double> errloc(2);
    int errtri;
    double err = rmeshes[i]->compute_mesh_error_saved(comptris,errtri,errloc,rsamples);
    if (err > merr) {
      merr = err;
      merrtri = errtri;
      merrloc = errloc;
    }
    if (debug > 0) std::cout << "Rank " << myrank << " Region Mesh " << mintab+i 
			     << " err " << err << " tri " << merrtri 
			     << " loc x " << merrloc[0] << " y " << merrloc[1]
			     << std::endl;
    //
    // get triangles/error from this mesh
    //
    std::vector<Tri> curtris;
    rmeshes[i]->get_tris(curtris);

    //
    // default fill save vectors first time through
    //
    if (i < 1) {
      trierrs.assign(curtris.size(),0.);
      trierrlocs.assign(2*curtris.size(),1./3.);
    }

    //
    // save triangles with max error above tolerance
    //
    for (std::size_t j=0;j<curtris.size();j++) {
      if (curtris[j].maxerr > trierrs[j]) {
	trierrs[j] = curtris[j].maxerr;
	trierrlocs[2*j] = curtris[j].maxerrloc[0];
	trierrlocs[2*j+1] = curtris[j].maxerrloc[1];
      }
    }
  }
  if (debug > 0) std::cout << " rank " << myrank << " maxerr " << merr 
			   << " maxerrtri " << merrtri 
			   << " maxerrloc x " << merrloc[0] << " y " << merrloc[1]
			   << std::endl;

  //
  // communicate errors
  //
  int mtris = trierrs.size();
  if (mpimaster) {
    std::vector<double> me(mpicommsize*mtris);
    MPI_Gather(trierrs.data(),mtris,MPI_DOUBLE,me.data(),mtris,MPI_DOUBLE,0,MPI_COMM_WORLD);
    std::vector<double> ml(2*mpicommsize*mtris);
    MPI_Gather(trierrlocs.data(),2*mtris,MPI_DOUBLE,ml.data(),2*mtris,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if (debug > 0) {
      for (std::size_t i=0;i<me.size();i++) {
	std::cout << " recrank " << i << " maxerr " << me[i] << " maxerrtri " << i%me.size()
		  << " maxerrloc x " << ml[2*i] << " y " << ml[2*i+1] << std::endl;
      }
    }

    //
    // find max error for each tri
    //
    for (int i=0;i<mtris;i++) {
      for (int j=1;j<mpicommsize;j++) {
	if (me[i+j*mtris] > me[i]) {
	  me[i] = me[i+j*mtris];
	  ml[2*i] = ml[2*(i+j*mtris)];
	  ml[2*i+1] = ml[2*(i+j*mtris)+1];
	}
      }
    }

    //
    // find overall max error
    //
    merrtri = 0;
    merr = me[0];
    for (int i=1;i<mtris;i++) {
      if (me[i] > merr) {
	merrtri = i;
	merr = me[i];
      }
    }
    merrloc[0] = ml[2*merrtri];
    merrloc[1] = ml[2*merrtri+1];
    if (debug > 0) std::cout << " fmaxerr " << merr << " maxerrtri " << merrtri
			     << " maxerrloc x " << merrloc[0] << " y " << merrloc[1]
			     << std::endl;

    //
    // get triangle information
    //
    std::vector<Tri> tv;
    rmeshes[0]->get_tris(tv);

    //
    // prepare triangles for sorting
    //
    std::vector<TriErrSorter> tesv;
    for (int i=0;i<mtris;i++) {
      if (me[i] > 1.) {
	//
	// check error location
	//
	double b1 = ml[2*i];
	double b2 = ml[2*i+1];
	int n1,n2;
	if (b1 < 1.e-10) {
	  n1 = tv[i].n[1];
	  n2 = tv[i].n[2];
	  if (n1 < n2) {
	    tesv.push_back(TriErrSorter(i,n1,n2,b2,me[i],b1,b2));
	  }
	  else {
	    tesv.push_back(TriErrSorter(i,n2,n1,1.-b1-b2,me[i],b1,b2));
	  }
	}
	else if (b2 < 1.e-10) {
	  n1 = tv[i].n[2];
	  n2 = tv[i].n[0];
	  if (n1 < n2) {
	    tesv.push_back(TriErrSorter(i,n1,n2,1.-b1-b2,me[i],b1,b2));
	  }
	  else {
	    tesv.push_back(TriErrSorter(i,n2,n1,b1,me[i],b1,b2));
	  }
	}
	else if (1.-b1-b2 < 1.e-10) {
	  n1 = tv[i].n[0];
	  n2 = tv[i].n[1];
	  if (n1 < n2) {
	    tesv.push_back(TriErrSorter(i,n1,n2,b1,me[i],b1,b2));
	  }
	  else {
	    tesv.push_back(TriErrSorter(i,n2,n1,b2,me[i],b1,b2));
	  }
	}
	else {
	  //
	  // not on boundary so just fake the entries
	  //
	  tesv.push_back(TriErrSorter(i,-1,-1,0.,me[i],b1,b2));
	}
      }
    }
    std::sort(tesv.begin(),tesv.end(),tesvduplicatefun);

    //
    // Process sorted errors, removing duplicates and boundary locations
    //
    int nbn = rmeshes[0]->get_num_boundary_nodes();
    std::vector<TriErrSorter> tesv2;
    for (std::size_t i=0;i<tesv.size();i++) {
      //
      // accept interior points
      //
      if (tesv[i].n1 < 0) {
	tesv2.push_back(tesv[i]);
	continue;
      }
      //
      // reject boundary states
      //
      if (tesv[i].n1 < nbn && tesv[i].n2 < nbn &&
	  ((tesv[i].n1+1)%nbn == tesv[i].n2 || (tesv[i].n2+1)%nbn == tesv[i].n1)) {
	if (debug > 0) std::cerr << "UTriTables::compute_region_error: warning -- found unacceptable maxerr on triangle boundary: " << tesv[i].e << " tri " << tesv[i].t << std::endl;
	continue;
      }
      //
      // skip duplicates
      //
      if (tesv2.size() > 0) {
	if (tesv[i].n1 == tesv2.back().n1 &&
	    tesv[i].n2 == tesv2.back().n2 &&
	    fabs(tesv[i].c-tesv2.back().c) < 1.e-10*fabs(tesv[i].c+tesv2.back().c)) continue;
	else tesv2.push_back(tesv[i]);
      }
      else tesv2.push_back(tesv[i]);
    }

    //
    // Now sort by maxerr and set up return vector
    //
    std::sort(tesv2.begin(),tesv2.end(),tesvmaxerrfun);
    for (std::size_t i=0;i<tesv2.size();i++) {
      maxerrtris.push_back(std::pair<int,std::vector<double> >(tesv2[i].t,tesv2[i].b));
    }
  }
  else {
    //
    // communicate errors
    //
    MPI_Gather(trierrs.data(),mtris,MPI_DOUBLE,NULL,mtris,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(trierrlocs.data(),2*mtris,MPI_DOUBLE,NULL,2*mtris,MPI_DOUBLE,0,MPI_COMM_WORLD);
  }

  //
  // broadcast max error
  //
  MPI_Bcast(&merr,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  
  return merr;
}

double UTriTables::compute_region_error_initial( std::vector<RMesh *> & rmeshes,
                                                 std::vector<std::pair<int,double> > & tev )
{
  //
  // create triangles on master mesh
  //
  if (mpimaster) {
    rmeshes[0]->find_tris();
  }

  //
  // compute errors on all meshes
  //
  double merr = -1.;
  int merrtri = -1;
  std::vector<double> trierrs;
  for (std::size_t i=0;i<rmeshes.size();i++) {
    std::vector<double> errloc(2);
    int errtri;
    double err = rmeshes[i]->compute_mesh_error_thread(errtri,errloc,rsamples,1);
    if (err > merr) {
      merr = err;
      merrtri = errtri;
    }
    if (debug > 0) std::cout << "Rank " << myrank << " Region Mesh " << mintab+i 
			     << " err " << err << " tri " << merrtri 
			     << std::endl;
    //
    // get triangles/error from this mesh
    //
    std::vector<Tri> curtris;
    rmeshes[i]->get_tris(curtris);

    //
    // default fill save vectors first time through
    //
    if (i < 1) {
      trierrs.assign(curtris.size(),0.);
    }

    //
    // save max error for triangles
    //
    for (std::size_t j=0;j<curtris.size();j++) {
      if (curtris[j].maxerr > trierrs[j]) {
	trierrs[j] = curtris[j].maxerr;
      }
    }
  }
  if (debug > 0) std::cout << " rank " << myrank << " maxerr " << merr 
			   << " maxerrtri " << merrtri 
			   << std::endl;

  //
  // communicate errors
  //
  int mtris = trierrs.size();
  if (mpimaster) {
    std::vector<double> me(mpicommsize*mtris);
    MPI_Gather(trierrs.data(),mtris,MPI_DOUBLE,me.data(),mtris,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if (debug > 0) {
      for (std::size_t i=0;i<me.size();i++) {
	std::cout << " recrank " << i << " maxerr " << me[i] << " maxerrtri " << i%me.size()
		  << std::endl;
      }
    }

    //
    // find max error for each tri
    //
    for (int i=0;i<mtris;i++) {
      for (int j=1;j<mpicommsize;j++) {
	if (me[i+j*mtris] > me[i]) {
	  me[i] = me[i+j*mtris];
	}
      }
    }

    //
    // find overall max error
    //
    merrtri = 0;
    merr = me[0];
    for (int i=1;i<mtris;i++) {
      if (me[i] > merr) {
	merrtri = i;
	merr = me[i];
      }
    }
    if (debug > 0) std::cout << " fmaxerr " << merr << " maxerrtri " << merrtri
			     << std::endl;

    //
    // get triangle information
    //
    std::vector<Tri> tv;
    rmeshes[0]->get_tris(tv);

    //
    // Now enter into return vector
    //
    tev.resize(mtris);
    for (int i=0;i<mtris;i++) {
      tev[i].first = i;
      tev[i].second = me[i];
    }
    for (std::size_t i=0;i<tev.size();i++) {
      //std::cout << " tev i " << i << " err " << tev[i].second << std::endl;
    }
  }
  else {
    //
    // communicate errors
    //
    MPI_Gather(trierrs.data(),mtris,MPI_DOUBLE,NULL,mtris,MPI_DOUBLE,0,MPI_COMM_WORLD);
  }

  //
  // broadcast max error
  //
  MPI_Bcast(&merr,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  
  return merr;
}

double UTriTables::compute_region_error_initial( std::vector<RMesh *> & rmeshes,
                                                 std::vector<int> & tnodes,
                                                 std::vector<std::pair<int,double> > & tev )
{
  //
  // create triangles on master mesh
  //
  if (mpimaster) {
    rmeshes[0]->find_tris();
  }

  //
  // get and share triangles to recompute
  //
  std::vector<int> comptris;
  if (mpimaster) {
    rmeshes[0]->find_unsaved_errors(tnodes,comptris);
  }
  int cnsize = comptris.size();
  MPI_Bcast(&cnsize,1,MPI_INT,0,MPI_COMM_WORLD);
  if (!mpimaster) comptris.resize(cnsize);
  MPI_Bcast(comptris.data(),cnsize,MPI_INT,0,MPI_COMM_WORLD);

  //
  // compute errors on all meshes
  //
  double merr = -1.;
  int merrtri = -1;
  std::vector<double> trierrs;
  for (std::size_t i=0;i<rmeshes.size();i++) {
    std::vector<double> errloc(2);
    int errtri;
    double err = rmeshes[i]->compute_mesh_error_thread_saved(comptris,errtri,errloc,rsamples,1);
    if (err > merr) {
      merr = err;
      merrtri = errtri;
    }
    if (debug > 0) std::cout << "Rank " << myrank << " Region Mesh " << mintab+i 
			     << " err " << err << " tri " << merrtri 
			     << std::endl;
    //
    // get triangles/error from this mesh
    //
    std::vector<Tri> curtris;
    rmeshes[i]->get_tris(curtris);

    //
    // default fill save vectors first time through
    //
    if (i < 1) {
      trierrs.assign(curtris.size(),0.);
    }

    //
    // save max error for triangles
    //
    for (std::size_t j=0;j<curtris.size();j++) {
      if (curtris[j].maxerr > trierrs[j]) {
	trierrs[j] = curtris[j].maxerr;
      }
    }
  }
  if (debug > 0) std::cout << " rank " << myrank << " maxerr " << merr 
			   << " maxerrtri " << merrtri 
			   << std::endl;

  //
  // communicate errors
  //
  int mtris = trierrs.size();
  if (mpimaster) {
    std::vector<double> me(mpicommsize*mtris);
    MPI_Gather(trierrs.data(),mtris,MPI_DOUBLE,me.data(),mtris,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if (debug > 0) {
      for (std::size_t i=0;i<me.size();i++) {
	std::cout << " recrank " << i << " maxerr " << me[i] << " maxerrtri " << i%me.size()
		  << std::endl;
      }
    }

    //
    // find max error for each tri
    //
    for (int i=0;i<mtris;i++) {
      for (int j=1;j<mpicommsize;j++) {
	if (me[i+j*mtris] > me[i]) {
	  me[i] = me[i+j*mtris];
	}
      }
    }

    //
    // find overall max error
    //
    merrtri = 0;
    merr = me[0];
    for (int i=1;i<mtris;i++) {
      if (me[i] > merr) {
	merrtri = i;
	merr = me[i];
      }
    }
    if (debug > 0) std::cout << " fmaxerr " << merr << " maxerrtri " << merrtri
			     << std::endl;

    //
    // get triangle information
    //
    std::vector<Tri> tv;
    rmeshes[0]->get_tris(tv);

    //
    // Now enter into return vector
    //
    tev.resize(mtris);
    for (int i=0;i<mtris;i++) {
      tev[i].first = i;
      tev[i].second = me[i];
    }
    for (std::size_t i=0;i<tev.size();i++) {
      //std::cout << " tev i " << i << " err " << tev[i].second << std::endl;
    }
  }
  else {
    //
    // communicate errors
    //
    MPI_Gather(trierrs.data(),mtris,MPI_DOUBLE,NULL,mtris,MPI_DOUBLE,0,MPI_COMM_WORLD);
  }

  //
  // broadcast max error
  //
  MPI_Bcast(&merr,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  
  return merr;
}

double UTriTables::optimize_nodes( std::vector<RMesh *> & rmeshes,
				   RMeshMMOptData & rmmmod,
				   std::vector<std::pair<int,double> > & trierrs,
                                   const bool accurate,
                                   const double lastmaxerr )
{
  double time0 = MPI_Wtime();
  //
  // setup optimizer
  //
  UTriMMOpt utmmo(mpimaster,rmeshes,rmmmod,logvars,myrank,mpicommsize,rsamples);

  MCMinParams mcmp;
  double mult = 1.0;
  if (accurate || rmeshes[0]->get_num_nodes()%100 == 0) mult *= 2.0;
  if (lastmaxerr > steperr) mult *= stepmult;
  mcmp.maxsteps = int(150*mult*rmmmod.opt_nodes.size());
  mcmp.bestCostStop = 0.9;
  mcmp.costTemp = 0.01;
  mcmp.paramTempInit = 1.e10;
  mcmp.paramTempMax = 1.e50;
  mcmp.paramTempMult = 1.001;
  utmmo.initialize(mcmp);

  double time1 = MPI_Wtime();
  //
  // optimize nodes on the search domain
  //
  std::vector<double> state;
  utmmo.getInitState(state);
  double minerr = utmmo.minimize(state);

  double time2 = MPI_Wtime();
  //
  // for now redo the whole mesh to double check the error calc
  //

  //
  // save optimum node locations and new triangles
  //
  if (mpimaster) {
    int totn = rmeshes[0]->get_num_nodes();
    for (std::size_t i=0;i<rmmmod.opt_nodes.size();i++) {
      if (rmmmod.opt_nodes[i] < totn) {
        if (logvars) rmeshes[0]->move_node(rmmmod.opt_nodes[i],exp(state[2*i]),exp(state[2*i+1]));
        else rmeshes[0]->move_node(rmmmod.opt_nodes[i],state[2*i],state[2*i+1]);
      }
      else {
        if (logvars) rmeshes[0]->add_node(exp(state[2*i]),exp(state[2*i+1]));
        else rmeshes[0]->add_node(state[2*i],state[2*i+1]);
      }
    }
  }
  double time3 = MPI_Wtime();

  double maxerr = compute_region_error_initial(rmeshes,rmmmod.opt_nodes,trierrs);
  double time4 = MPI_Wtime();
  std::cout << "UTriTables::optimize_nodes: nn " << rmeshes[0]->get_num_nodes()
            << " sme " << minerr << " fme " << maxerr
            << " tini " << time1-time0 << " tmin " << time2-time1
            << " tmov " << time3-time2 << " terr " << time4-time3;

  return maxerr;
}

bool tbsortfun( std::pair<int,int> a,
                std::pair<int,int> b )
{
  int a1,a2;
  if (a.first < a.second) {
    a1 = a.first;
    a2 = a.second;
  }
  else {
    a1 = a.second;
    a2 = a.first;
  }
  int b1,b2;
  if (b.first < b.second) {
    b1 = b.first;
    b2 = b.second;
  }
  else {
    b1 = b.second;
    b2 = b.first;
  }
  if (a1 < b1) return true;
  else if (a1 == b1 && a2 < b2) return true;
  else return false;
}

bool dnbsortfun( std::pair<int,int> a,
                 std::pair<int,int> b )
{
  if (a.first < b.first) return true;
  else if (a.first > b.first) return false;
  else {
    if (a.second < b.second) return true;
    else return false;
  }
}

void UTriTables::setup_opt( std::vector<RMesh *> & rmeshes,
                            const int curtri,
                            std::vector<int> & lasttrinodes,
                            RMeshMMOptData & rmmmod )
{
  if (mpimaster) {
    //
    // get the nodes for the current triangle
    //
    std::vector<int> onodes;
    rmeshes[0]->get_tri_nodes(curtri,onodes);

    int min = 0;
    if (onodes[1] < onodes[min]) min = 1;
    if (onodes[2] < onodes[min]) min = 2;

    //
    // check if this is the same triangle as the last one and then save it as last
    //
    bool same = true;
    for (int i=0;i<3;i++) {
      if (lasttrinodes[i] != onodes[(min+i)%3]) same = false;
      lasttrinodes[i] = onodes[(min+i)%3];
    }

    //
    // only optimize non-boundary nodes
    //
    rmmmod.opt_nodes.clear();
    int nbn = rmeshes[0]->get_num_boundary_nodes();
    for (std::size_t i=0;i<onodes.size();i++)
      if (onodes[i] >= nbn) rmmmod.opt_nodes.push_back(onodes[i]);

    //
    // add a node if the last triangle is the same or there are no opt nodes
    //
    if (same || rmmmod.opt_nodes.size() < 1) {
      //
      // get the next node number to add
      //
      int nextnode = rmeshes[0]->get_num_nodes();
      for (int iii=0;iii<addnodes;iii++) {
	rmmmod.opt_nodes.push_back(nextnode+iii);
	onodes.push_back(nextnode+iii);
      }

      //
      // add node(s) to current tri
      //
      rmeshes[0]->add_nodes(curtri,addnodes);

      //
      // recompute the mesh
      //
      rmeshes[0]->find_tris();
    }

    //
    // get nodes and triangles for the mesh
    //
    std::vector<Tri> t;
    rmeshes[0]->get_tris(t);
    std::vector<Node> n;
    rmeshes[0]->get_nodes(n);

    //
    // find all triangles touching the current nodes (domain for node movement)
    //
    std::vector<int> t1s;
    for (std::size_t i=0;i<onodes.size();i++) t1s.insert(t1s.end(),n[onodes[i]].tris.begin(),n[onodes[i]].tris.end());
    std::sort(t1s.begin(),t1s.end());
    std::vector<int> t1;
    t1.push_back(t1s[0]);
    for (std::size_t i=1;i<t1s.size();i++) if (t1s[i] != t1.back()) t1.push_back(t1s[i]);

    //
    // get nodes of domain
    //
    std::vector<int> n1;
    for (std::size_t i=0;i<t1.size();i++) n1.insert(n1.end(),t[t1[i]].n,t[t1[i]].n+3);
    std::sort(n1.begin(),n1.end());
    
    //
    // get triangles adjacent to domain
    //
    std::vector<int> t2s;
    t2s.insert(t2s.end(),n[n1[0]].tris.begin(),n[n1[0]].tris.end());
    for (std::size_t i=1;i<n1.size();i++) {
      if (n1[i] != n1[i-1]) t2s.insert(t2s.end(),n[n1[i]].tris.begin(),n[n1[i]].tris.end());
    }
    std::sort(t2s.begin(),t2s.end());
    std::vector<int> t2;
    t2.push_back(t2s[0]);
    for (std::size_t i=1;i<t2s.size();i++) if (t2s[i] != t2.back()) t2.push_back(t2s[i]);
    
    //
    // get boundary of extended domain for delaunay triangulation
    //
    std::vector<std::pair<int,int> > tb;
    for (std::size_t i=0;i<t2.size();i++)
      for (int j=0;j<3;j++)
	tb.push_back(std::pair<int,int>(t[t2[i]].n[j],t[t2[i]].n[(j+1)%3]));
    std::sort(tb.begin(),tb.end(),tbsortfun);
    std::vector<std::pair<int,int> > tb2;
    for (std::size_t i=0;i<tb.size();i++) {
      if (i+1 < tb.size() && tb[i].first == tb[i+1].second && tb[i].second == tb[i+1].first) i++;
      else tb2.push_back(tb[i]);
    }

    //
    // order the extended domain boundary
    //
    std::vector<int> tbloops;
    tbloops.push_back(0);
    int curtbloop(0);
    for (std::size_t i=1;i<tb2.size();i++) {
      //
      // find the segment that matches the last segment node
      //
      int next=-1;
      for (std::size_t j=i;j<tb2.size();j++) {
        if (tb2[i-1].second == tb2[j].first) {
          next = j;
          break;
        }
      }
      if (next < 0) {
        if (tb2[i-1].second == tb2[tbloops[curtbloop]].first) {
          tbloops.push_back(i);
          curtbloop++;
          continue;
        }
        else {
          rmeshes[0]->print_tris();
          std::cout << "curtri " << curtri << std::endl;
          for (std::size_t ii=0;ii<t1.size();ii++) std::cout << "t1 i " << ii << " n " << t1[ii] << std::endl;
          for (std::size_t ii=0;ii<n1.size();ii++) std::cout << "n1 i " << ii << " n " << n1[ii] << std::endl;
          for (std::size_t ii=0;ii<t2.size();ii++) std::cout << "t2 i " << ii << " n " << t2[ii] << std::endl;
          for (std::size_t ii=0;ii<tb.size();ii++) std::cout << "tb i " << ii << " f " << tb[ii].first << " s " << tb[ii].second << std::endl;
          for (std::size_t ii=0;ii<tb2.size();ii++) std::cout << "tb2 i " << ii << " f " << tb2[ii].first << " s " << tb2[ii].second << std::endl;
          throw std::runtime_error("UTriTables::setup_opt: failed to find next boundary segment or close loop");
        }
      }
      
      //
      // swap with current location to connect boundary
      //
      std::swap(tb2[i],tb2[next]);
    }
    if (tb2[tbloops[curtbloop]].first != tb2.back().second) throw std::runtime_error("UTriTables::setup_opt: failed to connect boundary endpoints");
    tbloops.push_back(tb2.size());

    //
    // find longest loop
    //
    int lloop = 0;
    int llsize = tbloops[1]-tbloops[0];
    for (std::size_t i=1;i<tbloops.size()-1;i++) {
      int csize = tbloops[i+1]-tbloops[i];
      if (csize > llsize) {
        llsize = csize;
        lloop = i;
      }
    }

    //
    // save longest loop
    //
    std::vector<std::pair<int,int> > tb3(tb2.begin()+tbloops[lloop],tb2.begin()+tbloops[lloop+1]);

    //
    // check for duplicate nodes (forming loops) on the boundary
    //
    std::vector<std::pair<int,int> > dnb;
    for (std::size_t i=0;i<tb3.size();i++) dnb.push_back(std::pair<int,int>(tb3[i].first,i));
    std::sort(dnb.begin(),dnb.end(),dnbsortfun);
    std::vector<std::pair<int,int> > clipbn;
    for (std::size_t i=1;i<dnb.size();i++) {
      if (dnb[i].first == dnb[i-1].first) clipbn.push_back(std::pair<int,int>(dnb[i-1].second,dnb[i].second));
    }

    //
    // remove loops from boundary -- for now choose the smaller
    // possible loop, this is not robust in general
    //
    for (std::size_t i=0;i<clipbn.size();i++) {
      int o1 = clipbn[i].second-clipbn[i].first;
      int o2 = clipbn[i].first-clipbn[i].second+dnb.size();
      if (o1 < o2) {
        for (int j=clipbn[i].first;j<clipbn[i].second;j++) tb3[j].first = -1;
      }
      else {
        for (int j=0;j<clipbn[i].first;j++) tb3[j].first = -1;
        for (int j=clipbn[i].second;j<int(dnb.size());j++) tb3[j].first = -1;
      }
    }

    //
    // make the boundary
    //
    rmmmod.boundary_nodes.clear();
    for (std::size_t i=0;i<tb3.size();i++) if (tb3[i].first >= 0) rmmmod.boundary_nodes.push_back(tb3[i].first);

    //
    // get all nodes for the triangulation
    //
    std::vector<int> alln;
    for (std::size_t i=0;i<t2.size();i++) alln.insert(alln.end(),t[t2[i]].n,t[t2[i]].n+3);
    std::sort(alln.begin(),alln.end());
    std::vector<int> alln2;
    alln2.push_back(alln[0]);
    for (std::size_t i=1;i<alln.size();i++)
      if (alln[i] != alln2.back()) alln2.push_back(alln[i]);
    
    //
    // remove boundary nodes and opt nodes from all node list
    //
    for (std::size_t i=0;i<rmmmod.boundary_nodes.size();i++) {
      for (std::size_t j=0;j<alln2.size();j++)
        if (alln2[j] == rmmmod.boundary_nodes[i]) alln2[j] = -1;
    }
    for (std::size_t i=0;i<rmmmod.opt_nodes.size();i++) {
      for (std::size_t j=0;j<alln2.size();j++)
        if (alln2[j] == rmmmod.opt_nodes[i]) alln2[j] = -1;
    }

    //
    // save interior nodes
    //
    rmmmod.interior_nodes.clear();
    for (std::size_t i=0;i<alln2.size();i++)
      if (alln2[i] >= 0) rmmmod.interior_nodes.push_back(alln2[i]);

    //
    // eliminate domain triangles from extended domain list
    //
    std::vector<int> t2e;
    for (std::size_t i=0,j=0;i<t2.size();i++) {
      if (j == t1.size()) {
        for (std::size_t k=i;k<t2.size();k++) t2e.push_back(t2[k]);
        break;
      }
      if (t2[i] < t1[j]) t2e.push_back(t2[i]);
      else if (t2[i] == t1[j]) j++;
      else {
        while (j < t1.size()) {
          if (t2[i] < t1[j]) {
            t2e.push_back(t2[i]);
            break;
          }
          else if (t2[i] == t1[j]) {
            j++;
            break;
          }
          else j++;
        }
      }
    }

    //
    // save data
    //
    rmmmod.domain_tris.swap(t1);
    rmmmod.delaunay_tris.swap(t2e);

    //
    // compute delaunay tris in submesh space
    //
    rmmmod.deltris.clear();
    for (std::size_t i=0;i<rmmmod.delaunay_tris.size();i++) rmmmod.deltris.push_back(t[rmmmod.delaunay_tris[i]]);
    std::map<int,int> nodemap;
    int k=0;
    for (std::size_t i=0;i<rmmmod.boundary_nodes.size();i++,k++) nodemap[rmmmod.boundary_nodes[i]] = k;
    for (std::size_t i=0;i<rmmmod.interior_nodes.size();i++,k++) nodemap[rmmmod.interior_nodes[i]] = k;
    for (std::size_t i=0;i<rmmmod.deltris.size();i++) {
      for (int j=0;j<3;j++) rmmmod.deltris[i].n[j] = nodemap[rmmmod.deltris[i].n[j]];
    }
    //
    // reorder nodes to ensure smallest number is always first
    //
    for (std::size_t i=0;i<rmmmod.deltris.size();i++) {
      Tri & ct = rmmmod.deltris[i];
      int min = 0;
      if (ct.n[1] < ct.n[min]) min = 1;
      if (ct.n[2] < ct.n[min]) min = 2;
      if (min == 1) {
        std::swap(ct.n[0],ct.n[1]);
        std::swap(ct.n[1],ct.n[2]);
      }
      else if (min == 2) {
        std::swap(ct.n[0],ct.n[1]);
        std::swap(ct.n[0],ct.n[2]);
      }
    }

    //
    // compute width of domain for optimization
    //
    double xmin(n[n1[0]].d.inputs[0]);
    double xmax(n[n1[0]].d.inputs[0]);
    double ymin(n[n1[0]].d.inputs[1]);
    double ymax(n[n1[0]].d.inputs[1]);
    for (std::size_t i=1;i<n1.size();i++) {
      // skip duplicates
      if (n1[i] == n1[i-1]) continue;
      if (n[n1[i]].d.inputs[0] < xmin) xmin = n[n1[i]].d.inputs[0];
      if (n[n1[i]].d.inputs[0] > xmax) xmax = n[n1[i]].d.inputs[0];
      if (n[n1[i]].d.inputs[1] < ymin) ymin = n[n1[i]].d.inputs[1];
      if (n[n1[i]].d.inputs[1] > ymax) ymax = n[n1[i]].d.inputs[1];
    }
    if (logvars) {
      rmmmod.xywidth[0] = log(xmax)-log(xmin);
      rmmmod.xywidth[1] = log(ymax)-log(ymin);
    }
    else {
      rmmmod.xywidth[0] = xmax-xmin;
      rmmmod.xywidth[1] = ymax-ymin;
    }

    //
    // compute the domain tris for use in check state function
    //
    rmmmod.domtris.clear();
    for (std::size_t i=0;i<rmmmod.domain_tris.size();i++) rmmmod.domtris.push_back(t[rmmmod.domain_tris[i]]);
    for (std::size_t i=0;i<rmmmod.domtris.size();i++)
      rmmmod.domtris[i].computeBcoords(n,logvars);
    
    
  }

  //
  // communicate results
  //
  rmmmod.opt_tri = curtri;
  MPI_Bcast(lasttrinodes.data(),3,MPI_INT,0,MPI_COMM_WORLD);
  int dsize = rmmmod.opt_nodes.size();
  MPI_Bcast(&dsize,1,MPI_INT,0,MPI_COMM_WORLD);
  if (!mpimaster) rmmmod.opt_nodes.resize(dsize);
  MPI_Bcast(rmmmod.opt_nodes.data(),dsize,MPI_INT,0,MPI_COMM_WORLD);
  dsize = rmmmod.boundary_nodes.size();
  MPI_Bcast(&dsize,1,MPI_INT,0,MPI_COMM_WORLD);
  if (!mpimaster) rmmmod.boundary_nodes.resize(dsize);
  MPI_Bcast(rmmmod.boundary_nodes.data(),dsize,MPI_INT,0,MPI_COMM_WORLD);
  dsize = rmmmod.interior_nodes.size();
  MPI_Bcast(&dsize,1,MPI_INT,0,MPI_COMM_WORLD);
  if (!mpimaster) rmmmod.interior_nodes.resize(dsize);
  MPI_Bcast(rmmmod.interior_nodes.data(),dsize,MPI_INT,0,MPI_COMM_WORLD);
  dsize = rmmmod.domain_tris.size();
  MPI_Bcast(&dsize,1,MPI_INT,0,MPI_COMM_WORLD);
  if (!mpimaster) rmmmod.domain_tris.resize(dsize);
  MPI_Bcast(rmmmod.domain_tris.data(),dsize,MPI_INT,0,MPI_COMM_WORLD);
  dsize = rmmmod.delaunay_tris.size();
  MPI_Bcast(&dsize,1,MPI_INT,0,MPI_COMM_WORLD);
  if (!mpimaster) rmmmod.delaunay_tris.resize(dsize);
  MPI_Bcast(rmmmod.delaunay_tris.data(),dsize,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(rmmmod.xywidth,2,MPI_DOUBLE,0,MPI_COMM_WORLD);

  if (debug > 1) {
    std::cout << "opt_nodes: ";
    for (std::size_t i=0;i<rmmmod.opt_nodes.size();i++) std::cout <<rmmmod. opt_nodes[i] << " ";
    std::cout << std::endl;
    std::cout << "boundary_nodes: ";
    for (std::size_t i=0;i<rmmmod.boundary_nodes.size();i++) std::cout << rmmmod.boundary_nodes[i] << " ";
    std::cout << std::endl;
    std::cout << "interior_nodes: ";
    for (std::size_t i=0;i<rmmmod.interior_nodes.size();i++) std::cout << rmmmod.interior_nodes[i] << " ";
    std::cout << std::endl;
    std::cout << "domain_tris: ";
    for (std::size_t i=0;i<rmmmod.domain_tris.size();i++) std::cout << rmmmod.domain_tris[i] << " ";
    std::cout << std::endl;
    std::cout << "delaunay_tris: ";
    for (std::size_t i=0;i<rmmmod.delaunay_tris.size();i++) std::cout << rmmmod.delaunay_tris[i] << " ";
    std::cout << std::endl;
    std::cout << "errtris: ";
    for (std::size_t i=0;i<rmmmod.errtris.size();i++) std::cout << rmmmod.errtris[i] << " ";
    std::cout << std::endl;
    std::cout << "xwidth " << rmmmod.xywidth[0] << " ywidth " << rmmmod.xywidth[1] << std::endl;
    std::cout << "opt_tri " << rmmmod.opt_tri << std::endl;
  }

}

double UTriTables::refine_regions_pt( std::vector<std::vector<TMesh *> > & rms,
                                      int region )
{
  //
  // the region to refine: do the last one if unspecified
  //
  int r = rms[0].size()-1;
  if (region >= 0 && region < r) r = region;

  //
  // check region meshes all have the same size
  //
  int rsizempi = rms[0].size();
  MPI_Bcast(&rsizempi,1,MPI_INT,0,MPI_COMM_WORLD);
  std::size_t rsize = rsizempi;
  for (std::size_t i=1;i<rms.size();i++) {
    if (rsize != rms[i].size()) {
      throw std::runtime_error("regions mesh vectors do not have the same sizes in UTriTables::refine_regions_pt");
    }
  }

  //
  // sanity checks -- only boundary nodes and the same number in the region to refine
  //
  rsizempi = rms[0][r]->get_num_boundary_nodes();
  MPI_Bcast(&rsizempi,1,MPI_INT,0,MPI_COMM_WORLD);
  for (std::size_t i=0;i<rms.size();i++) {
    if (rms[i][r]->get_num_boundary_nodes() != rsizempi) {
      throw std::runtime_error("UTriTables::refine_regions_pt method must start with identical numbers of boundary points");
    }
    if (rms[i][r]->get_num_nodes() != rsizempi) {
      throw std::runtime_error("UTriTables::refine_regions_pt method must start with only boundary nodes");
    }
  }

  //
  // Set up vector of RMeshes for desired region
  //
  std::vector<RMesh *> rmeshes;
  for (std::size_t i=0;i<rms.size();i++) rmeshes.push_back(dynamic_cast<RMesh *>(rms[i][r]));

  //
  // start evaluation threads
  //
  if (totthr > 0)
    for (std::size_t i=0;i<rms.size();i++) rmeshes[i]->start_eval_threads();

  //
  // Compute initial error on all mesh triangles
  //
  std::vector<std::pair<int,double> > trierrs;
  double maxerr = compute_region_error_initial(rmeshes,trierrs);
  double besterr(maxerr);
  bool use_multiplier(true);

  //
  // Loop to reduce max error below 1
  //
  std::vector<int> lasttrinodes(3,-1);
  while (maxerr > 1.) {
    double time1 = MPI_Wtime();
    //
    // Choose max err triangle
    //
    int curtri;
    if (mpimaster) {
      //
      // trierrs only valid on master
      //
      std::sort(trierrs.begin(),trierrs.end(),mtesortfun);
      curtri = trierrs[0].first;
    }
    MPI_Bcast(&curtri,1,MPI_INT,0,MPI_COMM_WORLD);
    
    double time2 = MPI_Wtime();
    //
    // setup the optimization data
    //
    RMeshMMOptData rmmmod;
    setup_opt(rmeshes,curtri,lasttrinodes,rmmmod);
    double time3 = MPI_Wtime();

    //
    // optimize node locations to minimize error
    //
    maxerr = optimize_nodes(rmeshes,rmmmod,trierrs,use_multiplier,maxerr);
    if (maxerr < besterr) {
      besterr = maxerr;
      use_multiplier = true;
    }
    else use_multiplier = false;
    double time4 = MPI_Wtime();
    std::cout << " bme " << besterr << " ttri " << time2-time1
              << " tset " << time3-time2 << " topt " << time4-time3 << std::endl;
  }

  //
  // stop evaluation threads
  //
  if (totthr > 0)
    for (std::size_t i=0;i<rms.size();i++) rmeshes[i]->stop_eval_threads();

  return maxerr;
}

double UTriTables::refine_regions( std::vector<std::vector<TMesh *> > & rms,
                                   int region )
{
  //
  // the region to refine: do the last one if unspecified
  //
  int r = rms[0].size()-1;
  if (region >= 0 && region < r) r = region;

  //
  // check region meshes all have the same size
  //
  int rsizempi = rms[0].size();
  MPI_Bcast(&rsizempi,1,MPI_INT,0,MPI_COMM_WORLD);
  std::size_t rsize = rsizempi;
  for (std::size_t i=1;i<rms.size();i++) {
    if (rsize != rms[i].size()) {
      throw std::runtime_error("regions mesh vectors do not have the same sizes in UTriTables::refine_regions");
    }
  }

  //
  // sanity checks -- only boundary nodes and the same number in the region to refine
  //
  rsizempi = rms[0][r]->get_num_boundary_nodes();
  MPI_Bcast(&rsizempi,1,MPI_INT,0,MPI_COMM_WORLD);
  for (std::size_t i=0;i<rms.size();i++) {
    if (rms[i][r]->get_num_boundary_nodes() != rsizempi) {
      throw std::runtime_error("UTriTables::refine_regions method must start with identical numbers of boundary points");
    }
    if (rms[i][r]->get_num_nodes() != rsizempi) {
      throw std::runtime_error("UTriTables::refine_regions method must start with only boundary nodes");
    }
  }

  //
  // Set up vector of RMeshes for desired region
  //
  std::vector<RMesh *> rmeshes;
  for (std::size_t i=0;i<rms.size();i++) rmeshes.push_back(dynamic_cast<RMesh *>(rms[i][r]));

  //
  // setup reference mesh and node index pointer for it
  //
  if (mpimaster) rmeshes[0]->find_tris();
  std::vector<int> tnodes;

  std::vector<std::pair<int,std::vector<double> > > maxerrtris;
  double maxerr = compute_region_error(rmeshes,tnodes,maxerrtris);
  double bmaxerr(maxerr);

  int k=0;
  if (mpimaster) std::cout << "region " << r << " nodes " << k << " err " << maxerr << " berr " << bmaxerr << std::endl;

  //
  // add nodes until error is below desired tolerance
  //
  while (maxerr > 1.) {
    int done = 0;

    //
    // nodes are added to master mesh
    //
    if (mpimaster) {
      if (maxerrtris.size() == 0) done = 1;

      //
      // only add a max of addnodes nodes per iteration
      //
      int startnode = rmeshes[0]->get_num_nodes();
      for (int i=0;i<std::min(addnodes,int(maxerrtris.size()));i++,k++) rmeshes[0]->add_node(maxerrtris[i].first,maxerrtris[i].second);
      int endnode = rmeshes[0]->get_num_nodes();
      tnodes.clear();
      for (int i=startnode;i<endnode;i++) tnodes.push_back(i);
    }

    MPI_Bcast(&done,1,MPI_INT,0,MPI_COMM_WORLD);
    
    if (done) {
      if (mpimaster) std::cerr << "UtriTables::refine_regions: failed to reach error target likely due to boundary error values" << std::endl;
      break;
    }

    //
    // recompute error
    //
    maxerr = compute_region_error(rmeshes,tnodes,maxerrtris);
    if (maxerr < bmaxerr) bmaxerr = maxerr;

    if (mpimaster) std::cout << "region " << r << " nodes " << k
                             << " err " << maxerr << " berr " << bmaxerr << std::endl;
  }

  return maxerr;
}

double UTriTables::setnodelocs( std::vector<std::vector<BMesh *> > & bms,
                                const int b,
                                double rtol,
                                int nnodes )
{
  //
  // clear all but edge nodes
  //
  for (std::size_t i=0;i<bms.size();i++) bms[i][b]->clear_nodes();

  //int debug = 1;
  int curn = 1;
  double err,curerr;
  double lastx = 0.;
  while(1) {
    //
    // compute error in remaining piece
    //
    err = 0;
    for (std::size_t j=0;j<bms.size();j++) {
      curerr = bms[j][b]->setnodeloc(curn+1,1.0,bsamples);
      if (curerr > err) err = curerr;
      if (debug > 0) std::cerr << "UTriTables::setnodelocs 1 curerr " << curerr << " err " << err << std::endl;
    }

    //
    // maximize on the master
    //
    if (mpimaster) MPI_Reduce(MPI_IN_PLACE,&err,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    else MPI_Reduce(&err,NULL,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    if (debug > 0) std::cerr << "UTriTables::setnodelocs 1 master err " << err << std::endl;

    //
    // done when the error tolerance or the max number of nodes is reached
    //
    int done = 0;
    if (mpimaster) {
      if ((err <= rtol && nnodes == 0) || (curn >= nnodes-1 && nnodes > 0)) done = 1;
    }
    MPI_Bcast(&done,1,MPI_INT,0,MPI_COMM_WORLD);
    if (done==1) break;

    //
    // add a node and move it to meet the given tolerance
    //

    double xmin = lastx;
    double xmax = 1.0;
    double llerr = err;
    double llx = xmax;
    double lx = (xmin+xmax)/2.;

    //
    // sync process for next setnodeloc calls
    //
    MPI_Bcast(&lx,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

    //
    // add the node at the midpoint of the remaining distance and get new error
    //
    double lerr = 0;
    for (std::size_t j=0;j<bms.size();j++) {
      curerr = bms[j][b]->setnodeloc(curn,lx,bsamples);
      if (curerr > lerr) lerr = curerr;
      if (debug > 0) std::cerr << "UTriTables::setnodelocs 2 curerr " << curerr << " err " << err << std::endl;
    }

    //
    // maximize on the master
    //
    if (mpimaster) MPI_Reduce(MPI_IN_PLACE,&lerr,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    else MPI_Reduce(&lerr,NULL,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    if (debug > 0) std::cerr << "UTriTables::setnodelocs 1 master err " << lerr << std::endl;

    //
    // update bounds
    //
    if (lerr > rtol) xmax = lx;
    else xmin = lx;

    int i=0;

    if (debug > 0) std::cerr << "curn " << curn << " i " << i << " xmin " << xmin << " xmax " << xmax << " lx " << lx << " llx " << llx << " lerr " << lerr << " llerr " << llerr << std::endl;

    //
    // move node in a bounded secant search
    // master does the logic, others only run setnodeloc
    //
    while (1) {
      done = 0;
      if (mpimaster) {
	if (fabs(lerr-rtol) < 1.e-7*fabs(lerr+rtol)) done = 1;
	else if (fabs(lx-llx) < 1.e-14*fabs(lx+llx) || i > 100) done = 1;
      }
      MPI_Bcast(&done,1,MPI_INT,0,MPI_COMM_WORLD);
      if (done == 1) break;

      i++;

      //
      // make new estimate and bisect if out of bounds
      //
      double x = 1.;
      if (mpimaster) {
	x = lx+(rtol-lerr)*(lx-llx)/(lerr-llerr);
	if (debug > 0) std::cerr << " x " << x << " lt " << (x<=xmin) << " gt " << (x>=xmax) << std::endl;
	if (x >= xmax || x <= xmin || (lerr-llerr) == 0.) {
	  if (fabs(lerr-llerr)<1.e-8*fabs(lerr+llerr) || (lerr > rtol && llerr > rtol)) x = (xmax+99.*xmin)/100.;
	  else x = (xmax+xmin)/2.;
	}
	if (debug > 0) std::cerr << " x " << x << std::endl;
      }
      MPI_Bcast(&x,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      //
      // adjust node location and compute new error
      //
      err = 0;
      for (std::size_t j=0;j<bms.size();j++) {
	curerr = bms[j][b]->setnodeloc(curn,x,bsamples);
	if (curerr > err) err = curerr;
	if (debug > 0) std::cerr << "UTriTables::setnodelocs 3 curerr " << curerr << " err " << err << std::endl;
      }

      //
      // maximize on the master
      //
      if (mpimaster) MPI_Reduce(MPI_IN_PLACE,&err,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
      else MPI_Reduce(&err,NULL,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
      if (debug > 0) std::cerr << "UTriTables::setnodelocs 1 master err " << err << std::endl;

      if (mpimaster) {
	//
	// update bounds
	//
	if (err > rtol) xmax = x;
	else xmin = x;

	//
	// update saved results for secant search
	//
	llerr = lerr;
	lerr = err;
	llx = lx;
	lx = x;

	if (debug > 0) std::cerr << "curn " << curn << " i " << i << " xmin " << xmin << " xmax " << xmax << " lx " << lx << " llx " << llx << " lerr " << lerr << " llerr " << llerr << " derr " << (lerr-llerr)/lerr << std::endl;
      }
    }

    lastx = lx;

    curn++;

  }

  return err;

}

double UTriTables::refine_boundaries( std::vector<std::vector<BMesh *> > & bms,
                                      bool uniformerr,
                                      int boundary,
                                      int debug )
{
  //
  // check boundary meshes all have the same size
  //
  int bsizempi = bms[0].size();
  MPI_Bcast(&bsizempi,1,MPI_INT,0,MPI_COMM_WORLD);
  std::size_t bsize = bsizempi;
  for (std::size_t i=1;i<bms.size();i++) {
    if (bsize != bms[i].size()) {
      throw std::runtime_error("boundary mesh vectors do not have the same sizes in UTriTables::refine_boundaries");
    }
  }

  //
  // the boundary to refine: do the last one if unspecified
  //
  int b = bms[0].size()-1;
  if (boundary >= 0 && boundary < b) b = boundary;

  //debug = 1;
  int numerrcalcs = 0;

  //
  // sanity checks
  //
  for (std::size_t i=0;i<bms.size();i++) {
    if (bms[i][b]->get_num_nodes() != 2) {
      throw std::runtime_error("UTriTables::refine_boundaries method must start with two points in each boundary");
    }
  }

  //
  // start evaluation threads
  //
  if (totthr > 0)
    for (std::size_t i=0;i<bms.size();i++) bms[i][b]->start_eval_threads();

  //
  // initialize node count to 0 to insert as many as needed
  //
  int numnodes = 0;

  //
  // compute initial node locations
  //
  double err = setnodelocs(bms,b,1.0,numnodes);

  if (debug > 0) std::cerr << "initial error " << err << " nodes " << numnodes << std::endl;

  if (err > 1.0) {
    throw std::runtime_error("UTriTables::refine_boundaries failed to meet tolerance with requested nodes");
  }

  //
  // fix number of nodes
  //
  numnodes = bms[0][b]->get_num_nodes();

  //
  // exit if we dont want to try for a uniform error across all intervals
  //
  if (!uniformerr) {
    //
    // stop evaluation threads
    //
    if (totthr > 0)
      for (std::size_t i=0;i<bms.size();i++) bms[i][b]->stop_eval_threads();

    return 1.0;
  }

  //
  // use secant/bisection algorithm to find uniform error placement of nodes
  //
  double tmin = err;
  double tmax = 1.0;
  double lldiff = err-1.0;
  double llt = 1.0;
  double lt = (tmin+tmax)/2.;

  //
  // best tolerance
  //
  double bestt = 1.0;
  double bestld = lldiff;

  //
  // sync processes on next node location
  //
  MPI_Bcast(&lt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

  //
  // compute node location for midpoint of tolerance bounds
  //
  err = setnodelocs(bms,b,lt,numnodes);
  double ldiff = err-lt;

  //
  // update best value
  //
  if (fabs(ldiff) < fabs(bestld)) {
    bestld = ldiff;
    bestt = lt;
  }

  //
  // update bounds
  //
  if (ldiff < 0.) tmax = lt;
  else tmin = lt;

  int i=0;
  
  if (debug > 0) std::cerr << "refine2 i " << i << " nn " << numnodes << " tmin " << tmin << " tmax " << tmax << " err " << err << " lt " << lt << " llt " << llt << " ldiff " << ldiff << " lldiff " << lldiff << std::endl;

  //
  // secant/bisection search for tolerance
  // master does the logic, others only run setnodelocs
  //
  while (1) {
    int done = 0;
    if (mpimaster) {
      done = 1;
      if (fabs(ldiff) > 1.e-6
	  && fabs(lt-llt) > 1.e-14*fabs(lt+llt) ) done = 0;
    }
    MPI_Bcast(&done,1,MPI_INT,0,MPI_COMM_WORLD);
    if (done == 1) break;

    i++;
    //
    // make new estimate and bisect if out of bounds
    //
    double t = 1.;
    if (mpimaster) {
      t = lt-ldiff*(lt-llt)/(ldiff-lldiff);
      if (t >= tmax || t <= tmin || (ldiff-lldiff) == 0.) t = (tmax+tmin)/2.;
      if (debug > 0) std::cerr << " t " << t << std::endl;
    }
    MPI_Bcast(&t,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

    //
    // compute new error
    //
    err = setnodelocs(bms,b,t,numnodes);

    if (mpimaster) {
      //
      // update bounds
      //
      if (err < t) tmax = t;
      else tmin = t;

      //
      // update saved results for secant search
      //
      lldiff = ldiff;
      ldiff = err-t;
      llt = lt;
      lt = t;

      //
      // update best value
      //
      if (fabs(ldiff) < fabs(bestld)) {
	bestld = ldiff;
	bestt = lt;
      }

      if (debug > 0) std::cerr << "refine2 i " << i << " nn " << numnodes << " tmin " << tmin << " tmax " << tmax << " err " << err << " lt " << lt << " llt " << llt << " ldiff " << ldiff << " lldiff " << lldiff << std::endl;
    }
  }

  //
  // move to best result in case secant search was noisy
  //
  MPI_Bcast(&bestt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  err = setnodelocs(bms,b,bestt,numnodes);

  if (mpimaster) {
    if (debug > 0) std::cerr << "UTriTables::refine_boundaries best diff " << bestld << " t " << bestt << " err " << err << std::endl;

    if (err > 1.0) std::cerr << "Warning: UTriTables::refine_boundaries failed to meet desired tolerance" << std::endl;

    if (debug > 0) std::cerr << "UTriTables::refine_boundaries numerrcalcs " << numerrcalcs << std::endl;
  }

  //
  // stop evaluation threads
  //
  if (totthr > 0)
    for (std::size_t i=0;i<bms.size();i++) bms[i][b]->stop_eval_threads();

  return std::max(err,bestt);

}

void UTriTables::mesh_tables_PX(const double tolerance)
{
  //
  // overestimate desired tolerance to hedge for inexact error calculations
  //
  double boundary_tolerance = tolerance*0.9;
  double region_tolerance   = tolerance*0.99;
  std::vector<int> globalvars;
  globalvars.push_back(0);
  globalvars.push_back(1);
  globalvars.push_back(2);
  globalvars.push_back(3);
  globalvars.push_back(4);
  globalvars.push_back(5);
  globalvars.push_back(7);
  globalvars.push_back(8);
  globalvars.push_back(11);
  ErrorMap boundary_EMap(boundary_tolerance,1.e-2,globalvars);
  //
  // hack -- find the critical point node for use in tolerance definition
  //
  int cpn = -1;
  for (std::size_t i=0;i<rtmeshes[0].nodes.size();i++)
    if (rtmeshes[0].nodes[i].d.outputs[0] == 322.0 &&
	rtmeshes[0].nodes[i].d.outputs[1] == 647.096) cpn = i;
  std::vector<int> critvars;
  critvars.push_back(0);
  critvars.push_back(1);
  critvars.push_back(2);
  critvars.push_back(3);
  critvars.push_back(4);
  critvars.push_back(5);
  if (cpn >= 0) {
    std::vector<double> & d(rtmeshes[0].nodes[cpn].d.inputs);
    boundary_EMap.add_error_region(d[0],d[1],d[0]*1.e-3,d[1]*1.e-3,boundary_tolerance,1.e-2,critvars);
  }

  //
  // refine boundaries
  //
  std::vector<double> btimes(rtmeshes[0].boundaries.size());
  if (mpimaster) std::cout << "Refining boundaries of PX mesh" << std::endl;
  for (std::size_t i=0;i<rtmeshes[0].boundaries.size();i++) {
    double starttime = MPI_Wtime();
    //
    // create boundaries for each table
    //
    std::vector<int> bl;
    bl.push_back(i);
    for (std::size_t j=0;j<rtmeshes.size();j++) {
      bmst[j].push_back(new BMesh(&rtmeshes[j],bl,meshtype,logvars,totthr,boundary_EMap));
    }

    //
    // refine all boundaries concurrently
    //
    double accuracy = refine_boundaries(bmst,false,-1,0);
    double endtime = MPI_Wtime();
    btimes[i] = endtime-starttime;
    if (mpimaster) std::cout << "boundary " << i << " nodes " << bmst[0][i]->get_num_nodes() << " accuracy " << accuracy*boundary_tolerance << " time " << btimes[i] << std::endl;
  }
  double btime(0.0);
  for (std::size_t i=0;i<btimes.size();i++) btime += btimes[i];
  std::cout << "bmesh time " << btime << std::endl;

  ErrorMap region_EMap(region_tolerance,1.e-2,globalvars);
  if (cpn >= 0) {
    std::vector<double> & d(rtmeshes[0].nodes[cpn].d.inputs);
    region_EMap.add_error_region(d[0],d[1],d[0]*1.3e-2,d[1]*1.3e-2,region_tolerance,1.e-2,critvars);
  }

  //
  // refine regions
  //
  std::vector<double> rtimes(rtmeshes[0].regions.size());
  if (mpimaster) std::cout << "Refining and smoothing PX mesh regions" << std::endl;
  for (std::size_t i=0;i<rtmeshes[0].regions.size();i++) {
    double starttime = MPI_Wtime();
    //
    // create regions for each table
    //
    for (std::size_t j=0;j<rtmeshes.size();j++) {
      rmst[j].push_back(new RMesh(&rtmeshes[j],i,bmst[j],meshtype,logvars,totthr,region_EMap,0));
    }

    //
    // refine all regions concurrently
    //
    double accuracy;
    if (nodealg == "optimize") accuracy = refine_regions_pt(rmst,-1);
    else accuracy = refine_regions(rmst,-1);
    double endtime = MPI_Wtime();
    rtimes[i] = endtime-starttime;
    if (mpimaster) std::cout << "region " << i
			     << " total nodes = " << rmst[0][i]->get_num_nodes()
			     << " boundary nodes = " << rmst[0][i]->get_num_boundary_nodes()
			     << " accuracy = " << accuracy*region_tolerance
			     << " time = " << rtimes[i] << std::endl;

  }
  double rtime(0.0);
  for (std::size_t i=0;i<rtimes.size();i++) rtime += rtimes[i];
  std::cout << "rmesh time " << rtime << std::endl;

  std::cout << "btimes:";
  for (std::size_t i=0;i<btimes.size();i++) std::cout << " " << i << " = " << btimes[i];
  std::cout << std::endl;
  std::cout << "rtimes:";
  for (std::size_t i=0;i<rtimes.size();i++) std::cout << " " << i << " = " << rtimes[i];
  std::cout << std::endl;
  std::cout << "btime = " << btime << " rtime = " << rtime << " total time = " << btime+rtime << std::endl;
}

UTriMMOpt::UTriMMOpt( const bool master,
                      std::vector<RMesh *> & regionmeshes,
                      RMeshMMOptData & meshdata,
                      const int logscale,
                      int mpirank,
                      int mpisize,
                      int region_samples )
  : MCMin(master), rmeshes(0), data(meshdata), logvars(logscale),
    myrank(mpirank), mpicommsize(mpisize), rsamples(region_samples), debug(1)
{
  //
  // generate sub-meshes for each region mesh
  //
  for (std::size_t i=0;i<regionmeshes.size();i++) rmeshes.push_back(new RMesh(regionmeshes[i],meshdata));

  //
  // generate boundary point mapping
  //
  nbp = data.boundary_nodes.size();
  for (int i=0;i<nbp;i++) boundarymap[data.boundary_nodes[i]] = i;

}

UTriMMOpt::~UTriMMOpt()
{
  for (std::size_t i=0;i<rmeshes.size();i++) delete rmeshes[i];
}

void UTriMMOpt::getInitState( std::vector<double> & state )
{
  if (mpimaster) {
    rmeshes[0]->getSubMeshState(data.opt_nodes.size(),state);
  }
  else throw std::runtime_error("UTriMMOpt::getInitState: not master process");

}

void UTriMMOpt::updateWidths( std::vector<double> & widths )
{
  for (std::size_t i=0;i<widths.size();i+=2) {
    widths[i] = data.xywidth[0];
    widths[i+1] = data.xywidth[1];
  }
}

void UTriMMOpt::printMeshState( const std::vector<double> & state,
                                std::ostream & ofile )
{
  if (mpimaster) {
    //
    // triangulate with passed node locations
    //
    rmeshes[0]->modify_submesh(data,state);
    
    //
    // print the boundary and the mesh
    //
    rmeshes[0]->print_boundary(ofile);
    rmeshes[0]->print_tris(ofile);
  }

}

double UTriMMOpt::f( const std::vector<double> & state )
{
  std::vector<double> times;

  // start time
  times.push_back(MPI_Wtime());

  int fail = 0;
  if (mpimaster) {
    //
    // triangulate with new node locations -- false return value means
    // Delaunay mesh not preserved
    //
    try {
      rmeshes[0]->modify_submesh(data,state);
    }
    catch (std::runtime_error & e) {
      std::ostringstream oss;
      oss << "UTriMMOpt::f: modify_submesh failed on master mesh with error: " << e.what() << " Marking mesh as invalid" << std::endl;
      fail = 1;
    }
  }
  MPI_Bcast(&fail,1,MPI_INT,0,MPI_COMM_WORLD);
  if (fail != 0) return std::numeric_limits<double>::max();

  // master triangulation and delaunay check time
  times.push_back(MPI_Wtime());

  //
  // get and share triangles to recompute
  //
  std::vector<int> comptris;
  if (mpimaster) {
    std::vector<int> tnodes(state.size()/2);
    int totn = rmeshes[0]->get_num_nodes();
    for (std::size_t i=0;i<tnodes.size();i++) tnodes[i] = totn-1-i;
    rmeshes[0]->find_unsaved_errors(tnodes,comptris);
  }
  int cnsize = comptris.size();
  MPI_Bcast(&cnsize,1,MPI_INT,0,MPI_COMM_WORLD);
  if (!mpimaster) comptris.resize(cnsize);
  MPI_Bcast(comptris.data(),cnsize,MPI_INT,0,MPI_COMM_WORLD);

  // mesh transfer time
  times.push_back(MPI_Wtime());

  //
  // compute errors on triangles
  //
  double merr = -1.;
  for (std::size_t i=0;i<rmeshes.size();i++) {
    int etri;
    std::vector<double> eloc;
    double err = rmeshes[i]->compute_mesh_error_thread_saved(comptris,etri,eloc,rsamples,1);
    if (err > merr) merr = err;
  }

  //
  // maximize error on master
  //
  if (mpimaster) MPI_Reduce(MPI_IN_PLACE,&merr,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  else MPI_Reduce(&merr,NULL,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Bcast(&merr,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

  return merr;
}

bool UTriMMOpt::checkState( const std::vector<double> & state )
{
  bool goodstate = true;

  //
  // check if each point in the state falls within one of the
  // reference triangles
  //
  std::vector<double> b(3);
  for (std::size_t j=0;j<state.size();j+=2) {
    int mytri = -1;
    for (std::size_t i=0;i<data.domtris.size();i++) {
      data.domtris[i].getBcoords(state[j],state[j+1],b);
      if (b[0] >= 0. && b[0] <= (1.+2.e-14) &&
          b[1] >= 0. && b[1] <= (1.+2.e-14) &&
          b[2] >= 0. && b[2] <= (1.+2.e-14)) {
        mytri = i;
        break;
      }
    }
    if (mytri < 0) {
      goodstate = false;
      break;
    }
    else {
      // reject state if the point is too close to the boundary
      std::vector<int> p(3,-1);

      for (int i=0;i<3;i++)
        if (boundarymap.count(data.domtris[mytri].n[i]) > 0)
          p[i] = boundarymap[data.domtris[mytri].n[i]];

      for (int i=0;i<3;i++) {
        if (p[i] >= 0 && p[(i+1)%3] >= 0 && (p[i]+1)%nbp == p[(i+1)%3]
            && rmeshes[0]->get_angle_sq(p[i],p[(i+1)%3],state[j],state[j+1]) < 1.e-8) {
          goodstate = false;
        }
        else if (p[i] >= 0 &&
                 (rmeshes[0]->get_angle_sq(p[i],(p[i]+1)%nbp,state[j],state[j+1]) < 1.e-8 ||
                  rmeshes[0]->get_angle_sq(p[i],(p[i]-1+nbp)%nbp,state[j],state[j+1]) < 1.e-8)) {
          goodstate = false;          
        }
      }

      if (!goodstate) break;
    }
  }

  return goodstate;
}
