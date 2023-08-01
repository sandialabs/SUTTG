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
#include <string>
#include <sstream>
#include <algorithm>
#include <limits>
#include <iomanip>

//
// EOS library defines
//
#include "EOSParam.H"
#include "EOSModel.H"
#include "EOSFactory.H"
#include "PhaseBoundaryInfo.H"

//
// build table defines
//
#include "Tables.H"
#include "TriMesh.H"
#include "RMesh.H"
#include "BMesh.H"
#include "mpi.h"

//
// table output defines
//
#include "utri_eos_support.h"
#include "utri_eos_mig.h"

Tables::Tables( std::vector<Parameters> & model_params,
		const std::string model_name,
		const std::string mesh_type,
		const std::vector<double> & xvar_bounds,
		const std::vector<double> & yvar_bounds,
		const int logxy,
		const int totalthreads,
		const int mpirank,
		const int mpisize )
  : Xbounds(xvar_bounds), Ybounds(yvar_bounds),
    mintab(0), maxtab(model_params.size()),
    mpimaster(false), myrank(mpirank), mpicommsize(mpisize),
    logvars(logxy), totthr(totalthreads)
{
  //
  // setup the master process
  //
  if (mpirank == 0) mpimaster = true;

  //
  // setup the range of tables
  //
  int psize = model_params.size();
  if (mpisize < 1) throw std::runtime_error("Tables::Tables: mpisize must be at least 1");
  else if (mpisize > psize) throw std::runtime_error("Tables::Tables: mpisize must be less than number of tables");
  int ntabs = psize/mpisize;
  int extra = psize%mpisize;
  mintab = (ntabs+1)*std::min(mpirank,extra)+ntabs*std::max(0,mpirank-extra);
  if (mpirank < extra) ntabs++;
  maxtab = std::min(mintab+ntabs,maxtab);

  //
  // check table bounds -- two and increasing
  //
  if (Xbounds.size() != 2) throw std::runtime_error("Tables::Tables: Xbounds must be of size 2");
  if (Ybounds.size() != 2) throw std::runtime_error("Tables::Tables: Ybounds must be of size 2");
  if (Xbounds[0] > Xbounds[1]) throw std::runtime_error("Tables::Tables: Xbounds must be in increasing order");
  if (Ybounds[0] > Ybounds[1]) throw std::runtime_error("Tables::Tables: Ybounds must be in increasing order");

  //
  // check mesh type
  //
  if (mesh_type == "PH") meshtype = TriMesh::PH_MESH;
  else if (mesh_type == "PS") meshtype = TriMesh::PS_MESH;
  else if (mesh_type == "PE") meshtype = TriMesh::PE_MESH;
  else throw std::runtime_error("Tables::Tables: invalid mesh type "+mesh_type);

  //
  // create an EOS factory
  //
  EOSFactory eosfactory;

  //
  // create the meshes
  //
  for (int i=mintab;i<maxtab;i++) {
    std::cout << "# Creating model " << i << std::endl;

    //
    // create the next model
    //
    EOSModel * model = eosfactory.getModel(model_params[i],model_name);

    //
    // calculate phase boundary information
    //
    model->calculatePhaseBoundaryInfo(Xbounds[0],Xbounds[1]);

    //
    // get the phase boundary information
    //
    std::vector<PhaseBoundaryInfo *> * pbi = model->getPhaseBoundaryInfo();

    //
    // set the minimum/maximum temperature
    //
    model->setTemperatureBounds(Ybounds[0]*0.5,Ybounds[1]*1.1);

    //
    // create the meshes
    //
    rtmeshes.push_back(TriMesh(pbi,model));
  }

  create_regions_PX();
  rmst.resize(maxtab-mintab);

}

Tables::~Tables()
{
  //
  // cleanup the models
  //
  for (std::size_t i=0;i<rtmeshes.size();i++) delete rtmeshes[i].model;

  //
  // cleanup any meshes
  //
  for (std::size_t i=0;i<rmst.size();i++)
    for (std::size_t j=0;j<rmst[i].size();j++)
      if (rmst[i][j] != 0) delete rmst[i][j];
  for (std::size_t i=0;i<rmse.size();i++)
    for (std::size_t j=0;j<rmse[i].size();j++)
      if (rmse[i][j] != 0) delete rmse[i][j];

}

void Tables::create_regions_PX()
{
  for (std::size_t i=0;i<rtmeshes.size();i++) {
    //
    // References to the current mesh
    //
    TriMesh & rtmesh = rtmeshes[i];

    //
    // Make the rough P-X mesh regions
    //

    //
    // hardcoding the iapwsif97 and iapws95 models for now
    //
    // require the vapor dome to intersect at most the constant P boundaries
    //
    
    //
    // Make the corner nodes
    //
    Node n(2,13);
    n.d.inputs[0] = Xbounds[0];
    n.d.inputs[1] = Ybounds[0];
    rtmesh.model->getThermalEOS_PT_D(n.d,3);
    rtmesh.nodes.push_back(n);
    n.d.inputs[0] = Xbounds[1];
    n.d.inputs[1] = Ybounds[0];
    rtmesh.model->getThermalEOS_PT_D(n.d,3);
    rtmesh.nodes.push_back(n);
    n.d.inputs[0] = Xbounds[1];
    n.d.inputs[1] = Ybounds[1];
    rtmesh.model->getThermalEOS_PT_D(n.d,3);
    rtmesh.nodes.push_back(n);
    n.d.inputs[0] = Xbounds[0];
    n.d.inputs[1] = Ybounds[1];
    rtmesh.model->getThermalEOS_PT_D(n.d,3);
    rtmesh.nodes.push_back(n);

    //
    // swap the corner nodes to the PH space
    //
    rtmesh.nodes_swap_PTtoPX(TriMesh::PH_MESH);

    //
    // check for intersection of vapor dome with lower pressure boundary
    //
    if (Xbounds[0] < 22064000.) {
      double lt,lh[2];
      rtmesh.model->getPhaseCoexistence(4,Xbounds[0],lt,lh);
      bool lowerboundary = false;
      if (lt >= Ybounds[0] && lt <= Ybounds[1]) {
	lowerboundary = true;
	//
	// add nodes for phase boundary
	//
	n.d.inputs[0] = Xbounds[0];
	n.d.inputs[1] = lh[0];
	rtmesh.model->getThermalEOS_PH_D(n.d,3);
	rtmesh.nodes.push_back(n);
	n.d.inputs[1] = lh[1];
	rtmesh.model->getThermalEOS_PH_D(n.d,3);
	rtmesh.nodes.push_back(n);
	n.d.inputs[1] = lh[0];
	rtmesh.model->getThermalEOS_PH_D(n.d,4);
	rtmesh.nodes.push_back(n);
	n.d.inputs[1] = lh[1];
	rtmesh.model->getThermalEOS_PH_D(n.d,4);
	rtmesh.nodes.push_back(n);
      }

      bool critinside = false;
      bool upperboundary = false;
      if (Xbounds[1] >= 22064000.) {
	//
	// check if critical point is in boundaries
	//
	if (Ybounds[0] < 647.096 && Ybounds[1] > 647.096
	    && Xbounds[0] < 22064000. && Xbounds[1] > 22064000.) {
	  critinside = true;
	  //
	  // add nodes for critical point
	  //
	  // pressure-temperature version does not work well when
	  // models are noisy due to very flat isotherms
	  // so use density-temperature version since that will be exact
	  //
	  n.d.inputs[0] = 322.0;
	  n.d.inputs[1] = 647.096;
	  rtmesh.model->getThermalEOS_RT_D(n.d,3);
	  std::swap(n.d.inputs[0],n.d.outputs[0]);
	  std::swap(n.d.inputs[1],n.d.outputs[1]);
	  rtmesh.nodes.push_back(n);
	  rtmesh.model->getThermalEOS_PH_D(n.d,4);
	  rtmesh.nodes.push_back(n);
	}
      }
      else {
	//
	// check if boundary in region
	//
	rtmesh.model->getPhaseCoexistence(4,Xbounds[1],lt,lh);
	if (lt >= Ybounds[0] && lt <= Ybounds[1]) {
	  upperboundary = true;
	  //
	  // add nodes for phase boundary
	  //
	  n.d.inputs[0] = Xbounds[1];
	  n.d.inputs[1] = lh[0];
	  rtmesh.model->getThermalEOS_PH_D(n.d,3);
	  rtmesh.nodes.push_back(n);
	  n.d.inputs[1] = lh[1];
	  rtmesh.model->getThermalEOS_PH_D(n.d,3);
	  rtmesh.nodes.push_back(n);
	  n.d.inputs[1] = lh[0];
	  rtmesh.model->getThermalEOS_PH_D(n.d,4);
	  rtmesh.nodes.push_back(n);
	  n.d.inputs[1] = lh[1];
	  rtmesh.model->getThermalEOS_PH_D(n.d,4);
	  rtmesh.nodes.push_back(n);
	}
      }

      if ((lowerboundary && !(critinside || upperboundary))
	  || (!lowerboundary && (critinside || upperboundary))) {
	throw std::runtime_error("Tables::create_regions_PX: liquid vapor boundary crosses fixed T boundary");
      }
    }

    //
    // swap all nodes to the proper space
    //
    rtmesh.nodes_swap_PHtoPX(meshtype);

    //
    // Create boundaries/regions based upon number of nodes
    //
    int nnodes = rtmesh.nodes.size();
    if (nnodes == 4) {
      //
      // Make the boundaries
      //
      std::vector<int> en1(1),en2(1);

      int curphase = 3;
      // low T
      en1[0] = 0;
      en2[0] = 1;
      rtmesh.boundaries.push_back(Boundary(1,en1,en2,Ybounds[0],curphase));
      // high T
      en1[0] = 3;
      en2[0] = 2;
      rtmesh.boundaries.push_back(Boundary(1,en1,en2,Ybounds[1],curphase));
      // low P
      en1[0] = 0;
      en2[0] = 3;
      rtmesh.boundaries.push_back(Boundary(3,en1,en2,Xbounds[0],curphase));
      // high P
      en1[0] = 1;
      en2[0] = 2;
      rtmesh.boundaries.push_back(Boundary(3,en1,en2,Xbounds[1],curphase));

      //
      // Make the region
      //
      Region r;
      r.bounds.push_back(0);
      r.bounds.push_back(3);
      r.bounds.push_back(1);
      r.bounds.push_back(2);
      r.phase = curphase;
      rtmesh.regions.push_back(r);
    }
    else if (nnodes == 10) {
      //
      // Make the boundaries
      //
      std::vector<int> en1(1),en2(1);

      int phase1 = 3;
      int phase2 = 4;
      // low T
      en1[0] = 0;
      en2[0] = 1;
      rtmesh.boundaries.push_back(Boundary(1,en1,en2,Ybounds[0],phase1));
      // high T
      en1[0] = 3;
      en2[0] = 2;
      rtmesh.boundaries.push_back(Boundary(1,en1,en2,Ybounds[1],phase1));
      // low P 1
      en1[0] = 0;
      en2[0] = 5;
      rtmesh.boundaries.push_back(Boundary(3,en1,en2,Xbounds[0],phase1));
      // low P 2
      en1[0] = 4;
      en2[0] = 3;
      rtmesh.boundaries.push_back(Boundary(3,en1,en2,Xbounds[0],phase1));
      // low P 3
      en1[0] = 7;
      en2[0] = 6;
      rtmesh.boundaries.push_back(Boundary(3,en1,en2,Xbounds[0],phase2));
      // high P
      en1[0] = 1;
      en2[0] = 2;
      rtmesh.boundaries.push_back(Boundary(3,en1,en2,Xbounds[1],phase1));

      en1.resize(2);
      en2.resize(2);

      // low H l-v
      en1[0] = 5; en1[1] = 7;
      en2[0] = 8; en2[1] = 9;
      rtmesh.boundaries.push_back(Boundary(0,en1,en2,0,2));
      // high H l-v
      en1[0] = 4; en1[1] = 6;
      en2[0] = 8; en2[1] = 9;
      rtmesh.boundaries.push_back(Boundary(0,en1,en2,0,0));

      //
      // Make the regions
      //
      Region r;
      r.bounds.push_back(0);
      r.bounds.push_back(5);
      r.bounds.push_back(1);
      r.bounds.push_back(3);
      r.bounds.push_back(7);
      r.bounds.push_back(6);
      r.bounds.push_back(2);
      r.phase = phase1;
      rtmesh.regions.push_back(r);

      r.bounds.clear();
      r.bounds.push_back(6);
      r.bounds.push_back(7);
      r.bounds.push_back(4);
      r.phase = phase2;
      rtmesh.regions.push_back(r);
    }
    else if (nnodes == 12) {
      //
      // Make the boundaries
      //
      std::vector<int> en1(1),en2(1);

      int phase1 = 3;
      int phase2 = 4;
      // low T
      en1[0] = 0;
      en2[0] = 1;
      rtmesh.boundaries.push_back(Boundary(1,en1,en2,Ybounds[0],phase1));
      // high T
      en1[0] = 3;
      en2[0] = 2;
      rtmesh.boundaries.push_back(Boundary(1,en1,en2,Ybounds[1],phase1));
      // low P 1
      en1[0] = 0;
      en2[0] = 5;
      rtmesh.boundaries.push_back(Boundary(3,en1,en2,Xbounds[0],phase1));
      // low P 2
      en1[0] = 4;
      en2[0] = 3;
      rtmesh.boundaries.push_back(Boundary(3,en1,en2,Xbounds[0],phase1));
      // low P 3
      en1[0] = 7;
      en2[0] = 6;
      rtmesh.boundaries.push_back(Boundary(3,en1,en2,Xbounds[0],phase2));
      // high P 1
      en1[0] = 1;
      en2[0] = 9;
      rtmesh.boundaries.push_back(Boundary(3,en1,en2,Xbounds[1],phase1));
      // high P 2
      en1[0] = 8;
      en2[0] = 2;
      rtmesh.boundaries.push_back(Boundary(3,en1,en2,Xbounds[1],phase1));
      // high P 3
      en1[0] = 11;
      en2[0] = 10;
      rtmesh.boundaries.push_back(Boundary(3,en1,en2,Xbounds[1],phase1));

      en1.resize(2);
      en2.resize(2);

      // low H l-v
      en1[0] = 5; en1[1] = 7;
      en2[0] = 9; en2[1] = 11;
      rtmesh.boundaries.push_back(Boundary(0,en1,en2,0,2));
      // high H l-v
      en1[0] = 4; en1[1] = 6;
      en2[0] = 8; en2[1] = 10;
      rtmesh.boundaries.push_back(Boundary(0,en1,en2,0,0));

      //
      // Make the regions
      //
      Region r;
      r.bounds.push_back(0);
      r.bounds.push_back(5);
      r.bounds.push_back(8);
      r.bounds.push_back(2);
      r.phase = phase1;
      rtmesh.regions.push_back(r);

      r.bounds.clear();
      r.bounds.push_back(6);
      r.bounds.push_back(1);
      r.bounds.push_back(3);
      r.bounds.push_back(9);
      r.phase = phase1;
      rtmesh.regions.push_back(r);

      r.bounds.clear();
      r.bounds.push_back(7);
      r.bounds.push_back(9);
      r.bounds.push_back(4);
      r.bounds.push_back(8);
      r.phase = phase2;
      rtmesh.regions.push_back(r);
    }
    else {
      throw std::runtime_error("Tables::create_regions_PX: invalid number of nodes");
    }
  }
}

void Tables::write_table( const std::string filename,
                          std::vector<TMesh *> & rms,
                          const TriMesh::MeshTypes mtype,
                          const std::vector<double> & defaults )
{
  //
  // fill the table structure
  //
  eos_table_utri_t table;

  //
  // set defaults
  //
  if ((table.defaults = (double *)malloc(sizeof(double)*5)) == NULL) {
    throw std::runtime_error("failed to allocate defaults memory in Tables::write_table");
  }
  for (int j=0;j<5;j++) table.defaults[j] = defaults[j];

  //
  // set flags
  //
  if (logvars == 0) table.flags[0] = 0;
  else table.flags[0] = 1;
  table.flags[1] = I_RT;
  if (mtype == TriMesh::RE_MESH) table.flags[1] = I_RE;
  else if (mtype == TriMesh::PH_MESH) table.flags[1] = I_PH;
  else if (mtype == TriMesh::PS_MESH) table.flags[1] = I_PS;
  else if (mtype == TriMesh::PE_MESH) table.flags[1] = I_PE;
    
  //
  // each mesh is in its own table file -> 1 mode
  //
  table.nmodes = 1;
  if ((table.modes = (eos_data_points_t **)malloc(sizeof(eos_data_points_t *))) == NULL) {
    throw std::runtime_error("failed to malloc memory for modes in Tables::write_table");
  }
  table.modes[0] = NULL;
  table.tris = NULL;
  table.np = 0;
  table.nt = 0;

  //
  // holder for boundary nodes
  //
  std::deque<Node> bnodes[2];

  //
  // add each region to table
  //
  for (std::size_t j=0;j<rms.size();j++) {
    TMesh & r = *rms[j];

    //
    // get the nodes
    //
    std::vector<Node> rnodes;
    r.get_nodes(rnodes);

    //
    // get the triangles
    //
    std::vector<Tri> rtris;
    r.get_tris(rtris);

    //
    // resize data holder
    //
    if ((table.modes[0] = (eos_data_points_t *) eospointsrealloc(table.modes[0],table.np+rnodes.size())) == NULL) {
      throw std::runtime_error("failed to reallocate node storage in Tables::write_table");
    }

    //
    // add node data
    //
    for (std::size_t k=0;k<rnodes.size();k++) {
      set_eosdatapoint(table.modes[0],table.np+k,mtype,rnodes[k]);
    }

    //
    // resize tri holder
    //
    if ((table.tris = (eos_table_tri_t *) realloc(table.tris,sizeof(eos_table_tri_t)*(table.nt+rtris.size()))) == NULL) {
      throw std::runtime_error("failed to allocate tri storage in Tables::write_table");
    }

    //
    // add triangle data with appropriate offset
    //
    for (std::size_t k=0;k<rtris.size();k++) {
      table.tris[table.nt+k].p[0] = rtris[k].n[0]+table.np;
      table.tris[table.nt+k].p[1] = rtris[k].n[1]+table.np;
      table.tris[table.nt+k].p[2] = rtris[k].n[2]+table.np;
    }

    //
    // update table totals
    //
    table.nt += rtris.size();
    table.np += rnodes.size();

    //
    // record boundary points if RE or PX grid
    //
    if (mtype == TriMesh::RE_MESH ||
	mtype == TriMesh::PH_MESH ||
	mtype == TriMesh::PS_MESH ||
	mtype == TriMesh::PE_MESH ) {

      bool last1(false),last2(false);
      for (int k=0;k<r.get_num_boundary_nodes();k++) {
        if (rnodes[k].d.outputs[1] == Ybounds[0]) {
          if (k == 0 && rnodes[1].d.outputs[1] != Ybounds[0]) last1 = true;
          // lower boundary goes on back
          else bnodes[0].push_back(rnodes[k]);
        }
        else if (rnodes[k].d.outputs[1] == Ybounds[1]) {
          if (k == 0 && rnodes[1].d.outputs[1] != Ybounds[1]) last2 = true;
          // upper boundary goes on front due to anti-clockwise ordering
          else bnodes[1].push_front(rnodes[k]);
        }
      }
      //
      // if this boundary started with only one point on the edge,
      // stick that point on now
      //
      if (last1) bnodes[0].push_back(rnodes[0]);
      if (last2) bnodes[1].push_front(rnodes[0]);
    }
  }

  //
  // add boundary points
  //

  //
  // resize data holder
  //
  if ((table.modes[0] = (eos_data_points_t *) eospointsrealloc(table.modes[0],table.np+2*(bnodes[0].size()+bnodes[1].size()))) == NULL) {
    throw std::runtime_error("failed to reallocate node storage in Tables::write_table");
  }

  //
  // add node data -- twice for each boundary
  //
  for (int j=0;j<2;j++) {
    table.nbp[j] = bnodes[j].size();
    table.nbi[j] = table.np;
    table.nbt[j] = 0;
    for (std::size_t k=0;k<bnodes[j].size();k++) {
      set_eosdatapoint(table.modes[0],table.np+k,mtype,bnodes[j][k]);
      set_eosdatapoint(table.modes[0],table.np+k+table.nbp[j],mtype,bnodes[j][k]);
      //
      // flag phase as extrapolation
      //
      table.modes[0]->phase[table.np+k] *= -1;
      table.modes[0]->phase[table.np+k+table.nbp[j]] *= -1;
      //
      // count need for new tris
      //
      if (k > 0) {
	if (table.modes[0]->phase[table.np+k] == table.modes[0]->phase[table.np+k-1])
	  table.nbt[j] += 2;
      }
    }
    table.np += 2*bnodes[j].size();
  }

  //
  // resize tri holder
  //
  int newt = table.nbt[0]+table.nbt[1];
  if ((table.tris = (eos_table_tri_t *) realloc(table.tris,sizeof(eos_table_tri_t)*(table.nt+newt))) == NULL) {
    throw std::runtime_error("failed to allocate tri storage in Tables::write_table");
  }

  //
  // initialize new triangle data
  //
  for (int k=0;k<newt;k++) {
    table.tris[table.nt+k].p[0] = -1;
    table.tris[table.nt+k].p[1] = -1;
    table.tris[table.nt+k].p[2] = -1;
  }

  //
  // create triangle offsets
  //
  table.nbt[1] = table.nt + table.nbt[0];
  table.nbt[0] = table.nt;

  table.nt += newt;

  //
  // initialize rptree pointers although the tree won't be created/written
  //
  table.rptree.splits = NULL;
  table.rptree.trilookup = NULL;
  table.rptree.trilist = NULL;

  std::cerr << "Saving table " << filename << std::endl;
  if (eostableutriwritenc(filename.c_str(),&table) != 0) {
    throw std::runtime_error("Failed to write utri EOS table");
  }
    
  //
  // cleanup
  //
  eos_table_utri_delete(&table);

}

void Tables::write_tables( const std::string basename,
                           const std::vector<double> & defaults )
{
  if (meshtype == TriMesh::PH_MESH || meshtype == TriMesh::PS_MESH || meshtype == TriMesh::PE_MESH) {
    //
    // PX tables
    //
    for (std::size_t i=0;i<rmst.size();i++) {
      //
      // make filename for this table
      //
      std::ostringstream filename;
      filename << basename;
      if (meshtype == TriMesh::PH_MESH) filename << "-ph.nc";
      else if (meshtype == TriMesh::PH_MESH) filename << "-ps.nc";
      else filename << "-pe.nc";
      
      std::cerr << "Building table " << filename.str() << std::endl;
      
      write_table(filename.str(),rmst[i],meshtype,defaults);
    }
  }
}

void Tables::set_eosdatapoint( eos_data_points_t *d,
                               int p,
                               TriMesh::MeshTypes mtype,
                               Node & n )
{

  d->phase[p] = int(n.d.outputs[4]);
  d->rho[p]   = n.d.inputs[0];
  d->p[p]     = n.d.outputs[0];
  if (mtype == TriMesh::RT_MESH) {
    d->T[p]   = n.d.inputs[1];
    d->E[p]   = n.d.outputs[1];
  }
  else {
    d->T[p]   = n.d.outputs[1];
    d->E[p]   = n.d.inputs[1];
  }
  d->S[p]     = n.d.outputs[2];
  d->F[p]     = n.d.outputs[3];

  //
  // copy derivative values if available
  //
  if (n.d.outputs.size() > 9) {
    d->G[p]     = 0.;
    d->H[p]     = 0.;
    d->dpdr[p]  = n.d.outputs[5];
    d->dpdT[p]  = n.d.outputs[6];
    d->dEdr[p]  = n.d.outputs[7];
    d->dEdT[p]  = n.d.outputs[8];
    d->cs[p]    = n.d.outputs[9];
    d->dcsdr[p] = n.d.outputs[10];
  }
  else {
    d->G[p]     = 0.;
    d->H[p]     = 0.;
    d->dpdr[p]  = 0.;
    d->dpdT[p]  = 0.;
    d->dEdr[p]  = 0.;
    d->dEdT[p]  = 0.;
    d->cs[p]    = 0.;
    d->dcsdr[p] = 0.;
  }

  if (mtype == TriMesh::PH_MESH || mtype == TriMesh::PS_MESH || mtype == TriMesh::PE_MESH) {
    d->phase[p] = int(n.d.outputs[10]);
    d->p[p] = n.d.inputs[0];
    d->rho[p] = n.d.outputs[0];
    d->T[p] = n.d.outputs[1];
    d->G[p] = 0.;
    d->H[p] = 0.;
    d->dEdr[p] = n.d.outputs[4]; // really G
    d->F[p] = n.d.outputs[5];
    d->dpdT[p] = n.d.outputs[6]; // really C_P
    d->dEdT[p] = n.d.outputs[7];
    d->cs[p] = n.d.outputs[8];
    d->dpdr[p] = n.d.outputs[9]; // really K_T
    if (mtype == TriMesh::PH_MESH) {
      d->dcsdr[p] = n.d.inputs[1]; // really H
      d->E[p] = n.d.outputs[2];
      d->S[p] = n.d.outputs[3];
    }
    else if (mtype == TriMesh::PS_MESH) {
      d->S[p] = n.d.inputs[1];
      d->E[p] = n.d.outputs[2];
      d->dcsdr[p] = n.d.outputs[3]; // really H
    }
    else { // TriMesh::PE_MESH
      d->E[p] = n.d.inputs[1];
      d->dcsdr[p] = n.d.outputs[2]; // really H
      d->S[p] = n.d.outputs[3];
    }
  }

  //
  // transport properties
  //
  if (n.d.outputs.size() > 10) {
    d->visc[p] = n.d.outputs[11];
    d->thcon[p] = n.d.outputs[12];
  }
  else {
    d->visc[p] = 0.;
    d->thcon[p] = 0.;
  }    
}
