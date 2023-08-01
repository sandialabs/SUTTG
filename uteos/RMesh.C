/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SUTTG, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <stack>
#include <limits>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "RMesh.H"
#include "BoundDelaunayTriangulation.H"
#include "qd/dd_real.h"

//
// Constructor
//
RMesh::RMesh( TriMesh * trimesh,
              int region,
              std::vector<BMesh *> & boundarymeshes,
              TriMesh::MeshTypes meshtype,
              int logscale,
              int tot_threads,
              const ErrorMap & error_map,
              int debug_level ) :
  TMesh(trimesh,region,meshtype,logscale,tot_threads,error_map,debug_level), a_n(0,0)
{
  //
  // Put in boundary nodes in anti-clockwise order -- this assumes the
  // first boundary is ordered correctly
  //
  std::vector<int> & bounds = trimesh->regions[region].bounds;
  int phase = trimesh->regions[region].phase;
  boundarymeshes[bounds[0]]->get_nodes(phase,nodes);
  for (std::size_t i=1;i<bounds.size();i++) {
    std::vector<Node> bn;
    boundarymeshes[bounds[i]]->get_nodes(phase,bn);
    if (bn.front() == nodes.back()) {
      nodes.insert(nodes.end(),bn.begin()+1,bn.end());
    }
    else if (bn.back() == nodes.back()) {
      nodes.insert(nodes.end(),bn.rbegin()+1,bn.rend());
    }
    else {
      std::cerr << "RMesh::RMesh: boundary " << i << " nodes " << bn.size() << std::endl;
      std::cerr << "RMesh::RMesh: n1 x1 " << nodes.back().d.inputs[0] << " x2 " << nodes.back().d.inputs[1]
		<< " y1 " << nodes.back().d.outputs[0] << " y4 " << nodes.back().d.outputs[1]
		<< " y3 " << nodes.back().d.outputs[2] << " y5 " << nodes.back().d.outputs[3]
		<< " y5 " << nodes.back().d.outputs[4] << " y6 " << nodes.back().d.outputs[5]
		<< " y7 " << nodes.back().d.outputs[6] << " y8 " << nodes.back().d.outputs[7]
		<< " y9 " << nodes.back().d.outputs[8] << " y10 " << nodes.back().d.outputs[9]
		<< " y11 " << nodes.back().d.outputs[10] << std::endl
		<< "              b1 x1 " << bn.front().d.inputs[0] << " x2 " << bn.front().d.inputs[1]
		<< " y1 " << bn.front().d.outputs[0] << " y4 " << bn.front().d.outputs[1]
		<< " y3 " << bn.front().d.outputs[2] << " y5 " << bn.front().d.outputs[3]
		<< " y5 " << bn.front().d.outputs[4] << " y6 " << bn.front().d.outputs[5]
		<< " y7 " << bn.front().d.outputs[6] << " y8 " << bn.front().d.outputs[7]
		<< " y9 " << bn.front().d.outputs[8] << " y10 " << bn.front().d.outputs[9]
		<< " y11 " << bn.front().d.outputs[10] << std::endl
		<< "              b1 x1 " << bn.back().d.inputs[0] << " x2 " << bn.back().d.inputs[1]
		<< " y1 " << bn.back().d.outputs[0] << " y4 " << bn.back().d.outputs[1]
		<< " y3 " << bn.back().d.outputs[2] << " y5 " << bn.back().d.outputs[3]
		<< " y5 " << bn.back().d.outputs[4] << " y6 " << bn.back().d.outputs[5]
		<< " y7 " << bn.back().d.outputs[6] << " y8 " << bn.back().d.outputs[7]
		<< " y9 " << bn.back().d.outputs[8] << " y10 " << bn.back().d.outputs[9]
		<< " y11 " << bn.back().d.outputs[10] << std::endl;

      throw std::runtime_error("boundary meshes do not have matching end nodes");
    }
  }

  //
  // Check that the boundary loop was closed and remove the redundant node
  //
  if (nodes.front() == nodes.back()) {
    nodes.pop_back();
  }
  else {
    throw std::runtime_error("boundary meshes do not form a loop");
  }

  nbnodes = nodes.size();

  //
  // determine bounds
  //
  xbounds[0] = 9.e99; xbounds[1] = 0.;
  ybounds[0] = 9.e99; ybounds[1] = 0.;
  for (std::size_t i=0;i<nodes.size();i++) {
    std::pair<double,double> xy = get_independent_vars(nodes[i]);
    if (xy.first > xbounds[1]) xbounds[1] = xy.first;
    if (xy.first < xbounds[0]) xbounds[0] = xy.first;
    if (xy.second > ybounds[1]) ybounds[1] = xy.second;
    if (xy.second < ybounds[0]) ybounds[0] = xy.second;
  }
  if (logvars) {
    xbounds[0] = log(xbounds[0]);
    xbounds[1] = log(xbounds[1]);
    ybounds[0] = log(ybounds[0]);
    ybounds[1] = log(ybounds[1]);
  }
}

//
// Constructor for submesh
//
RMesh::RMesh(RMesh * pmesh, const RMeshMMOptData & rmd)
  : TMesh(pmesh->tm,pmesh->tmregion,pmesh->type,pmesh->logvars,pmesh->totthr,pmesh->emap,0/*pmesh->debug*/), a_n(0,0)
{
  //
  // push on boundary nodes first
  //
  for (std::size_t i=0;i<rmd.boundary_nodes.size();i++)
    nodes.push_back(pmesh->nodes[rmd.boundary_nodes[i]]);
  nbnodes = nodes.size();

  //
  // put in other nodes
  //
  for (std::size_t i=0;i<rmd.interior_nodes.size();i++)
    nodes.push_back(pmesh->nodes[rmd.interior_nodes[i]]);
  for (std::size_t i=0;i<rmd.opt_nodes.size();i++) {
    if (rmd.opt_nodes[i] >= int(pmesh->nodes.size())) {
      std::cout << "RMesh::RMesh: add node" << std::endl;
      pmesh->add_node(rmd.opt_tri);
    }
    nodes.push_back(pmesh->nodes[rmd.opt_nodes[i]]);
  }
  
  //
  // copy bounds from parent
  //
  xbounds[0] = pmesh->xbounds[0];
  xbounds[1] = pmesh->xbounds[1];
  ybounds[0] = pmesh->ybounds[0];
  ybounds[1] = pmesh->ybounds[1];

  //
  // start thread evaluators
  //
  if (totthr > 0) start_eval_threads();

}

//
// Destructor
//
RMesh::~RMesh()
{

}

double RMesh::compute_submesh_error(const std::vector<int> & tlist,const int samples)
{
  double maxerr = 0.;
  for (std::size_t i=0;i<tlist.size();i++) {
    compute_errorb(i,samples);
    if (debug > 0) std::cerr << "t " << i << " err " << tris[i].maxerr << " loc " << tris[i].maxerrloc[0] << " " << tris[i].maxerrloc[1] << std::endl;
    if (tris[i].maxerr > maxerr) maxerr = tris[i].maxerr;
  }
  if (debug > 0) {
    std::cerr << "# t maxerr " << maxerr << std::endl << std::endl << std::endl;
  }

  return maxerr;
}

void RMesh::compute_mesh_errors(const int samples)
{
  if (ncache.size() == 0) initialize_error_points(samples);

  double maxerr = 0.;
  double minerr = 9.e99;
  compute_tris_error();
  for (std::size_t i=0;i<tris.size();i++) {
    if (debug > 0) std::cerr << "t " << i << " err " << tris[i].maxerr << " loc " << tris[i].maxerrloc[0] << " " << tris[i].maxerrloc[1] << std::endl;
    if (tris[i].maxerr > maxerr) {
      maxerr = tris[i].maxerr;
    }
    if (tris[i].maxerr < minerr) minerr = tris[i].maxerr;
  }
  if (debug > 0) {
    std::cerr << "# t minerr " << minerr << " maxerr " << maxerr << std::endl << std::endl << std::endl;
  }

}

double RMesh::compute_mesh_error_saved(const std::vector<int> & comptris,int & tri,std::vector<double> & loc,const int samples)
{
  loc.resize(2);
  double maxerr = 0.;
  double minerr = 9.e99;
  loc[0] = -9.e99;
  loc[1] = -9.e99;
  tri = -1;
  std::size_t cti = 0;
  for (std::size_t i=0;i<tris.size();i++) {
    if (cti < comptris.size() && comptris[cti] == int(i)) {
      //
      // recompute this tri
      //
      compute_errorb(i,samples);
      cti++;
    }
    else {
      //
      // get saved error values for this tri
      //
      TMeshTriData d(tris[i]);
      std::vector<TMeshTriData>::iterator tmtdi = std::lower_bound(tmtd.begin(),tmtd.end(),d);
      if (tmtdi != tmtd.end()) {
        if (tmtdi->n0 == tris[i].n[0] && tmtdi->n1 == tris[i].n[1] && tmtdi->n2 == tris[i].n[2]) {
          tris[i].maxerr = tmtdi->maxerr;
          tris[i].maxerrloc[0] = tmtdi->maxerrloc0;
          tris[i].maxerrloc[1] = tmtdi->maxerrloc1;
        }
        else throw std::runtime_error("RMesh::compute_mesh_error_saved: failed to find a tri in the saved vector");
      }
      else throw std::runtime_error("RMesh::compute_mesh_error_saved: failed to find a tri in the saved vector (off end)");
    }
    
    if (debug > 0) std::cerr << "t " << i << " err " << tris[i].maxerr << " loc " << tris[i].maxerrloc[0] << " " << tris[i].maxerrloc[1] << std::endl;
    if (tris[i].maxerr > maxerr) {
      tri = i;
      maxerr = tris[i].maxerr;
      loc[0] = tris[i].maxerrloc[0];
      loc[1] = tris[i].maxerrloc[1];
    }
    if (tris[i].maxerr < minerr) minerr = tris[i].maxerr;
  }
  if (debug > 1) {
    std::cerr << "# t minerr " << minerr << " maxerr " << maxerr << std::endl << std::endl << std::endl;
  }

  //
  // now save these new tris
  //
  tmtd.clear();
  for (std::size_t i=0;i<tris.size();i++) tmtd.push_back(TMeshTriData(tris[i]));

  return maxerr;
}

double RMesh::compute_mesh_error(int & tri,std::vector<double> & loc,const int samples)
{
  loc.resize(2);
  double maxerr = 0.;
  double minerr = 9.e99;
  loc[0] = -9.e99;
  loc[1] = -9.e99;
  tri = -1;
  for (std::size_t i=0;i<tris.size();i++) {
    // first check that area is positive
    if (calc_area(tris[i]) <= 0.) {
      tris[i].maxerr = 0.;//9.e99;
      tris[i].maxerrloc[0] = 1./3.;
      tris[i].maxerrloc[1] = 1./3.;
    }
    else compute_errorb(i,samples);
    if (debug > 0) std::cerr << "t " << i << " err " << tris[i].maxerr << " loc " << tris[i].maxerrloc[0] << " " << tris[i].maxerrloc[1] << std::endl;
    if (tris[i].maxerr > maxerr) {
      tri = i;
      maxerr = tris[i].maxerr;
      loc[0] = tris[i].maxerrloc[0];
      loc[1] = tris[i].maxerrloc[1];
    }
    if (tris[i].maxerr < minerr) minerr = tris[i].maxerr;
  }
  if (debug > 1) {
    std::cerr << "# t minerr " << minerr << " maxerr " << maxerr << std::endl << std::endl << std::endl;
  }

  return maxerr;
}

//
// calculate area of triangle (negative if nodes oriented clockwise)
//
double RMesh::calc_area(Tri & t)
{
  std::pair<double,double> xy[3];
  for (int i=0;i<3;i++) xy[i] = get_independent_vars(nodes[t.n[i]]);
  if (logvars == 1) {
    for (int i=0;i<3;i++) {
      xy[i].first = log(xy[i].first);
      xy[i].second = log(xy[i].second);
    }
  }

  if (debug > 1) {
    std::cerr << std::setprecision(15);
    for (int i=0;i<3;i++) std::cerr << "# ca p " << i << " xy " << xy[i].first << " " << xy[i].second << std::endl;
  }

  double x10 = xy[1].first-xy[0].first;
  double y10 = xy[1].second-xy[0].second;
  double x20 = xy[2].first-xy[0].first;
  double y20 = xy[2].second-xy[0].second;

  return (x10*y20-x20*y10)/2.;

}

void RMesh::printTris(std::ostream & out, const std::vector<Tri> & t, const std::vector<Node> & n)
{
  std::pair<double,double> xy[3];

  out << std::setprecision(15);
  for (std::size_t ti=0;ti<t.size();ti++) {
    for (int i=0;i<3;i++) xy[i] = get_independent_vars(n[t[ti].n[i]]);
    out << "# t " << ti << " n1 " << t[ti].n[0] << " n2 " << t[ti].n[1] << " n3 " << t[ti].n[2] << std::endl;
    if (logvars == 1) {
      out << "x1 " << log(xy[0].first) << " y1 " << log(xy[0].second) << std::endl;
      out << "x2 " << log(xy[1].first) << " y2 " << log(xy[1].second) << std::endl;
      out << "x3 " << log(xy[2].first) << " y3 " << log(xy[2].second) << std::endl;
      out << "x1 " << log(xy[0].first) << " y1 " << log(xy[0].second) << std::endl;
      out << std::endl;
    }
    else {
      out << "x1 " << xy[0].first << " y1 " << xy[0].second << std::endl;
      out << "x2 " << xy[1].first << " y2 " << xy[1].second << std::endl;
      out << "x3 " << xy[2].first << " y3 " << xy[2].second << std::endl;
      out << "x1 " << xy[0].first << " y1 " << xy[0].second << std::endl;
      out << std::endl;
    }
  }
  out << std::endl;
}

template<>
double convertRealAtoReal<double,dd_real>( const dd_real & var )
{
  return to_double(var);
}

//
// Find the boundary constained Delaunay triangulation of the
// region. Clears and then fills the tris vector and clears and fills
// the tris vector in each node.
//
void RMesh::find_tris()
{
  //
  // Clear the tris vector
  //
  tris.clear();

  //
  // Setup data for the delaunay triangulation routine
  //
  std::vector<double> px(nodes.size());
  std::vector<double> py(nodes.size());
  for (std::size_t i=0;i<nodes.size();i++) {
    std::pair<double,double> xy = get_independent_vars(nodes[i]);
    if (logvars == 1) {
      px[i] = log(xy.first);
      py[i] = log(xy.second);
    }
    else {
      px[i] = xy.first;
      py[i] = xy.second;
    }
    // clear the node tri list and error tallies
    nodes[i].tris.clear();
    nodes[i].maxerr = 0.;
    nodes[i].minerr = 9.99e99;
  }
  int nt = 0;
  unsigned int old_cw;
  fpu_fix_start(&old_cw);
  try {

    BoundDelaunayTriangulation<double,dd_real> bdt(px,py,nbnodes,xbounds[1]-xbounds[0],ybounds[1]-ybounds[0]);
    bdt.triangulate();
    const std::vector<int> & t = bdt.getTriangles();
    nt = t.size()/3;

    int debug = 0;
    if (debug > 1) std::cout << std::endl << std::endl << "# trioutstart" << std::endl << "# nt " << nt << std::endl;
    for (int i=0;i<nt;i++) {
      tris.push_back(Tri(t[3*i+0],t[3*i+1],t[3*i+2]));
    }
    std::sort(tris.begin(),tris.end());
    if (debug > 1) std::cout << std::endl << std::endl << "# trioutstart" << std::endl << "# nt " << nt << std::endl;
    for (std::size_t i=0;i<tris.size();i++) {
      nodes[tris[i].n[0]].tris.push_back(i);
      nodes[tris[i].n[1]].tris.push_back(i);
      nodes[tris[i].n[2]].tris.push_back(i);
      if (debug > 1) {
        std::cout << "### tri " << i << std::endl;
        print_tri(i);
      }
    }
    if (debug > 1) std::cout << "# trioutend" << std::endl << std::endl;
  }
  catch (std::runtime_error & e) {
    fpu_fix_end(&old_cw);
    std::cerr << e.what() << std::endl;
    std::cerr << std::setprecision(20) << "nn " << nodes.size() << " nbnodes " << nbnodes
              << " xb " << xbounds[0] << " " << xbounds[1]
              << " yb " << ybounds[0] << " " << ybounds[1] << std::endl;
    for (std::size_t i=0;i<nodes.size();i++) {
      std::cerr << "node " << i << " x " << px[i] << " y " << py[i] << std::endl;
    }
    throw std::runtime_error("BoundDelaunayTriangulation failed in RMesh::find_tris");
  }
  fpu_fix_end(&old_cw);
  if (nt < 1) {
    throw std::runtime_error("BoundDelaunayTriangulation failed to find any triangles");
  }

}

//
// return a vector of triangle indices that cannot reuse a prior error computation
//
void RMesh::find_unsaved_errors(const std::vector<int> & modnodes,std::vector<int> & unsavedtris)
{
  //
  // initialize return vector
  //
  unsavedtris.clear();

  //
  // find all triangles touching the modified nodes -- these must be recomputed
  //
  std::vector<int> comptris;
  for (std::size_t i=0;i<modnodes.size();i++) comptris.insert(comptris.end(),nodes[modnodes[i]].tris.begin(),nodes[modnodes[i]].tris.end());

  //
  // find all new triangles not in the saved set
  //
  for (std::size_t i=0;i<tris.size();i++) {
    TMeshTriData d(tris[i]);

    if (!std::binary_search(tmtd.begin(),tmtd.end(),d)) comptris.push_back(i);
  }

  //
  // sort and remove duplicate triangles for recomputation
  //
  std::sort(comptris.begin(),comptris.end());
  if (comptris.size() > 0) unsavedtris.push_back(comptris[0]);
  for (std::size_t i=1;i<comptris.size();i++) if (comptris[i] != unsavedtris.back()) unsavedtris.push_back(comptris[i]);

}

void RMesh::gettribccoords(Tri & t,double px,double py,double *l)
{
  std::pair<double,double> n[3];
  for (int i=0;i<3;i++) n[i] = get_independent_vars(nodes[t.n[i]]);

  double a[7];
  a[0] = n[0].first - n[2].first;
  a[1] = n[1].first - n[2].first;
  a[2] = n[2].first;
  a[3] = n[0].second - n[2].second;
  a[4] = n[1].second - n[2].second;
  a[5] = n[2].second;
  a[6] = 1./(a[0]*a[4]-a[1]*a[3]);

  l[0] = (a[1]*(a[5]-py)-(a[2]-px)*a[4])*a[6];
  l[1] = (a[3]*(a[2]-px)-(a[5]-py)*a[0])*a[6];
  l[2] = 1.-l[0]-l[1];
  if (fabs(l[0]) < 2.e-14) l[0] = 0.;
  if (fabs(l[1]) < 2.e-14) l[1] = 0.;
  if (fabs(l[2]) < 2.e-14) l[2] = 0.;
}

void RMesh::calc_errors_on_region(std::vector<int> & rnodes,double & maxerr,double & minerr)
{
  //
  // create triangles
  //
  find_tris();

  //
  // calculate errors on rnodes
  //
  maxerr = 0.;
  minerr = 9.e99;
  for (std::size_t i=0;i<tris.size();i++) {
    tris[i].maxerr = -1.;
    for (int j=0;j<3;j++) {
      for (std::size_t k=0;k<rnodes.size();k++) {
	if (tris[i].n[j] == rnodes[k]) {
	  if (tris[i].maxerr < 0.) compute_error(i);
	  if (debug > 1) std::cerr << "t " << i << " err " << tris[i].maxerr << std::endl;
	  if (tris[i].maxerr > maxerr) maxerr = tris[i].maxerr;
	  if (tris[i].maxerr < minerr) minerr = tris[i].maxerr;
	}
      }
    }
  }
  if (debug > 1) {
    std::cerr << "# t minerr " << minerr << " maxerr " << maxerr << std::endl << std::endl << std::endl;
  }

}

void RMesh::calc_errors_on_tris(std::vector<int> & t,double & maxerr,double & minerr)
{
  //
  // calculate errors on triangles
  //
  maxerr = 0.;
  minerr = 9.e99;
  for (std::size_t i=0;i<t.size();i++) {
    compute_error(t[i]);
    if (tris[t[i]].maxerr > maxerr) maxerr = tris[t[i]].maxerr;
    if (tris[t[i]].maxerr < minerr) minerr = tris[t[i]].maxerr;
  }
  if (debug > 1) {
    std::cerr << "# t minerr " << minerr << " maxerr " << maxerr << std::endl << std::endl << std::endl;
  }

}

void RMesh::adjust_node_maxerr(double x,double y,double & maxerr)
{
  //
  // reference to node that is being moved
  //
  Node & nn = nodes[a_ni];

  //
  // restore original node location for domain check
  //
  nn = a_n;

  //
  // check if new point in domain
  //
  bool inside = false;
  for (std::size_t i=0;i<a_t.size();i++) {
    double l[3];
    gettribccoords(a_t[i],x,y,l);
    if (l[0] >= 0.0 && l[0] <= 1.0 &&
	l[1] >= 0.0 && l[1] <= 1.0 &&
	l[2] >= 0.0 && l[2] <= 1.0 &&
	l[a_tbside[i]] > 0.01) {
      inside = true;
      break;
    }
  }
  if (!inside) {
    maxerr = 9.e99;
    return;
  }

  //
  // update node location
  //
  move_node(a_ni,x,y);

  double minerr;
  calc_errors_on_region(a_minreg,maxerr,minerr);

}

void RMesh::adjust_node_maxerr_top(double x,double y,double & maxerr)
{
  //
  // reference to node that is being moved
  //
  Node & nn = nodes[a_ni];

  //
  // restore original node location for domain check
  //
  nn = a_n;

  //
  // check if new point in domain
  //
  bool inside = false;
  for (std::size_t i=0;i<a_t.size();i++) {
    double l[3];
    gettribccoords(a_t[i],x,y,l);
    if (l[0] >= 0.0 && l[0] <= 1.0 &&
	l[1] >= 0.0 && l[1] <= 1.0 &&
	l[2] >= 0.0 && l[2] <= 1.0 &&
	l[a_tbside[i]] > 0.01) {
      inside = true;
      break;
    }
  }
  if (!inside) {
    maxerr = 9.e99;
    return;
  }

  //
  // update node location
  //
  move_node(a_ni,x,y);

  //
  // check that all tris are still anti-clockwise
  //
  for (std::size_t i=0;i<a_t.size();i++) {
    if (calc_area(a_t[i]) <= 0.) {
      maxerr = 8.e88;
      return;
    }
  }

  double minerr;
  calc_errors_on_tris(a_minreg,maxerr,minerr);

}

bool Trisortfun(Tri a,Tri b)
{
  if (a.n[0] < b.n[0]) return true;
  else if (a.n[0] == b.n[0]) {
    if (a.n[1] < b.n[1]) return true;
    else if (a.n[1] == b.n[1]) {
      if (a.n[2] < b.n[2]) return true;
      else return false;
    }
    else return false;
  }
  else return false;
}

void RMesh::getSubMeshState(const int ns,std::vector<double> & state)
{
  if (ns > int(nodes.size())) throw std::runtime_error("RMesh::getSubMeshState: requested more nodes than in mesh");

  //
  // ensure state vector is of proper length
  //
  state.resize(2*ns);

  //
  // get state vector info
  //
  int ntot = nodes.size();
  for (int i=0;i<ns;i++) {
    std::pair<double,double> iv = get_independent_vars(nodes[ntot-ns+i]);
    if (logvars) {
      state[2*i] = log(iv.first);
      state[2*i+1] = log(iv.second);
    }
    else {
      state[2*i] = iv.first;
      state[2*i+1] = iv.second;
    }
  }
}

bool RMesh::modify_submesh(RMeshMMOptData & d,const std::vector<double> & state)
{
  //
  // check if state size matches data
  //
  if (state.size() != 2*d.opt_nodes.size()) {
    throw std::runtime_error("RMesh::modify_submesh: state vector has different size than number of nodes to optimize");
  }

  //
  // update node locations with passed state
  //
  int offset = nodes.size()-d.opt_nodes.size();
  for (std::size_t i=0;i<d.opt_nodes.size();i++) {
    if (logvars) move_node(offset+i,exp(state[2*i]),exp(state[2*i+1]));
    else move_node(offset+i,state[2*i],state[2*i+1]);
  }
  
  //
  // triangulate
  //
  find_tris();

  //
  // Check Delaunay preservation outside domain of interest
  //
  std::vector<int> ntd;
  for (std::size_t i=0;i<tris.size();i++) {
    Tri & ct = tris[i];
    int min = 0;
    if (ct.n[1] < ct.n[min]) min = 1;
    if (ct.n[2] < ct.n[min]) min = 2;
    for (std::size_t j=0;j<d.deltris.size();j++) {
      if (ct.n[min] == d.deltris[j].n[0] &&
          ct.n[(min+1)%3] == d.deltris[j].n[1] &&
          ct.n[(min+2)%3] == d.deltris[j].n[2]) {
        ntd.push_back(j);
        break;
      }
    }
  }
  if (ntd.size() != d.deltris.size()) return false;

  //
  // save only new triangles in domain for error calc -- these are all
  // triangles not in the delaunay set ntd
  //
  int last = -1;
  std::sort(ntd.begin(),ntd.end());
  d.errtris.clear();
  for (std::size_t i=0;i<ntd.size();i++) {
    for (int j=last+1;j<ntd[i];j++) d.errtris.push_back(j);
    last = ntd[i];
  }
  for (int j=last+1;j<int(tris.size());j++) d.errtris.push_back(j);

  return true;
}
