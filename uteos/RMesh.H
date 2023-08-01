/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SUTTG, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#ifndef RMESH_H
#define RMESH_H

#include <vector>
#include <cmath>
#include <iostream>

#include "Node.H"
#include "TriMesh.H"
#include "BMesh.H"
#include "TMesh.H"
#include "ErrorMap.H"

class RMeshMMOptData
{
public:
  std::vector<int> opt_nodes;
  std::vector<int> boundary_nodes;
  std::vector<int> interior_nodes;
  std::vector<int> domain_tris;
  std::vector<int> delaunay_tris;
  std::vector<Tri> deltris;
  std::vector<int> errtris;
  std::vector<Tri> domtris;
  double xywidth[2];
  int opt_tri;
};

class RMesh : public TMesh
{
public:
  //
  // Constructor
  //
  RMesh( TriMesh * trimesh,
         int region,
         std::vector<BMesh *> & boundarymeshes,
         TriMesh::MeshTypes meshtype,
         int logscale,
         int tot_threads,
         const ErrorMap & error_map,
         int debug_level = 0 );

  //
  // Constructor for submesh
  //
  RMesh( RMesh * pmesh,
         const RMeshMMOptData & rmd );

  //
  // Destructor
  //
  ~RMesh();

  void printTris( std::ostream & out,
                  const std::vector<Tri> & t,
                  const std::vector<Node> & n );

  //
  // call back function for node optimization
  //
  void adjust_node_maxerr( double x,
                           double y,
                           double & maxerr );

  //
  // call back function for node optimization with fixed topology
  //
  void adjust_node_maxerr_top( double x,
                               double y,
                               double & maxerr );

  //
  // Method to compute the (boundary constrained) Delaunay
  // triangulation of the region
  //
  void find_tris();

  //
  // Method to compute the max error on the mesh and return the
  // associated triangle and barycentric coordinate location for the
  // max error
  //
  double compute_mesh_error( int & tri,
                             std::vector<double> & loc,
                             const int samples = 2 );

  //
  // Method to compute the max error on the mesh and return the
  // associated triangle and barycentric coordinate location for the
  // max error recomputing errors only for the passed triangles
  //
  double compute_mesh_error_saved( const std::vector<int> & comptris,
                                   int & tri,
                                   std::vector<double> & loc,
                                   const int samples = 2 );

  //
  // Method to compute and return the max error on the submesh
  // specified by the triangle list
  //
  double compute_submesh_error( const std::vector<int> & tlist,
                                const int samples = 2 );

  //
  // Method to compute the errors on the mesh
  //
  void compute_mesh_errors( const int samples );

  //
  // Modify a portion of the mesh described by the given data and state
  //
  bool modify_submesh( RMeshMMOptData & d,
                       const std::vector<double> & state );

  //
  // fill the state variable with the coordinates for the last ns nodes
  //
  void getSubMeshState( const int ns,
                        std::vector<double> & state );

  //
  // return a vector of triangle indices that need errors recomputated
  //
  void find_unsaved_errors( const std::vector<int> & modnodes,
                            std::vector<int> & unsavedtris );

private:
  //
  // calculate the errors for triangles touching a given set of nodes
  //
  void calc_errors_on_region( std::vector<int> & rnodes,
                              double & maxerr,
                              double & minerr );

  //
  // calculate the errors for a set of triangles
  //
  void calc_errors_on_tris( std::vector<int> & t,
                            double & maxerr,
                            double & minerr );

  //
  // get barycentric coordinates in a triangle
  //
  void gettribccoords( Tri & t,
                       double px,
                       double py,
                       double *l );

  //
  // calculate area of triangle (negative if nodes oriented clockwise)
  //
  double calc_area( Tri & t );

  //
  // variables to be used in the adjust_node_maxerr callback method
  //
  Node a_n;
  int a_ni;
  std::vector<Tri> a_t;
  std::vector<int> a_tbside;
  std::vector<int> a_minreg;

};

#endif
