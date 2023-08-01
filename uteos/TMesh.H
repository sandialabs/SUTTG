/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SUTTG, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#ifndef TMESH_H
#define TMESH_H

#include <vector>
#include <cmath>
#include <iostream>

#include <pthread.h>
#include <semaphore.h>

#include "Node.H"
#include "TriMesh.H"
#include "BMesh.H"

struct TMeshTriData;

class TMesh {
public:
  //
  // Constructor
  //
  TMesh( TriMesh * trimesh,
         int region,
         TriMesh::MeshTypes meshtype,
         int logscale,
         int total_threads,
         const ErrorMap & error_map,
         int debug_level = 0 );

  //
  // Destructor
  //
  virtual ~TMesh();

  //
  // method to compute error on a triangle
  //
  void compute_error( int ti,
                      int samples = 4 );

  //
  // method to compute error on a triangle
  //
  void compute_errorb( int ti,
                       int samples = 4 );

  //
  // method to compute the mesh error using threaded error calcs
  //
  double compute_mesh_error_thread( int & tri,
                                    std::vector<double> & loc,
                                    const int samples,
                                    const int type );

  //
  // method to compute the mesh error using threaded error calcs and saved calculations
  //
  double compute_mesh_error_thread_saved( const std::vector<int> & comptris,
                                          int & tri,
                                          std::vector<double> & loc,
                                          const int samples,
                                          const int type );

  //
  // method to compute error over all triangles with precomputed samples
  //
  void compute_tris_error();

  //
  // add an empty node to the mesh -- must call move_node on it before
  // any other computations
  //
  void add_node_empty();

  //
  // add a node to the mesh at the given location
  //
  void add_node( double x,
                 double y );

  //
  // add a node to the mesh in the specified triangle and barycentric coordinates
  //
  void add_node( const int ti,
                 const std::vector<double> & bc );

  //
  // add a node to the mesh in the center of the specified triangle
  //
  void add_node( const int ti );

  //
  // add nodes to the mesh in the specified triangle in a grid pattern
  //
  void add_nodes( const int ti,
                  const int numn );

  //
  // add a node to the mesh in the center of the specified triangle
  // appending the location to the state vector
  //
  void add_node_state( const int ti,
                       std::vector<double> & states );

  //
  // move a node to the given location
  //
  void move_node( int node,
                  double x,
                  double y );

  //
  // move a node toward the given location
  //
  void move_node_toward( int node,
                         double x,
                         double y,
                         double d );

  //
  // update interior nodes with state array
  //
  void set_interior_nodes( const std::vector<double> & state );

  //
  // get the nodes in the mesh
  //
  void get_nodes( std::vector<Node> & n ) {
    n = nodes;
  }

  //
  // get the tris in the mesh
  //
  void get_tris( std::vector<Tri> & t ) {
    t = tris;
  }

  //
  // get the square of the angle between a point and the segment given
  // by the two specified nodes
  //
  double get_angle_sq( const int n1,
                       const int n2,
                       const double x,
                       const double y );

  //
  // get the nodes numbers of a triangle
  //
  void get_tri_nodes( const int nt,
                      std::vector<int> & n ) {
    if (nt >= int(tris.size())) throw std::runtime_error("TMesh::get_tri_nodes: nt out of bounds");
    n.resize(3);
    for (int i=0;i<3;i++) n[i] = tris[nt].n[i];
  }

  //
  // return number of nodes in the mesh
  //
  int get_num_nodes() {
    return nodes.size();
  }

  //
  // return number of boundary nodes in the mesh
  //
  int get_num_boundary_nodes() {
    return nbnodes;
  }

  //
  // return true for log scaling of parameters
  //
  bool log_scaling() {
    return logvars;
  }

  //
  // clears (erases) all interior nodes
  //
  void clear_interior_nodes() {
    nodes.erase(nodes.begin()+nbnodes,nodes.end());
  }

  //
  // initialize reusable points for error calculations
  //
  void initialize_error_points();

  //
  // initialize reusable points for error calculations with a given
  // sample size on each triangle
  //
  void initialize_error_points( const int samples );

  //
  // compute barycentric coords for the triangles
  //
  void compute_barycentric_coords();

  //
  // clear all points and memory from the precomputed cache
  //
  void clear_node_cache() { std::vector<EOSData>().swap(ncache); }

  //
  // print the mesh boundary for debugging
  //
  void print_boundary( std::ostream & ofile = std::cout );

  //
  // print the mesh triangles for debugging
  //
  void print_tris( std::ostream & ofile = std::cout ) {
    for (std::size_t i=0;i<tris.size();i++) print_tri(i,ofile);
    ofile << std::endl;
  }

  //
  // start and stop the threads for evaluation of errors
  //
  void start_eval_threads();
  void stop_eval_threads();

protected:
  //
  // nodes in the mesh
  //
  std::vector<Node> nodes;

  //
  // number of nodes on boundary of mesh (loaded in first part of nodes vector)
  //
  int nbnodes;

  //
  // evaluate a node
  //
  void evaluate_node( Node & node );

  //
  // evaluate a data point
  //
  void evaluate_point( EOSData & d );

  //
  // bounds on the mesh nodes in logvars coordinates
  //
  double xbounds[2],ybounds[2];

  //
  // triangles in the mesh
  //
  std::vector<Tri> tris;

  //
  // print the ith triangle independent variable coordinates
  //
  void print_tri( int i,
                  std::ostream & ofile = std::cout );

  //
  // return the tri number containing the point or -1 if not in any tris
  //
  int find_triangle( const double x,
                     const double y,
                     std::vector<double> & b,
                     const int debug = 0 );

  //
  // type of independent mesh vars (R,T) or (R,E)
  //
  TriMesh::MeshTypes type;

  //
  // interpolation function
  //
  void node_interp( double b1,
                    double b2,
                    Node & n1,
                    Node & n2,
                    Node & n3,
                    Node & result );

  //
  // set the independent variables
  //
  void set_independent_vars( double x,
                             double y,
                             Node & n );   // x->n
  void set_independent_vars( Node & n1,
                             Node & n2 ); // n1->n2

  //
  // get the independent variables
  //
  std::pair<double,double> get_independent_vars( const Node & n );

  //
  // debug output flag
  //
  int debug;

  //
  // flag for using log transformation of independent variables
  //
  int logvars;

  //
  // TriMesh this Mesh approximates
  //
  TriMesh * tm;

  //
  // Lookup number for the region of this mesh in tm
  //
  int tmregion;

  //
  // Method to compute relative error of one node to another
  //
  double compute_rel_node_error( Node & n1,
                                 Node & n2 );

  //
  // Error evaluation map
  //
  const ErrorMap & emap;

  //
  // vector of saved eos evaluations
  //
  std::vector<EOSData> ncache;

  //
  // total number of threads to start
  //
  int totthr;

  //
  // array of data for saving triangle error computations
  //
  std::vector<TMeshTriData> tmtd;

private:
  
  //
  // thread variables
  //
  int pthrid;
  std::vector<pthread_t> pthrs;

  //
  // thread synchronization support
  //
  pthread_mutex_t pthridmut;
  pthread_cond_t startcondb;
  pthread_mutex_t startmutb;
  sem_t * startsem;
  sem_t * endsem;

  //
  // nodes for evaluation and comparison
  //
  std::vector<Node> compnodes;

  //
  // load balancing support
  //
  std::vector<int> taskstatus;
  std::vector<pthread_mutex_t> taskmut;

  //
  // thread start function
  //
  static void *threadfun( void * arg );

  //
  // thread work function
  //
  void nodes_compute();
  void nodes_computeb();

  //
  // error flag for worker threads
  //
  int werror;

  //
  // error messages for each worker
  //
  std::vector<std::string> werrors;
  
};

struct TMeshTriData
{
  TMeshTriData() : maxerr(0.), maxerrloc0(0.), maxerrloc1(0.), n0(0), n1(0), n2(0) {}
  TMeshTriData( const Tri & t ) : maxerr(t.maxerr), 
                                  maxerrloc0(t.maxerrloc[0]), maxerrloc1(t.maxerrloc[1]),
                                  n0(t.n[0]), n1(t.n[1]), n2(t.n[2]) {}

  double maxerr;
  double maxerrloc0,maxerrloc1;
  int n0,n1,n2;
};

inline bool operator < ( const TMeshTriData & lhs, 
                         const TMeshTriData & rhs )
{
  if (lhs.n0 < rhs.n0) return true;
  else if (lhs.n0 == rhs.n0) {
    if (lhs.n1 < rhs.n1) return true;
    else if (lhs.n1 == rhs.n1) {
      if (lhs.n2 < rhs.n2) return true;
      else return false;
    }
    else return false;
  }
  else return false;
}

inline bool operator == ( const TMeshTriData & lhs,
                          const TMeshTriData & rhs )
{
  if (lhs.n0 == rhs.n0 && lhs.n1 == rhs.n1 && lhs.n2 == rhs.n2) return true;
  else return false;
}

inline bool operator != ( const TMeshTriData & lhs,
                          const TMeshTriData & rhs ) { return !operator==(lhs,rhs); }
inline bool operator >  ( const TMeshTriData & lhs,
                          const TMeshTriData & rhs ) { return  operator< (rhs,lhs); }
inline bool operator <= ( const TMeshTriData & lhs,
                          const TMeshTriData & rhs ) { return !operator> (lhs,rhs); }
inline bool operator >= ( const TMeshTriData & lhs,
                          const TMeshTriData & rhs ) { return !operator< (lhs,rhs); }

#endif
