/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SUTTG, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#ifndef BMESH_H
#define BMESH_H

#include <vector>
#include <cmath>

#include <pthread.h>
#include <semaphore.h>

#include "Node.H"
#include "TriMesh.H"
#include "ErrorMap.H"

class BMesh {
public:
  //
  // Constructor
  //
  BMesh( TriMesh * trimesh,
         std::vector<int> & boundarylist,
         TriMesh::MeshTypes meshtype,
         int logscale,
         int total_threads,
         const ErrorMap & error_map,
         int debug_level = 0 );

  //
  // Destructor
  //
  ~BMesh();

  //
  // method to compute error between given nodes
  //
  void compute_error( int node,
                      double & maxerr,
                      int samples = 10000 );

  //
  // method to compute error between given nodes with threaded node evaluations
  //
  void compute_error_thread( int node,
                             double & maxerr,
                             int samples = 10000 );

  //
  // clear all nodes except two fixed edge nodes
  //
  void clear_nodes() { 
    nodes1.erase(nodes1.begin()+1,nodes1.end()-1);
    if (nodes2.size() > 0) nodes2.erase(nodes2.begin()+1,nodes2.end()-1);
  }

  //
  // add a node to the mesh at the given location
  //
  void add_node( double x );

  //
  // move a node to the given location
  //
  void move_node( int node,
                  double x );

  //
  // get the nodes in the mesh
  //
  void get_nodes( const int phase,
                  std::vector<Node> & n );

  //
  // set node number n to the location x and return the error in the
  // interval between nodes numbered n and n-1, calculated given the
  // passed number of uniformly spaced samples
  // preconditions:
  // there must already exist at least n+1 nodes (including the two fixed boundary nodes)
  // n and samples must be positive
  // x must lie in the interval (x_{n-1},x_{n+1}) where x_n is in [0,1] and x_n < x_{n+1}
  // special case:
  // if x >= 1 and n == nodes1.size(), no node movement is performed and the error of the
  // last interval is returned
  //
  double setnodeloc( int n,
                     double x,
                     int samples = 100 );

  //
  // return number of nodes in the mesh
  //
  int get_num_nodes() const { return nodes1.size(); }

  //
  // start and stop the threads for evaluation of errors
  //
  void start_eval_threads();
  void stop_eval_threads();

private:

  //
  // nodes in the mesh
  //
  std::vector<Node> nodes1,nodes2;

  //
  // evaluate a node on the region associated with one of the boundaries
  //
  void evaluate_node(Node & node,
                     int nodevec = 1 );

  //
  // evaluate a node on the exact boundary
  //
  void evaluate_node_on_boundary( Node & node,
                                  int nodevec = 1 );

  //
  // boundaries this mesh represents
  //
  std::vector<int> boundaries;

  //
  // type of independent mesh vars (R,T) == 0, (R,E) == 1
  //
  TriMesh::MeshTypes type;

  //
  // interpolation function
  //
  void node_interp( double x,
                    Node & n1,
                    Node & n2,
                    Node & result );

  //
  // get the independent variable according to the type of boundary
  //
  double get_independent_var( Node & n );

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
  // Method to compute relative error of one node to another
  //
  double compute_rel_node_error( Node & n1,
                                 Node & n2 );

  //
  // Error evaluation map
  //
  const ErrorMap & emap;

  //
  // counter for performance testing
  //
  int numerrcalcs;

  //
  // total number of threads to start
  //
  int totthr;

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
  std::vector<std::pair<Node,int> > compnodes;

  //
  // thread start function
  //
  static void *threadfun( void * arg );

  //
  // thread work function
  //
  void nodes_compute();

  //
  // error flag for worker threads
  //
  int werror;

  //
  // error messages for each worker
  //
  std::vector<std::string> werrors;

};

#endif
