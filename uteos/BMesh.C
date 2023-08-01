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
#include <sstream>

#include "BMesh.H"

//
// Constructor
//
BMesh::BMesh( TriMesh * trimesh,
              std::vector<int> & boundarylist,
              TriMesh::MeshTypes meshtype,
              int logscale,
              int total_threads,
              const ErrorMap & error_map,
              int debug_level ) :
  boundaries(boundarylist),type(meshtype),debug(debug_level),
  logvars(logscale),tm(trimesh),emap(error_map),
  totthr(total_threads),pthrid(0),pthrs(0),werrors(0)
{
  //
  // don't use threads
  //
  if (total_threads < 1) {
    totthr = 0;
  }

  //
  // Put in edge nodes
  //
  nodes1.push_back(tm->nodes[tm->boundaries[boundaries[0]].endnodes1[0]]);
  nodes1.push_back(tm->nodes[tm->boundaries[boundaries[0]].endnodes2[0]]);
  nodes1.front().x = 0;
  nodes1.back().x = 1;
  if (tm->boundaries[boundaries[0]].endnodes2.size() > 1) {
    nodes2.push_back(tm->nodes[tm->boundaries[boundaries[0]].endnodes1[1]]);
    nodes2.push_back(tm->nodes[tm->boundaries[boundaries[0]].endnodes2[1]]);
    nodes2.front().x = 0;
    nodes2.back().x = 1;
  }

}

//
// Destructor
//
BMesh::~BMesh()
{
  if (pthrs.size() > 0) stop_eval_threads();
}

void BMesh::start_eval_threads()
{
  //
  // error if we have not created this to do threads
  //
  if (totthr == 0) throw std::runtime_error("BMesh::start_eval_threads: at least one thread must be requested at BMesh construction");

  //
  // nothing to do if they are already started
  //
  if (pthrs.size() > 0) return;

  //
  // resize arrays
  //
  pthrs.resize(totthr);
  werrors.resize(totthr);

  //  
  // initialize syncronization vars
  //
  pthread_mutex_init(&pthridmut,NULL);
  pthread_mutex_init(&startmutb,NULL);
  pthread_cond_init(&startcondb,NULL);
#ifdef SEM_NAMED
  // Need to unlink first in case destructor did not get called at the end of
  // previous run, and the named semaphores did not get unlinked.
  sem_unlink("/bmeshstartsem");
  if ((startsem = sem_open("/bmeshstartsem",O_CREAT|O_EXCL,S_IRWXU,0)) == SEM_FAILED)
    throw std::runtime_error("sem_open failed on bmeshstartsem");
  sem_unlink("/bmeshendsem");
  if ((endsem = sem_open("/bmeshendsem",O_CREAT|O_EXCL,S_IRWXU,0)) == SEM_FAILED)
    throw std::runtime_error("sem_open failed on bmeshendsem");
#else
  startsem = new sem_t;
  endsem = new sem_t;
  if (sem_init(startsem,0,0) == -1)
    throw std::runtime_error("sem_init failed on bmeshstartsem");
  if (sem_init(endsem,0,0) == -1)
    throw std::runtime_error("sem_init failed on bmeshendsem");
#endif

  //
  // start worker threads
  //
  for (int i=0;i<totthr;i++) {
    if (pthread_create(&pthrs[i],NULL,&threadfun,this) != 0) {
      std::ostringstream oss;
      oss << "Error in BMesh: can't initialize thread " << i;
      throw std::runtime_error(oss.str());
    }
  }

}

void BMesh::stop_eval_threads()
{
  //
  // only stop threads if they have been started
  //
  if (pthrs.size() > 0) {
    for (int i=0;i<totthr;i++) {
      pthread_cancel(pthrs[i]);
    }
    for (int i=0;i<totthr;i++) {
      pthread_join(pthrs[i],NULL);
    }
    pthread_cond_destroy(&startcondb);
    pthread_mutex_destroy(&startmutb);
    pthread_mutex_destroy(&pthridmut);
#ifdef SEM_NAMED
    if (sem_close(startsem) == -1)
      throw std::runtime_error("sem_close failed on bmeshstartsem");
    if (sem_close(endsem) == -1)
      throw std::runtime_error("sem_close failed on bmeshendsem");
    if (sem_unlink("/bmeshstartsem") == -1)
      throw std::runtime_error("sem_unlink failed on bmeshstartsem");
    if (sem_unlink("/bmeshendsem") == -1)
      throw std::runtime_error("sem_unlink failed on bmeshendsem");
#else
    sem_destroy(startsem);
    sem_destroy(endsem);
    delete startsem;
    delete endsem;
#endif

    pthrs.clear();
    werrors.clear();
  }
}

void *BMesh::threadfun( void * arg )
{
  BMesh *p = (BMesh *)arg;
  p->nodes_compute();
  return NULL;
}

void BMesh::nodes_compute()
{
  int id;

  // get thread id
  pthread_mutex_lock(&pthridmut);
  id = pthrid;
  pthrid++;
  pthread_mutex_unlock(&pthridmut);

  while (1) {
    pthread_cleanup_push((void (*)(void *))pthread_mutex_unlock,(void *)&startmutb);
    pthread_mutex_lock(&startmutb);
    sem_post(startsem);
    pthread_cond_wait(&startcondb,&startmutb);
    pthread_mutex_unlock(&startmutb);
    pthread_cleanup_pop(0);

    Node neval(2,13);
    for (std::size_t i=id;i<compnodes.size()&&werror==0;i+=totthr) {
      neval.d.inputs = compnodes[i].first.d.inputs;
      try {
	evaluate_node(neval,compnodes[i].second);
      }
      catch (std::exception & e) {
	std::ostringstream oss;
	oss << "BMesh::nodes_compute: thread " << id 
            << " caught exception in evaluate_node, what(): " << e.what() << std::endl;
	werrors[id] += oss.str();
	werror = 1;
	break;
      }
      // save error in node x value for return to master routine
      compnodes[i].first.x = compute_rel_node_error(neval,compnodes[i].first);
    }
    
    sem_post(endsem);
  }

}

//
// method to compute error between given nodes with threaded node evaluations
//
void BMesh::compute_error_thread( int node,
                                  double & maxerr,
                                  int samples )
{
  
  if (totthr < 1) {
    //
    // use normal routine if not using threads
    //
    compute_error(node,maxerr,samples);

    return;
  }

  if (pthrs.size() == 0) throw std::runtime_error("BMesh::compute_error_thread: threads not initialized");

  numerrcalcs++;

  compnodes.clear();

  // node for interpolation
  Node n2(2,13);

  //
  // setup nodes for thread evaluation
  //
  for (int i=0;i<samples;i++) {
    //
    // the sample location between the two nodes
    //
    double x = (i+1.)/(samples+1.);

    //
    // n2 contains the interpolation result
    //
    node_interp(x,nodes1[node],nodes1[node+1],n2);

    //
    // put it on the list of nodes to compute
    //
    compnodes.push_back(std::pair<Node,int>(n2,1));
  }

  //
  // do it again for the second node vector
  //
  if (nodes2.size() != 0) {
    for (int i=0;i<samples;i++) {
      double x = (i+1.)/(samples+1.);

      node_interp(x,nodes2[node],nodes2[node+1],n2);

      compnodes.push_back(std::pair<Node,int>(n2,2));
    }
  }

  //
  // reset error info
  //
  werror = 0;
  for (int i=0;i<totthr;i++) werrors[i].clear();

  for (int j=0;j<totthr;j++) {
    sem_wait(startsem);
  }

  // start threads
  pthread_mutex_lock(&startmutb);
  pthread_cond_broadcast(&startcondb);
  pthread_mutex_unlock(&startmutb);

  for (int j=0;j<totthr;j++) {
    sem_wait(endsem);
  }

  //
  // check for errors in worker threads
  //
  if (werror != 0) {
    std::ostringstream oss;
    for (int i=0;i<totthr;i++) oss << werrors[i];
    if (debug > 0) std::cerr << oss.str() << std::endl;
    maxerr = 9.e99;
    return;
  }

  //
  // find max error
  //
  maxerr = compnodes[0].first.x;
  for (std::size_t i=1;i<compnodes.size();i++) {
    if (compnodes[i].first.x > maxerr) maxerr = compnodes[i].first.x;
  }

}

//
// return maximum error for all variables compared
// computes error in node 2 relative to 1
//
double BMesh::compute_rel_node_error( Node & n1,
                                      Node & n2 )
{
  return emap.compute_error(n1.d,n2.d,debug);
}

//
// method to compute error between given nodes
//
void BMesh::compute_error( int node,
                           double & maxerr,
                           int samples )
{

  numerrcalcs++;

  // zero out error
  maxerr = 0;

  // nodes for error computations
  Node n1(2,13),n2(2,13);

  for (int i=0;i<samples;i++) {
    //
    // the sample location between the two nodes
    //
    double x = (i+1.)/(samples+1.);

    //
    // n2 contains the interpolation result
    //
    node_interp(x,nodes1[node],nodes1[node+1],n2);

    //
    // n1 is the exact result evaluated at the interpolation point
    //
    n1.d.inputs = n2.d.inputs;
    try {
      evaluate_node(n1);
    }
    catch (std::runtime_error & e) {
      if (debug > 0) std::cerr << e.what() << std::endl;
      maxerr = 9.e99;
      continue;
    }

    //
    // Calculate relative error over all node quantities.  Regardless
    // of dependent variable type, it will give 0 relative error and
    // thus not contribute.
    //
    double relerr = emap.compute_error(n1.d,n2.d,debug);
    if (relerr > maxerr) maxerr = relerr;
  }

  if (nodes2.size() == 0) return;

  for (int i=0;i<samples;i++) {
    //
    // the sample location between the two nodes
    //
    double x = (i+1.)/(samples+1.);

    //
    // n2 contains the interpolation result
    //
    node_interp(x,nodes2[node],nodes2[node+1],n2);

    //
    // n1 is the exact result evaluated at the interpolation point
    //
    n1.d.inputs = n2.d.inputs;
    try {
      evaluate_node(n1,2);
    }
    catch (std::runtime_error & e) {
      if (debug > 0) std::cerr << e.what() << std::endl;
      maxerr = 9.e99;
      continue;
    }

    //
    // Calculate relative error over all node quantities.  Regardless
    // of dependent variable type, it will give 0 relative error and
    // thus not contribute.
    //
    double relerr = emap.compute_error(n1.d,n2.d,debug);
    if (relerr > maxerr) maxerr = relerr;
  }

  
}

//
// add a node to the mesh at the given location
//
void BMesh::add_node( double x )
{
  Node n(2,13);

  n.x = x;
  evaluate_node_on_boundary(n);
  
  nodes1.insert(nodes1.end()-1,n);

  if (nodes2.size() > 0) {
    evaluate_node_on_boundary(n,2);
  
    nodes2.insert(nodes2.end()-1,n);
  }
}

//
// move a node to the given location
//
void BMesh::move_node( int node,
                       double x )
{
  nodes1[node].x = x;
  evaluate_node_on_boundary(nodes1[node]);
  if (nodes2.size() > 0) {
    nodes2[node].x = x;
    evaluate_node_on_boundary(nodes2[node],2);
  }
}

//
// get the nodes in the mesh
//
void BMesh::get_nodes( const int phase,
                       std::vector<Node> & n )
{
  // decide which set of nodes to return
  Boundary & mb = tm->boundaries[boundaries[0]];
  if (mb.pbmnum < 0) {
    // phase is set in boundary
    if (mb.phase & phase) n = nodes1;
    else if (mb.phase2 & phase) n = nodes2;
    else throw std::runtime_error("BMesh::get_nodes: invalid phase in boundary");
  }
  else {
    // lookup phase in pbm
    int p1,pm,p2,pl,ph;
    tm->pbis[0][mb.pbmnum]->getPhases(p1,pm,p2,pl,ph);
    if ((mb.pbmside == 0 && phase & p1) || 
	(mb.pbmside == 2 && phase & p2)) n = nodes1;
    else if (phase & pm) n = nodes2;
    else throw std::runtime_error("BMesh::get_nodes: phase not found in pbm");
  }    
}

//
// evaluate a node
//
void BMesh::evaluate_node( Node & node,
                           int nodevec )
{
  // find the proper phase
  int bphase;
  Boundary & mb = tm->boundaries[boundaries[0]];
  if (mb.pbmnum < 0) {
    // phase is set in boundary
    if (nodevec == 1) bphase = mb.phase;
    else bphase = mb.phase2;
  }
  else {
    // lookup phase in pbm
    int p1,pm,p2,pl,ph;
    tm->pbis[0][mb.pbmnum]->getPhases(p1,pm,p2,pl,ph);
    if (nodevec == 1) {
      if (mb.pbmside == 0) bphase = p1;
      else bphase = p2;
    }
    else bphase = pm;
  }

  tm->evaluate_region(bphase,type,node.d);
}

//
// evaluate a node on the boundary
//
void BMesh::evaluate_node_on_boundary( Node & node,
                                       int nodevec )
{
  tm->evaluate_boundary(boundaries[0],nodevec,type,node.d,node.x);
}

//
// interpolation function
//
// Performs linear interpolation on all node quantities given the
// normalized distance x from node 1 to node 2 (i.e. 0<=x<=1)
// This assumes that all passed nodes have identically sized input and output vectors
//
void BMesh::node_interp( double x,
                         Node & n1,
                         Node & n2,
                         Node & result )
{
  if (logvars == 1) {
    // log scale the inputs
    for (std::size_t i=0;i<n1.d.inputs.size();i++)
      result.d.inputs[i] = exp(x*log(n2.d.inputs[i])+(1.-x)*log(n1.d.inputs[i]));
  }
  else {
    for (std::size_t i=0;i<n1.d.inputs.size();i++)
      result.d.inputs[i] = x*n2.d.inputs[i]+(1.-x)*n1.d.inputs[i];
  }

  for (std::size_t i=0;i<n1.d.outputs.size();i++)
    result.d.outputs[i] = x*n2.d.outputs[i]+(1.-x)*n1.d.outputs[i];

}

double BMesh::setnodeloc( int n,
                          double x,
                          int samples )
{
  //
  // sanity checks
  //
  if (n > get_num_nodes()) throw std::runtime_error("BMesh::setnodeloc: input n is too large");
  else if (n < 1)  throw std::runtime_error("BMesh::setnodeloc: input n is too small");
  if (samples < 1) throw std::runtime_error("BMesh::setnodeloc: at least one sample required");

  if (debug > 0) std::cerr << "BMesh::setnodeloc: n " << n << " x " << x << std::endl;

  double err = 9.e99;

  //
  // check that x is within needed bounds
  //
  if (n == get_num_nodes()) {
    //
    // check error at end node
    //
    if (x < 1.) throw std::runtime_error("BMesh::setnodeloc: tried to reset location of boundary node");
    
    compute_error_thread(n-2,err,samples);
  }
  else if (n == get_num_nodes()-1) {
    //
    // new node requested
    //

    if (x <= nodes1[n-1].x || x >= 1.) throw std::runtime_error("BMesh::setnodeloc: x value out of bounds for new node");

    add_node(x);

    compute_error_thread(n-1,err,samples);
  }
  else {
    //
    // moving a node
    //

    if (x <= nodes1[n-1].x || x >= nodes1[n+1].x) {
      std::cerr << "x " << x << " n " << n << " nx " << nodes1[n-1].x << " " << nodes1[n+1].x << std::endl;
      throw std::runtime_error("BMesh::setnodeloc: x value out of bounds for moved node");
    }

    move_node(n,x);

    compute_error_thread(n-1,err,samples);
  }    

  // print out nodes
  if (debug > 0) {
    std::cout << "BMesh::setnodeloc: nodes" << std::endl;
    for (std::size_t i=0;i<nodes1.size();i++) {
      std::cout << "node " << i << " x " << nodes1[i].x << " iv1 " << nodes1[i].d.inputs[0] << " iv2 " << nodes1[i].d.inputs[1] << std::endl;
    }
  }

  return err;
}
