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

#include "mpi.h"

#ifdef SEM_NAMED
#include <sys/stat.h>
#include <fcntl.h>
#endif

#include "TMesh.H"

//
// Constructor
//
TMesh::TMesh( TriMesh * trimesh,
              int region,
              TriMesh::MeshTypes meshtype,
              int logscale,
              int total_threads,
              const ErrorMap & error_map,
              int debug_level ) :
  type(meshtype),debug(debug_level),logvars(logscale),tm(trimesh),
  tmregion(region),emap(error_map),
  totthr(total_threads),pthrid(0),pthrs(0),taskmut(0),werrors(0)
{
  //
  // don't use threads
  //
  if (total_threads < 1) {
    totthr = 0;
  }

}

//
// Destructor
//
TMesh::~TMesh()
{
  if (pthrs.size() > 0) stop_eval_threads();
}

void TMesh::start_eval_threads()
{
  //
  // error if we have not created this to do threads
  //
  if (totthr == 0) throw std::runtime_error("TMesh::start_eval_threads: at least one thread must be requested at TMesh construction");

  //
  // nothing to do if they are already started
  //
  if (pthrs.size() > 0) return;

  //
  // resize arrays
  //
  pthrs.resize(totthr);
  werrors.resize(totthr);
  taskmut.resize(totthr);

  //  
  // initialize syncronization vars
  //
  pthread_mutex_init(&pthridmut,NULL);
  pthread_mutex_init(&startmutb,NULL);
  pthread_cond_init(&startcondb,NULL);
#ifdef SEM_NAMED
  // Need to unlink first in case destructor did not get called at the end of
  // previous run, and the named semaphores did not get unlinked.
  sem_unlink("/tmeshstartsem");
  if ((startsem = sem_open("/tmeshstartsem",O_CREAT|O_EXCL,S_IRWXU,0)) == SEM_FAILED)
    throw std::runtime_error("sem_open failed on tmeshstartsem");
  sem_unlink("/tmeshendsem");
  if ((endsem = sem_open("/tmeshendsem",O_CREAT|O_EXCL,S_IRWXU,0)) == SEM_FAILED)
    throw std::runtime_error("sem_open failed on tmeshendsem");
#else
  startsem = new sem_t;
  endsem = new sem_t;
  if (sem_init(startsem,0,0) == -1)
    throw std::runtime_error("sem_init failed on tmeshstartsem");
  if (sem_init(endsem,0,0) == -1)
    throw std::runtime_error("sem_init failed on tmeshendsem");
#endif

  //
  // start worker threads
  //
  for (int i=0;i<totthr;i++) {
    if (pthread_create(&pthrs[i],NULL,&threadfun,this) != 0) {
      std::ostringstream oss;
      oss << "Error in TMesh: can't initialize thread " << i;
      throw std::runtime_error(oss.str());
    }
  }

  //
  // initialize load balancing mutexes
  //
  for (int i=0;i<totthr;i++) pthread_mutex_init(&taskmut[i],NULL);

}

void TMesh::stop_eval_threads()
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
      throw std::runtime_error("sem_close failed on tmeshstartsem");
    if (sem_close(endsem) == -1)
      throw std::runtime_error("sem_close failed on tmeshendsem");
    if (sem_unlink("/tmeshstartsem") == -1)
      throw std::runtime_error("sem_unlink failed on tmeshstartsem");
    if (sem_unlink("/tmeshendsem") == -1)
      throw std::runtime_error("sem_unlink failed on tmeshendsem");
#else
    sem_destroy(startsem);
    sem_destroy(endsem);
    delete startsem;
    delete endsem;
#endif

    for (int i=0;i<totthr;i++) {
      pthread_mutex_destroy(&taskmut[i]);
    }

    pthrs.clear();
    werrors.clear();
    taskmut.clear();
  }
}

void *TMesh::threadfun( void * arg )
{
  TMesh *p = (TMesh *)arg;
  p->nodes_compute();
  return NULL;
}

void TMesh::nodes_computeb()
{
  int id;

  // get thread id
  pthread_mutex_lock(&pthridmut);
  id = pthrid;
  pthrid++;
  pthread_mutex_unlock(&pthridmut);

  // thread skip increment
  int tinc(1);
  if (totthr > 2) tinc = totthr/2+1;

  while (1) {
    pthread_cleanup_push((void (*)(void *))pthread_mutex_unlock,(void *)&startmutb);
    pthread_mutex_lock(&startmutb);
    sem_post(startsem);
    pthread_cond_wait(&startcondb,&startmutb);
    pthread_mutex_unlock(&startmutb);
    pthread_cleanup_pop(0);

    int curi(0);
    int curt(id);
    int xnum(compnodes.size());
    // only process if there are enough points
    if (curt+curi*totthr < xnum) {
      Node neval(2,13);
      while (1) {
        // find unprocessed point
        int done = 0;
        while (1) {
          // lock the points
          pthread_mutex_lock(&taskmut[curt]);
          while (1) {
            // search for untaken point
            if (curt+curi*totthr >= xnum) {
              break;
            }
            else if (taskstatus[curt+curi*totthr] == 0) {
              done = 1;
              // take the point
              taskstatus[curt+curi*totthr]++;
              break;
            }
            curi++;
          }
          // unlock the points
          pthread_mutex_unlock(&taskmut[curt]);
        
          if (done == 1) break;
          else {
            // go to next set of points
            curi = 0;
            curt = (curt+tinc)%totthr;
            if (curt == id) break;
          }
        }
      
        if (done == 0) {
          // no more points to process
          break;
        }
      
        // process point
        int ii = curt+curi*totthr;
        set_independent_vars(compnodes[ii],neval);
        try {
          evaluate_node(neval);
        }
        catch (std::exception & e) {
          std::ostringstream oss;
          oss << "TMesh::nodes_compute: thread " << id << " caught exception in evaluate_node, what(): " << e.what() << std::endl;
          werrors[id] += oss.str();
          werror = 1;
          break;
        }
        // save error in node x value for return to master routine
        compnodes[ii].x = compute_rel_node_error(neval,compnodes[ii]);

        curi++;
      }
    }
    
    sem_post(endsem);
  }

}

void TMesh::nodes_compute()
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
      set_independent_vars(compnodes[i],neval);
      try {
	evaluate_node(neval);
      }
      catch (std::exception & e) {
	std::ostringstream oss;
	oss << "TMesh::nodes_compute: thread " << id << " caught exception in evaluate_node, what(): " << e.what() << std::endl;
	werrors[id] += oss.str();
	werror = 1;
	break;
      }
      // save error in node x value for return to master routine
      compnodes[i].x = compute_rel_node_error(neval,compnodes[i]);
    }
    
    sem_post(endsem);
  }

}

int TMesh::find_triangle( const double x,
                          const double y,
                          std::vector<double> & b,
                          const int debug )
{
  int mytri = -1;
  for (std::size_t i=0;i<tris.size();i++) {
    tris[i].getBcoords(x,y,b);
    if (debug) print_tri(i);
    if (debug) std::cout << "# t " << i << " x " << x << " y " << y << " b " << b[0] << " " << b[1] << " " << b[2] << std::endl;
    if (b[0] >= 0. && b[0] <= (1.+2.e-14) &&
        b[1] >= 0. && b[1] <= (1.+2.e-14) &&
        b[2] >= 0. && b[2] <= (1.+2.e-14)) {
      mytri = i;
      break;
    }
  }
  return mytri;
}

void TMesh::initialize_error_points()
{
  if (debug > 1) std::cout << "TMesh::initialize_error_points: entry" << std::endl;
  ncache.clear();

  int totx = 100;
  int toty = 100;

  //
  // set point mesh increment
  //
  double dx = (xbounds[1]-xbounds[0])/(totx+1.);
  double dy = (ybounds[1]-ybounds[0])/(toty+1.);

  //
  // create points checking that each one lies in the mesh
  //
  EOSData d(2,13);
  std::vector<double> b(3);
  for (int j=0;j<toty;j++) {
    double y = ybounds[0]+dy*(j+1);
    for (int i=0;i<totx;i++) {
      double x = xbounds[0]+dx*(i+1);
      int myt = find_triangle(x,y,b);
      if (myt >= 0) {
        //
        // point is in a triangle and thus in mesh, so save
        //
        d.inputs[0] = x;
        d.inputs[1] = y;
        if (logvars) {
          d.inputs[0] = exp(x);
          d.inputs[1] = exp(y);
        }
        ncache.push_back(d);
      }
    }
    if (debug > 1) std::cout << "TMesh::initialize_error_points: j = " << j << " done" << std::endl;
  }

  //
  // evaluate all points
  //
  for (std::size_t i=0;i<ncache.size();i++) {
    evaluate_point(ncache[i]);
  }
  
}

void TMesh::initialize_error_points( const int samples )
{
  if (debug > 1) std::cout << "TMesh::initialize_error_points(samples): entry" << std::endl;
  ncache.clear();

  //
  // put in each interior node as a sample point
  //
  for (std::size_t i=nbnodes;i<nodes.size();i++) {
    ncache.push_back(nodes[i].d);
  }

  //
  // loop over each triangle putting in sample points except along
  // boundaries of region
  //
  for (std::size_t i=0;i<tris.size();i++) {
    // the triangle
    Tri & t = tris[i];
    
    // tri nodes
    Node & n1 = nodes[t.n[0]];
    Node & n2 = nodes[t.n[1]];
    Node & n3 = nodes[t.n[2]];
    Node nn(2,13);

    //
    // calculate sides
    //
    for (int j=0;j<samples;j++) {
      //
      // barycentric coords for a side
      //
      double b1 = (j+1)/(samples+1.);
      double b2 = 1.-b1;

      //
      // only add side points if it is not a boundary and second node
      // is higher number than first, to avoid double insert
      //
      if (t.n[1] > t.n[0] && (t.n[0] >= nbnodes || t.n[1] >= nbnodes)) {
        node_interp(b1,b2,n1,n2,n3,nn);
        ncache.push_back(nn.d);
      }
      if (t.n[2] > t.n[1] && (t.n[1] >= nbnodes || t.n[2] >= nbnodes)) {
        node_interp(b1,b2,n2,n3,n1,nn);
        ncache.push_back(nn.d);
      }
      if (t.n[0] > t.n[2] && (t.n[2] >= nbnodes || t.n[0] >= nbnodes)) {
        node_interp(b1,b2,n3,n1,n2,nn);
        ncache.push_back(nn.d);
      }
    }

    //
    // calculate interior
    //
    for (int k=0;k<samples-1;k++) {
      //
      // first barycentric coordinate
      //
      double b1 = (k+1)/(samples+1.);
      for (int j=0;j<samples-k-1;j++) {
        //
        // second barycentric coordinate
        //
        double b2 = (j+1)/(samples+1.);

        //
        // n21 contains the interpolation result
        //
        node_interp(b1,b2,n1,n2,n3,nn);
        ncache.push_back(nn.d);
      }
    }
  }

  //
  // evaluate all points
  //
  for (std::size_t i=0;i<ncache.size();i++) {
    evaluate_point(ncache[i]);
  }
  
}

//
// return maximum error for all variables compared
// computes error in node 2 relative to 1
//
double TMesh::compute_rel_node_error( Node & n1,
                                      Node & n2 )
{
  return emap.compute_error(n1.d,n2.d,debug);
}

//
// Method to compute error on a triangle formed by three nodes. Nodes
// are assumed to be in anti-clockwise order. Error is computed using
// a L-infinity norm over a uniform mesh on the triangle.
//
void TMesh::compute_error( int ti,
                           int samples )
{
  // get the triangle
  Tri & t = tris[ti];

  // zero out error and unflag location
  t.maxerr = 0;
  t.maxerrloc[0] = -9.e99;
  t.maxerrloc[1] = -9.e99;

  // tri nodes
  Node & n1 = nodes[t.n[0]];
  Node & n2 = nodes[t.n[1]];
  Node & n3 = nodes[t.n[2]];

  if (debug > 1) {
    std::pair<double,double> tmp = get_independent_vars(n1);
    std::cerr << "# n1 x " << tmp.first << " y " << tmp.second << std::endl;
    tmp = get_independent_vars(n2);
    std::cerr << "# n2 x " << tmp.first << " y " << tmp.second << std::endl;
    tmp = get_independent_vars(n3);
    std::cerr << "# n3 x " << tmp.first << " y " << tmp.second << std::endl;
  }

  // nodes for error computations
  Node n11(2,13),n21(2,13);
  Node n12(2,13),n22(2,13);
  Node n13(2,13),n23(2,13);

  //
  // calculate sides
  //
  for (int i=0;i<samples;i++) {
    //
    // barycentric coords for a side
    //
    double b1 = (i+1)/(samples+1.);
    double b2 = 1.-b1;

    //
    // n2[1-3] contains the interpolation result
    //
    node_interp(b1,b2,n1,n2,n3,n21);
    node_interp(b1,b2,n2,n3,n1,n22);
    node_interp(b1,b2,n3,n1,n2,n23);

    //
    // n1[1-3] is the exact result evaluated at the interpolation point
    //
    set_independent_vars(n21,n11);
    set_independent_vars(n22,n12);
    set_independent_vars(n23,n13);
    evaluate_node(n11);
    evaluate_node(n12);
    evaluate_node(n13);

    //
    // Calculate relative error over all node quantities.  Regardless
    // of dependent variable type, it will give 0 relative error and
    // thus not contribute.
    //
    if (debug > 1) std::cerr << "i " << i << " s 1 b1 " << b1 << " b2 " << b2 << std::endl;
    double relerr = compute_rel_node_error(n11,n21);
    if (relerr > t.maxerr) {
      t.maxerr = relerr;
      std::pair<double,double> xy = get_independent_vars(n21);
      t.maxerrloc[0] = xy.first;
      t.maxerrloc[1] = xy.second;
    }
    if (debug > 1) std::cerr << "i " << i << " s 2 b1 " << b1 << " b2 " << b2 << std::endl;
    relerr = compute_rel_node_error(n12,n22);
    if (relerr > t.maxerr) {
      t.maxerr = relerr;
      std::pair<double,double> xy = get_independent_vars(n22);
      t.maxerrloc[0] = xy.first;
      t.maxerrloc[1] = xy.second;
    }
    if (debug > 1) std::cerr << "i " << i << " s 3 b1 " << b1 << " b2 " << b2 << std::endl;
    relerr = compute_rel_node_error(n13,n23);
    if (relerr > t.maxerr) {
      t.maxerr = relerr;
      std::pair<double,double> xy = get_independent_vars(n23);
      t.maxerrloc[0] = xy.first;
      t.maxerrloc[1] = xy.second;
    }
  }

  //
  // perform same error calculations for interior nodes
  //
  for (int i=0;i<samples-1;i++) {
    //
    // first barycentric coordinate
    //
    double b1 = (i+1)/(samples+1.);
    for (int j=0;j<samples-i-1;j++) {
      //
      // second barycentric coordinate
      //
      double b2 = (j+1)/(samples+1.);

      //
      // n21 contains the interpolation result
      //
      node_interp(b1,b2,n1,n2,n3,n21);

      //
      // n11 is the exact result evaluated at the interpolation point
      //
      set_independent_vars(n21,n11);
      evaluate_node(n11);

      //
      // Calculate relative error over all node quantities.  Regardless
      // of dependent variable type, it will give 0 relative error and
      // thus not contribute.
      //
      if (debug > 1) std::cerr << "i " << i << " j " << j << " b1 " << b1 << " b2 " << b2 << std::endl;
      double relerr = compute_rel_node_error(n11,n21);
      if (relerr > t.maxerr) {
	t.maxerr = relerr;
	std::pair<double,double> xy = get_independent_vars(n21);
	t.maxerrloc[0] = xy.first;
	t.maxerrloc[1] = xy.second;
      }
    }
  }

  if (debug > 1) {
    print_tri(ti);
    std::cerr << "# compute_error " << ti << " e " << t.maxerr << std::endl;
  }
}

//
// Method to compute error on a triangle formed by three nodes. Nodes
// are assumed to be in anti-clockwise order. Error is computed using
// a L-infinity norm over a uniform mesh on the triangle. Error
// location stored using barycentric location instead of real location
// as in compute_error.
//
void TMesh::compute_errorb( int ti,
                            int samples )
{
  // get the triangle
  Tri & t = tris[ti];

  // zero out error and unflag location
  t.maxerr = 0;
  t.maxerrloc[0] = -9.e99;
  t.maxerrloc[1] = -9.e99;

  // tri nodes
  Node & n1 = nodes[t.n[0]];
  Node & n2 = nodes[t.n[1]];
  Node & n3 = nodes[t.n[2]];

  if (debug > 1) {
    std::pair<double,double> tmp = get_independent_vars(n1);
    std::cerr << "# n1 x " << tmp.first << " y " << tmp.second << std::endl;
    tmp = get_independent_vars(n2);
    std::cerr << "# n2 x " << tmp.first << " y " << tmp.second << std::endl;
    tmp = get_independent_vars(n3);
    std::cerr << "# n3 x " << tmp.first << " y " << tmp.second << std::endl;
  }

  // nodes for error computations
  Node n11(2,13),n21(2,13);
  Node n12(2,13),n22(2,13);
  Node n13(2,13),n23(2,13);

  //
  // calculate sides
  //
  for (int i=0;i<samples;i++) {
    //
    // barycentric coords for a side
    //
    double b1 = (i+1)/(samples+1.);
    double b2 = 1.-b1;

    //
    // n2[1-3] contains the interpolation result
    //
    node_interp(b1,b2,n1,n2,n3,n21);
    node_interp(b1,b2,n2,n3,n1,n22);
    node_interp(b1,b2,n3,n1,n2,n23);

    //
    // n1[1-3] is the exact result evaluated at the interpolation point
    //
    set_independent_vars(n21,n11);
    set_independent_vars(n22,n12);
    set_independent_vars(n23,n13);
    try {
      evaluate_node(n11);
      evaluate_node(n12);
      evaluate_node(n13);
    }
    catch (std::runtime_error & e) {
      if (debug > 0) std::cerr << e.what() << std::endl;
      t.maxerr = 9.e99;
      t.maxerrloc[0] = 1./3.;
      t.maxerrloc[1] = 1./3.;
      return;
    }

    //
    // Calculate relative error over all node quantities.  Regardless
    // of dependent variable type, it will give 0 relative error and
    // thus not contribute.
    //
    if (debug > 1) std::cerr << "i " << i << " s 1 b1 " << b1 << " b2 " << b2 << std::endl;
    double relerr = compute_rel_node_error(n11,n21);
    if (relerr > t.maxerr) {
      t.maxerr = relerr;
      t.maxerrloc[0] = b1;
      t.maxerrloc[1] = b2;
    }
    if (debug > 1) std::cerr << "i " << i << " s 2 b1 " << b1 << " b2 " << b2 << std::endl;
    relerr = compute_rel_node_error(n12,n22);
    if (relerr > t.maxerr) {
      t.maxerr = relerr;
      t.maxerrloc[0] = 1.-b1-b2;
      t.maxerrloc[1] = b1;
    }
    if (debug > 1) std::cerr << "i " << i << " s 3 b1 " << b1 << " b2 " << b2 << std::endl;
    relerr = compute_rel_node_error(n13,n23);
    if (relerr > t.maxerr) {
      t.maxerr = relerr;
      t.maxerrloc[0] = b2;
      t.maxerrloc[1] = 1.-b1-b2;
    }
  }

  //
  // perform same error calculations for interior nodes
  //
  for (int i=0;i<samples-1;i++) {
    //
    // first barycentric coordinate
    //
    double b1 = (i+1)/(samples+1.);
    for (int j=0;j<samples-i-1;j++) {
      //
      // second barycentric coordinate
      //
      double b2 = (j+1)/(samples+1.);

      //
      // n21 contains the interpolation result
      //
      node_interp(b1,b2,n1,n2,n3,n21);

      //
      // n11 is the exact result evaluated at the interpolation point
      //
      set_independent_vars(n21,n11);
      try {
        evaluate_node(n11);
      }
      catch (std::runtime_error & e) {
        if (debug > 0) std::cerr << e.what() << std::endl;
        t.maxerr = 9.e99;
        t.maxerrloc[0] = 1./3.;
        t.maxerrloc[1] = 1./3.;
        return;
      }

      //
      // Calculate relative error over all node quantities.  Regardless
      // of dependent variable type, it will give 0 relative error and
      // thus not contribute.
      //
      if (debug > 1) std::cerr << "i " << i << " j " << j << " b1 " << b1 << " b2 " << b2 << std::endl;
      double relerr = compute_rel_node_error(n11,n21);
      if (relerr > t.maxerr) {
	t.maxerr = relerr;
	t.maxerrloc[0] = b1;
	t.maxerrloc[1] = b2;
      }
    }
  }

  if (debug > 1) {
    print_tri(ti);
    std::cerr << "# compute_error " << ti << " e " << t.maxerr << std::endl;
  }
}

double TMesh::compute_mesh_error_thread( int & tri,
                                         std::vector<double> & loc,
                                         const int samples,
                                         const int type )
{
  //
  // set up return variables
  //
  loc.resize(2);
  double maxerr = 0.;
  double minerr = 9.e99;
  loc[0] = -9.e99;
  loc[1] = -9.e99;
  tri = -1;

  if (totthr < 1) {
    //
    // normal loop if no threads
    //
    for (std::size_t i=0;i<tris.size();i++) {
      if (type == 0) compute_error(i,samples);
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
    if (debug > 0) {
      std::cerr << "# t minerr " << minerr << " maxerr " << maxerr << std::endl << std::endl << std::endl;
    }
  }
  else {
    if (pthrs.size() == 0) throw std::runtime_error("TMesh::compute_mesh_error_thread: threads not initialized");

    //
    // setup nodes for thread evaluation
    //
    int nsamp = (samples+2)*(samples+3)/2-3;
    compnodes.clear();//resize(nsamp*tris.size());

    //
    // set comparison nodes for the triangle and store the location in
    // minerr and maxerr
    //
    for (std::size_t ii=0;ii<tris.size();ii++) {
      //
      // shortcut to current triangle
      //
      Tri & t = tris[ii];

      //
      // zero out error and unflag location
      //
      t.maxerr = 0;
      t.maxerrloc[0] = -9.e99;
      t.maxerrloc[1] = -9.e99;

      //
      // tri nodes
      //
      Node & n1 = nodes[t.n[0]];
      Node & n2 = nodes[t.n[1]];
      Node & n3 = nodes[t.n[2]];

      if (debug > 1) {
	std::pair<double,double> tmp = get_independent_vars(n1);
	std::cerr << "# n1 x " << tmp.first << " y " << tmp.second << std::endl;
	tmp = get_independent_vars(n2);
	std::cerr << "# n2 x " << tmp.first << " y " << tmp.second << std::endl;
	tmp = get_independent_vars(n3);
	std::cerr << "# n3 x " << tmp.first << " y " << tmp.second << std::endl;
      }

      //
      // nodes for error computations
      //
      Node n21(2,13);
      Node n22(2,13);
      Node n23(2,13);

      //
      // calculate sides
      //
      for (int i=0;i<samples;i++) {
	//
	// barycentric coords for a side
	//
	double b1 = (i+1)/(samples+1.);
	double b2 = 1.-b1;

	//
	// n2[1-3] contains the interpolation result
	//
	node_interp(b1,b2,n1,n2,n3,n21);
	node_interp(b1,b2,n2,n3,n1,n22);
	node_interp(b1,b2,n3,n1,n2,n23);

	//
	// save location in nodes
	//
	if (type == 0) {
	  std::pair<double,double> xy = get_independent_vars(n21);
	  n21.minerr = xy.first;
	  n21.maxerr = xy.second;
	  xy = get_independent_vars(n22);
	  n22.minerr = xy.first;
	  n22.maxerr = xy.second;
	  xy = get_independent_vars(n23);
	  n23.minerr = xy.first;
	  n23.maxerr = xy.second;
	}
	else {
	  n21.minerr = b1;
	  n21.maxerr = b2;
	  n22.minerr = 1.-b1-b2;
	  n22.maxerr = b1;
	  n23.minerr = b2;
	  n23.maxerr = 1.-b1-b2;
	}

	//
	// save to the computation nodes
	//
	compnodes.push_back(n21);
	compnodes.push_back(n22);
	compnodes.push_back(n23);
      }

      //
      // calculate interior nodes
      //
      for (int i=0;i<samples-1;i++) {
	//
	// first barycentric coordinate
	//
	double b1 = (i+1)/(samples+1.);
	for (int j=0;j<samples-i-1;j++) {
	  //
	  // second barycentric coordinate
	  //
	  double b2 = (j+1)/(samples+1.);

	  //
	  // n21 contains the interpolation result
	  //
	  node_interp(b1,b2,n1,n2,n3,n21);

	  //
	  // save location in node
	  //
	  if (type == 0) {
	    std::pair<double,double> xy = get_independent_vars(n21);
	    n21.minerr = xy.first;
	    n21.maxerr = xy.second;
	  }
	  else {
	    n21.minerr = b1;
	    n21.maxerr = b2;
	  }

	  //
	  // save to computation nodes
	  //
	  compnodes.push_back(n21);
	}
      }
    }

    //
    // reset error info
    //
    werror = 0;
    for (int i=0;i<totthr;i++) werrors[i].clear();

    //
    // reset node computation status
    //
    taskstatus.assign(compnodes.size(),0);

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
      throw std::runtime_error(oss.str());
    }

    //
    // set triangle errors
    //
    for (std::size_t i=0;i<tris.size();i++) {
      //
      // index to first evaluation for this triangle
      //
      int k = i*nsamp;

      //
      // find max error on this triangle
      //
      tris[i].maxerr = compnodes[k].x;
      tris[i].maxerrloc[0] = compnodes[k].minerr;
      tris[i].maxerrloc[1] = compnodes[k].maxerr;
      for (int j=1;j<nsamp;j++) {
	if (compnodes[k+j].x > tris[i].maxerr) {
	  tris[i].maxerr = compnodes[k+j].x;
	  tris[i].maxerrloc[0] = compnodes[k+j].minerr;
	  tris[i].maxerrloc[1] = compnodes[k+j].maxerr;
	}
      }
      
      //
      // update mesh max error
      //
      if (debug > 0) std::cerr << "t " << i << " err " << tris[i].maxerr << " loc " << tris[i].maxerrloc[0] << " " << tris[i].maxerrloc[1] << std::endl;
      if (tris[i].maxerr > maxerr) {
	tri = i;
	maxerr = tris[i].maxerr;
	loc[0] = tris[i].maxerrloc[0];
	loc[1] = tris[i].maxerrloc[1];
      }
      if (tris[i].maxerr < minerr) minerr = tris[i].maxerr;
    }
    if (debug > 0) {
      std::cerr << "# t minerr " << minerr << " maxerr " << maxerr << std::endl << std::endl << std::endl;
    }
  }

  return maxerr;

}

double TMesh::compute_mesh_error_thread_saved( const std::vector<int> & comptris,
                                               int & tri,
                                               std::vector<double> & loc,
                                               const int samples,
                                               const int type )
{
  //
  // set up return variables
  //
  loc.resize(2);
  double maxerr = 0.;
  double minerr = 9.e99;
  loc[0] = -9.e99;
  loc[1] = -9.e99;
  tri = -1;

  if (totthr < 1) {
    //
    // normal loop if no threads
    //
    std::size_t cti = 0;
    for (std::size_t i=0;i<tris.size();i++) {
      if (cti < comptris.size() && comptris[cti] == int(i)) {
	//
	// recompute this tri
	//
	if (type == 0) compute_error(i,samples);
	else compute_errorb(i,samples);
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
	  else throw std::runtime_error("TMesh::compute_mesh_error_thread_saved: failed to find a tri in the saved vector");
	}
	else throw std::runtime_error("TMesh::compute_mesh_error_thread_saved: failed to find a tri in the saved vector (off end)");
      }

      if (debug > 1) std::cerr << "ts " << i << " err " << tris[i].maxerr << " loc " << tris[i].maxerrloc[0] << " " << tris[i].maxerrloc[1] << std::endl;
      if (tris[i].maxerr > maxerr) {
	tri = i;
	maxerr = tris[i].maxerr;
	loc[0] = tris[i].maxerrloc[0];
	loc[1] = tris[i].maxerrloc[1];
      }
      if (tris[i].maxerr < minerr) minerr = tris[i].maxerr;
    }
    if (debug > 0) {
      std::cerr << "# t minerr " << minerr << " maxerr " << maxerr << std::endl << std::endl << std::endl;
    }
  }
  else {
    if (pthrs.size() == 0) throw std::runtime_error("TMesh::compute_mesh_error_thread: threads not initialized");

    //
    // setup nodes for thread evaluation
    //
    int nsamp = (samples+2)*(samples+3)/2-3;
    compnodes.clear();//resize(nsamp*tris.size());

    //
    // set comparison nodes for the triangle and store the location in
    // minerr and maxerr
    //
    std::size_t cti = 0;
    for (std::size_t ii=0;ii<tris.size();ii++) {
      //
      // shortcut to current triangle
      //
      Tri & t = tris[ii];

      //
      // zero out error and unflag location
      //
      t.maxerr = 0;
      t.maxerrloc[0] = -9.e99;
      t.maxerrloc[1] = -9.e99;

      //
      // skip this one if error saved
      //
      if (cti < comptris.size() && comptris[cti] == int(ii)) {
	//
	// recompute this tri
	//
	cti++;
      }
      else {
	//
	// get saved error values for this tri
	//
	TMeshTriData d(t);
	std::vector<TMeshTriData>::iterator tmtdi = std::lower_bound(tmtd.begin(),tmtd.end(),d);
	if (tmtdi != tmtd.end()) {
	  if (tmtdi->n0 == t.n[0] && tmtdi->n1 == t.n[1] && tmtdi->n2 == t.n[2]) {
	    t.maxerr = tmtdi->maxerr;
	    t.maxerrloc[0] = tmtdi->maxerrloc0;
	    t.maxerrloc[1] = tmtdi->maxerrloc1;
	  }
	  else throw std::runtime_error("TMesh::compute_mesh_error_thread_saved: failed to find a tri in the saved vector");
	}
	else throw std::runtime_error("TMesh::compute_mesh_error_thread_saved: failed to find a tri in the saved vector (off end)");
        continue;
      }
      
      //
      // tri nodes
      //
      Node & n1 = nodes[t.n[0]];
      Node & n2 = nodes[t.n[1]];
      Node & n3 = nodes[t.n[2]];

      if (debug > 1) {
	std::pair<double,double> tmp = get_independent_vars(n1);
	std::cerr << "# n1 x " << tmp.first << " y " << tmp.second << std::endl;
	tmp = get_independent_vars(n2);
	std::cerr << "# n2 x " << tmp.first << " y " << tmp.second << std::endl;
	tmp = get_independent_vars(n3);
	std::cerr << "# n3 x " << tmp.first << " y " << tmp.second << std::endl;
      }

      //
      // nodes for error computations
      //
      Node n21(2,13);
      Node n22(2,13);
      Node n23(2,13);

      //
      // calculate sides
      //
      for (int i=0;i<samples;i++) {
	//
	// barycentric coords for a side
	//
	double b1 = (i+1)/(samples+1.);
	double b2 = 1.-b1;

	//
	// n2[1-3] contains the interpolation result
	//
	node_interp(b1,b2,n1,n2,n3,n21);
	node_interp(b1,b2,n2,n3,n1,n22);
	node_interp(b1,b2,n3,n1,n2,n23);

	//
	// save location in nodes
	//
	if (type == 0) {
	  std::pair<double,double> xy = get_independent_vars(n21);
	  n21.minerr = xy.first;
	  n21.maxerr = xy.second;
	  xy = get_independent_vars(n22);
	  n22.minerr = xy.first;
	  n22.maxerr = xy.second;
	  xy = get_independent_vars(n23);
	  n23.minerr = xy.first;
	  n23.maxerr = xy.second;
	}
	else {
	  n21.minerr = b1;
	  n21.maxerr = b2;
	  n22.minerr = 1.-b1-b2;
	  n22.maxerr = b1;
	  n23.minerr = b2;
	  n23.maxerr = 1.-b1-b2;
	}

	//
	// save to the computation nodes
	//
	compnodes.push_back(n21);
	compnodes.push_back(n22);
	compnodes.push_back(n23);
      }

      //
      // calculate interior nodes
      //
      for (int i=0;i<samples-1;i++) {
	//
	// first barycentric coordinate
	//
	double b1 = (i+1)/(samples+1.);
	for (int j=0;j<samples-i-1;j++) {
	  //
	  // second barycentric coordinate
	  //
	  double b2 = (j+1)/(samples+1.);

	  //
	  // n21 contains the interpolation result
	  //
	  node_interp(b1,b2,n1,n2,n3,n21);

	  //
	  // save location in node
	  //
	  if (type == 0) {
	    std::pair<double,double> xy = get_independent_vars(n21);
	    n21.minerr = xy.first;
	    n21.maxerr = xy.second;
	  }
	  else {
	    n21.minerr = b1;
	    n21.maxerr = b2;
	  }

	  //
	  // save to computation nodes
	  //
	  compnodes.push_back(n21);
	}
      }
    }

    //
    // reset error info
    //
    werror = 0;
    for (int i=0;i<totthr;i++) werrors[i].clear();

    //
    // reset node computation status
    //
    taskstatus.assign(compnodes.size(),0);

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
      throw std::runtime_error(oss.str());
    }

    //
    // set triangle errors
    //
    for (std::size_t i=0;i<comptris.size();i++) {
      //
      // index to first evaluation for this triangle
      //
      int k = i*nsamp;

      //
      // this triangle
      //
      Tri & t = tris[comptris[i]];

      //
      // find max error on this triangle
      //
      t.maxerr = compnodes[k].x;
      t.maxerrloc[0] = compnodes[k].minerr;
      t.maxerrloc[1] = compnodes[k].maxerr;
      for (int j=1;j<nsamp;j++) {
	if (compnodes[k+j].x > t.maxerr) {
	  t.maxerr = compnodes[k+j].x;
	  t.maxerrloc[0] = compnodes[k+j].minerr;
	  t.maxerrloc[1] = compnodes[k+j].maxerr;
	}
      }
    }

    for (std::size_t i=0;i<tris.size();i++) {
      //
      // update mesh max error
      //
      if (debug > 1) std::cerr << "tt " << i << " err " << tris[i].maxerr << " loc " << tris[i].maxerrloc[0] << " " << tris[i].maxerrloc[1] << std::endl;
      if (tris[i].maxerr > maxerr) {
	tri = i;
	maxerr = tris[i].maxerr;
	loc[0] = tris[i].maxerrloc[0];
	loc[1] = tris[i].maxerrloc[1];
      }
      if (tris[i].maxerr < minerr) minerr = tris[i].maxerr;
    }
    if (debug > 0) {
      std::cerr << "# t minerr " << minerr << " maxerr " << maxerr << std::endl << std::endl << std::endl;
    }
  }

  //
  // now save these new tris
  //
  tmtd.clear();
  for (std::size_t i=0;i<tris.size();i++) tmtd.push_back(TMeshTriData(tris[i]));

  return maxerr;

}

//
// Method to compute error on each triangle in mesh using a L-infinity
// norm over a precomputes set of sample points.
//
void TMesh::compute_tris_error()
{
  //
  // loop over all sample points
  //
  for (std::size_t i=0;i<ncache.size();i++) {
    //
    // find the triangle for this point
    //
    double x = ncache[i].inputs[0];
    double y = ncache[i].inputs[1];
    if (logvars) {
      x = log(x);
      y = log(y);
    }
    std::vector<double> b(3);
    int t = find_triangle(x,y,b);
    if (t < 0) {
      find_triangle(x,y,b,1);
      throw std::runtime_error("TMesh::compute_tris_error: sample point not in any triangle");
    }

    //
    // Interpolation result for this point
    //
    Node ni(2,13);
    node_interp(b[0],b[1],nodes[tris[t].n[0]],nodes[tris[t].n[1]],nodes[tris[t].n[2]],ni);
   
    //
    // calculate error at this point
    // 
    if (debug > 1) std::cerr << "point i " << i << " t " << t << " b1 " << b[0] << " b2 " << b[1] << std::endl;
    double relerr = emap.compute_error(ncache[i],ni.d,debug);

    //
    // update triangle quantities
    //
    if (relerr > tris[t].maxerr) {
      tris[t].maxerr = relerr;
      tris[t].maxerrloc[0] = b[0];
      tris[t].maxerrloc[1] = b[1];
    }
  }
}

//
// add a node to the mesh in the specified triangle and barycentric coordinates
//
void TMesh::add_node( const int ti,
                      const std::vector<double> & bc )
{
  int tsize = tris.size();
  if (ti < 0 || ti > tsize) throw std::runtime_error("TMesh::add_node: invalid triangle number");

  //
  // get the triangle
  //
  Tri & t = tris[ti];

  //
  // tri nodes
  //
  Node & n1 = nodes[t.n[0]];
  Node & n2 = nodes[t.n[1]];
  Node & n3 = nodes[t.n[2]];

  //
  // new node
  //
  Node n(2,13);

  //
  // interpolate for the new node location
  //
  node_interp(bc[0],bc[1],n1,n2,n3,n);

  //
  // finish evaluation and insert
  //
  evaluate_node(n);
  nodes.push_back(n);
}

//
// add a node to the mesh in the center of the specified triangle
//
void TMesh::add_node( const int ti )
{
  std::vector<double> bc(2);
  bc[0] = 1./3.;
  bc[1] = 1./3.;
  add_node(ti,bc);
}

//
// add nodes to the mesh in the specified triangle in a grid pattern
//
void TMesh::add_nodes( const int ti,
                       const int numn )
{
  std::vector<double> bc(2);

  int sn = 3;
  int tot = (sn-1)*(sn-2)/2;
  while (numn > tot) {
    sn++;
    tot = (sn-1)*(sn-2)/2;
  }
  int k=0;
  for (int i=0;i<sn-2;i++) {
    bc[0] = (i+1.)/sn;
    for (int j=i;j<sn-2;j++) {
      bc[1] = (j-i+1.)/sn;
      add_node(ti,bc);
      k++;
      if (k >= numn) return;
    }
  }
}

//
// add a node to the mesh in the center of the specified triangle
// appending the location to the state vector
//
void TMesh::add_node_state( const int ti,
                            std::vector<double> & state )
{
  std::vector<double> bc(2);
  bc[0] = 1./3.;
  bc[1] = 1./3.;
  add_node(ti,bc);
  if (logvars == 1) {
    state.push_back(log(nodes.back().d.inputs[0]));
    state.push_back(log(nodes.back().d.inputs[1]));
  }
  else {
    state.push_back(nodes.back().d.inputs[0]);
    state.push_back(nodes.back().d.inputs[1]);
  }
}

//
// add a node to the mesh at the given location
//
void TMesh::add_node( double x,
                      double y )
{
  Node n(2,13);

  set_independent_vars(x,y,n);
  evaluate_node(n);

  nodes.push_back(n);
}

//
// add an empty node to the mesh -- must call move_node on it before
// any other computations
//
void TMesh::add_node_empty()
{
  nodes.push_back(Node(2,13));
}

//
// move a node to the given location
//
void TMesh::move_node( int node,
                       double x,
                       double y )
{
  set_independent_vars(x,y,nodes[node]);
  evaluate_node(nodes[node]);
}

//
// move a node toward the given location by the factor d
//
void TMesh::move_node_toward( int node,
                              double x,
                              double y,
                              double d )
{
  std::pair<double,double> c = get_independent_vars(nodes[node]);
  
  if (logvars) set_independent_vars(exp(log(x)*d+log(c.first)*(1.-d)),
                                    exp(log(y)*d+log(c.second)*(1.-d)),nodes[node]);
  else set_independent_vars(x*d+c.first*(1.-d),y*d+c.second*(1.-d),nodes[node]);
  evaluate_node(nodes[node]);
}

//
// add interior nodes to the mesh with locations given by consecutive
// values in the state array
//
void TMesh::set_interior_nodes( const std::vector<double> & state )
{
  if (state.size()%2 == 1) throw std::runtime_error("TMesh::set_interior_nodes: odd number of values in state vector");

  nodes.erase(nodes.begin()+nbnodes,nodes.end());

  for (std::size_t i=0;i<state.size();i+=2) {
    if (logvars) add_node(exp(state[i]),exp(state[i+1]));
    else add_node(state[i],state[i+1]);
  }
}

//
// evaluate a node
//
void TMesh::evaluate_node( Node & node )
{
  tm->evaluate_region(tm->regions[tmregion].phase,type,node.d);
}

//
// evaluate a data point
//
void TMesh::evaluate_point( EOSData & d )
{
  tm->evaluate_region(tm->regions[tmregion].phase,type,d);
}

//
// print the boundary
//
void TMesh::print_boundary( std::ostream & ofile )
{
  std::pair<double,double> xy;
  ofile << std::setprecision(15);
  for (int i=0;i<nbnodes+1;i++) {
    xy = get_independent_vars(nodes[i%nbnodes]);
    ofile << "bp " << i;
    if (logvars == 1) {
      ofile << " x " << log(xy.first) << " y " << log(xy.second) << std::endl;
    }
    else {
      ofile << " x " << xy.first << " y " << xy.second << std::endl;
    }
  }
  ofile << std::endl << std::endl;

}

//
// print a triangle
//
void TMesh::print_tri( int i,
                       std::ostream & ofile )
{
  Tri & t = tris[i];
  std::pair<double,double> xy[3];
  for (int i=0;i<3;i++) xy[i] = get_independent_vars(nodes[t.n[i]]);
  ofile << std::setprecision(15);
  if (logvars == 1) {
    ofile << "# n1 " << t.n[0] << " n2 " << t.n[1] << " n3 " << t.n[2] << std::endl;
    ofile << "x1 " << log(xy[0].first) << " y1 " << log(xy[0].second) << std::endl;
    ofile << "x2 " << log(xy[1].first) << " y2 " << log(xy[1].second) << std::endl;
    ofile << "x3 " << log(xy[2].first) << " y3 " << log(xy[2].second) << std::endl;
    ofile << "x1 " << log(xy[0].first) << " y1 " << log(xy[0].second) << std::endl;
    ofile << std::endl;
  }
  else {
    ofile << "# n1 " << t.n[0] << " n2 " << t.n[1] << " n3 " << t.n[2] << std::endl;
    ofile << "x1 " << xy[0].first << " y1 " << xy[0].second << std::endl;
    ofile << "x2 " << xy[1].first << " y2 " << xy[1].second << std::endl;
    ofile << "x3 " << xy[2].first << " y3 " << xy[2].second << std::endl;
    ofile << "x1 " << xy[0].first << " y1 " << xy[0].second << std::endl;
    ofile << std::endl;
  }
}

void TMesh::compute_barycentric_coords()
{
  double px[3],py[3];

  for (std::size_t i=0;i<tris.size();i++) {
    // compute barycentric coordinate helpers
    Tri & t = tris[i];
    if (logvars) {
      for (int j=0;j<3;j++) {
	px[j] = log(nodes[t.n[j]].d.inputs[0]);
	py[j] = log(nodes[t.n[j]].d.inputs[1]);
      }
    }
    else {
      for (int j=0;j<3;j++) {
	px[j] = nodes[t.n[j]].d.inputs[0];
	py[j] = nodes[t.n[j]].d.inputs[1];
      }
    }
    t.a[0] = px[0]-px[2];
    t.a[1] = px[1]-px[2];
    t.a[2] = px[2];
    t.a[3] = py[0]-py[2];
    t.a[4] = py[1]-py[2];
    t.a[5] = py[2];
    double invdet = 1./(t.a[0]*t.a[4]-t.a[1]*t.a[3]);
    t.a[0] *= invdet;
    t.a[1] *= invdet;
    t.a[3] *= invdet;
    t.a[4] *= invdet;
  }

}

//
// interpolation function
//
// Performs linear interpolation on all node quantities given the
// barycentric coordinates b1,b2 formed from the three nodes. b1 is
// assumed to correspond to node n1, b2 to n2, and b3=1-b1-b2 to n3.
// This assumes that all passed nodes have identically sized input and
// output vectors.
//
void TMesh::node_interp( double b1,
                         double b2,
                         Node & n1,
                         Node & n2,
                         Node & n3,
                         Node & result )
{
  if (logvars == 1) {
    // log scale the inputs
    for (std::size_t i=0;i<n1.d.inputs.size();i++)
      result.d.inputs[i] = exp(b1*log(n1.d.inputs[i])+b2*log(n2.d.inputs[i])+(1.-b1-b2)*log(n3.d.inputs[i]));
  }
  else {
    for (std::size_t i=0;i<n1.d.inputs.size();i++)
      result.d.inputs[i] = b1*n1.d.inputs[i]+b2*n2.d.inputs[i]+(1.-b1-b2)*n3.d.inputs[i];
  }

  for (std::size_t i=0;i<n1.d.outputs.size();i++)
    result.d.outputs[i] = b1*n1.d.outputs[i]+b2*n2.d.outputs[i]+(1.-b1-b2)*n3.d.outputs[i];

}

std::pair<double,double> TMesh::get_independent_vars( const Node & n )
{
  return std::pair<double,double>(n.d.inputs[0],n.d.inputs[1]);
}

//
// set the independent variable according to the type of boundary
//
// this version transfers value x into node n
//
void TMesh::set_independent_vars( double x,
                                  double y,
                                  Node & n )
{
  n.d.inputs[0] = x;
  n.d.inputs[1] = y;
}
//
// this version transfers values x and y in node n1 into node n2
//
void TMesh::set_independent_vars( Node & n1,
                                  Node & n2 )
{
  set_independent_vars(n1.d.inputs[0],n1.d.inputs[1],n2);
}

//
// get the square of the smallest angle in the triangle defined by the
// two specified nodes and point, with locations normalized by bounds
//
double TMesh::get_angle_sq( const int n1,
                            const int n2,
                            const double x,
                            const double y )
{
  int nn(nodes.size());
  if (n1 < 0 || n2 < 0 || n1 >= nn || n2 >= nn || n1 == n2)
    throw std::runtime_error("TMesh::get_angle: invalid node numbers");

  double px[3];
  double py[3];
  double norm[3];

  px[0] = x;
  py[0] = y;
  px[1] = nodes[n1].d.inputs[0];
  py[1] = nodes[n1].d.inputs[1];
  px[2] = nodes[n2].d.inputs[0];
  py[2] = nodes[n2].d.inputs[1];
  if (logvars) {
    px[1] = log(px[1]);
    py[1] = log(py[1]);
    px[2] = log(px[2]);
    py[2] = log(py[2]);
  }
  double dxi = 1./(xbounds[1]-xbounds[0]);
  double dyi = 1./(ybounds[1]-ybounds[0]);
  for (int i=0;i<3;i++) {
    px[i] *= dxi;
    py[i] *= dyi;
  }
  norm[0] = (px[2]-px[1])*(px[2]-px[1])+(py[2]-py[1])*(py[2]-py[1]);
  norm[1] = (px[0]-px[2])*(px[0]-px[2])+(py[0]-py[2])*(py[0]-py[2]);
  norm[2] = (px[1]-px[0])*(px[1]-px[0])+(py[1]-py[0])*(py[1]-py[0]);
  int minn = 0;
  if (norm[1] < norm[minn]) minn = 1;
  if (norm[2] < norm[minn]) minn = 2;

  int p0 = minn;
  int p1 = (minn+1)%3;
  int p2 = (minn+2)%3;

  double cp1 = (px[p1]-px[p0])*(py[p2]-py[p0]);
  double cp2 = (px[p2]-px[p0])*(py[p1]-py[p0]);
  return (cp1-cp2)*(cp1-cp2)/norm[p1]/norm[p2];
}
