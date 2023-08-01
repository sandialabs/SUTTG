/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SUTTG, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#ifndef BOUND_DELAUNAY_TRIANGULATION_H
#define BOUND_DELAUNAY_TRIANGULATION_H

//
// BoundDelaunayTriangulation.H:
//
// A templated class for Delaunay triangulation with a boundary
// constraint.
//

//
// Data precision:
//
// The class accepts double precision floats for the input location of
// points for triangulation. The class template variable "Real" is
// used for most internal floating point calculations. The template
// variable "RealA" is used for calculating the hull angles. RealA
// should generally be a larger precision than Real. For meshes that
// have very high aspect ratio triangles, specifing "Real" to use
// high-precision floating point class may result in correct
// triangulation of a mesh that otherwise would fail on a standard
// precision double variable.
//
// When built-in types are not used for RealA or Real, then some
// template specializations are required. In particular, one must
// specialize the "convertRealAtoReal" function for the specific class
// of interest. Also, the getEps function is used on Real and RealA to
// get the epsilon value. The default implementation uses the
// std::numeric_limits class. If this is not available for the desired
// classes, a specialization must be provided. For RealA, the getPI
// function calculates pi using an overload of the atan function. If
// this is not available, then the getPI function must also be
// specialized.
//

//
// algorithms:
//
// The sweep hull algorithm of David Sinclair is used for computing
// the Delaunay triangulation of the point set. The boundary
// conditions are imposed by flipping triangles to ensure neighboring
// boundary points reside on a triangle edge. Once such an edge is in
// place, the associated triangles may no longer be flipped across
// that edge, regardless of the criteria for Delaunay triangles. Once
// the constrained triangulation has been completed, triangles lying
// outside the boundary are pruned from the results.
//

#include <stdexcept> // for std::runtime_error
#include <vector>    // for std::vector
#include <list>      // for std::list
#include <stack>     // for std::stack
#include <cmath>     // for fabs
#include <limits>    // for std::numeric_limits

#include "BoundaryInfo.H" // for BoundaryInfo class

//
// convert from high precision number back to normal
//
template< class Real, class RealA >
Real convertRealAtoReal( const RealA & var )
{
  return var;
}

//
// internal helper structs for points and triangles
//
template < class Real >
struct Point 
{
  Point() : x(0.0), y(0.0), id(-1), rc2(0.0), theta(0.0), ti(-1) {}

  Real x,y;
  int id;
  Real rc2;
  Real theta;
  int ti;
};

//
// comparison function for testing equivalence of points
//
template < class Real >
bool PointCompare( Point<Real> p1,
                   Point<Real> p2 )
{
  if (p1.x < p2.x) return true;
  else if (p1.x > p2.x) return false;
  else {
    if (p1.y < p2.y) return true;
    else return false;
  }
}

//
// comparison function for point angles
//
template < class Real >
bool PointAngleCompare( Point<Real> p1,
                        Point<Real> p2 )
{
  return p1.theta < p2.theta;
}

//
// comparison function for point radii
//
template < class Real >
bool PointRadiusCompare( Point<Real> p1,
                         Point<Real> p2 )
{
  return p1.rc2 < p2.rc2;
}

template < class Real >
struct Triangle
{
  Point<Real> p[3];
  Point<Real> cc;
  int n[3];
  bool h[3];
  typename std::list<Point<Real> >::iterator hi[3];
};

//
// the main class
//

template < class Real, class RealA >
class BoundDelaunayTriangulation
{
public:

  //
  // Constructor
  //
  // The class constructor enters the input points with coordinates
  // (px,py) into the internal data structures, scaling the px and py
  // coordinates by xscale and yscale respectively. The size of px and
  // py must be identical. If nb is positive, then the first nb points
  // represent the boundary constraint, in anti-clockwise
  // ordering. The size of px and py must be at least nb.
  //
  // input:
  //
  //   px     -- vector of x coordinates for points
  //   py     -- vector of y coordinates for points
  //   nb     -- integer number of boundary points
  //   xscale -- scaling factor for x coordinates
  //   yscale -- scaling factor for y coordinates
  //
  BoundDelaunayTriangulation( const std::vector<double> & px,
                              const std::vector<double> & py,
                              const int nb,
                              const double xscale,
                              const double yscale );

  //
  // Destructor
  //
  BoundDelaunayTriangulation();

  //
  // triangulation method
  //
  void triangulate();

  //
  // method to obtain the resulting triangles
  //
  // The returned vector contains the node numbers for the triangles
  // as indexes to the original px and py vectors. The number of
  // triangles is 1/3 the size of the returned vector. Triangles nodes
  // are consecutive and in anti-clockwise order. For example, the
  // three nodes of the third triangle would be ordered (n[6],n[7],n[8]).
  //
  const std::vector<int> & getTriangles() const { return tris; }

private:

  //
  // the sweep hull algorithm (possibly recursive through boundary enforcement)
  //
  int SweepHullDelaunay( const std::vector<Point<Real> > & points,
                         const int nbp,
                         const std::vector<int> & bl,
                         std::vector<Triangle<Real> > & triangles,
                         const bool fixNeighbors = true );
  
  //
  // boundary enforcement routine for sweep hull algorithm
  //
  int enforceBoundary( const int id,
                       std::vector<Triangle<Real> > & t,
                       std::vector<int> & bplf,
                       const int nb,
                       const BoundaryInfo & binfo );

  //
  // recursively flip triangles to meet the Delaunay criteria
  //
  int DelaunayFlip( std::vector<Triangle<Real> > & t,
                    const int nb,
                    const BoundaryInfo & binfo,
                    const int tri = -1 );

  //
  // sort points by distance to a given point
  //
  void sortCenter( const Point<Real> & c,
                   std::vector<Point<Real> > & v );

  //
  // compute the circumcircle for a triangle
  // returns -1 on failure (points colinear at RealA precision)
  //
  int findCircumcircle( Triangle<Real> & t );

  //
  // Check if the given triangle has nodes (p1,p2,p3) oriented
  // anti-clockwise and if not flip the last two nodes
  //
  void setAntiTri( Triangle<Real> & t );

  //
  // calculate angle from x axis with a given point as origin
  //
  void calcCenterAngle( const Point<Real> & c,
                        Point<Real> & p );

  //
  // helper for debugging
  //
  void printTriangles( const std::vector<Triangle<Real> > & t);

  void ringIncrementHullIterator( std::list<Point<Real> > & hull,
                                  typename std::list<Point<Real> >::iterator & i,
                                  const bool up )
  {
    if (up) {
      i++;
      if (i == hull.end()) i = hull.begin();
    }
    else {
      if (i == hull.begin()) i = hull.end();
      i--;
    }
  }

  //
  // update neighbor information
  //
  void updateTriangleNeighbor( Triangle<Real> & t,
                               const int pid,
                               const int tnum )
  {
    if (t.p[0].id == pid) t.n[0] = tnum;
    else if (t.p[1].id == pid) t.n[1] = tnum;
    else t.n[2] = tnum;
  }

  //
  // number of boundary points
  //
  int numBP;

  //
  // nodes of the triangles
  //
  std::vector<int> tris;

  //
  // points for triangulation
  //
  std::vector<Point<Real> > points;

  //
  // debug flag for hull algorithm
  //
  int debug;

  //
  // debug flag for boundary enforcement
  //
  int ebdebug;

  //
  // area tolerance for triangles
  //
  Real areaCutoff;

  //
  // angle cutoff for new triangles
  //
  Real angleCutoff;

  //
  // high precision pi
  //
  RealA PI_RealA;
};

//
// get the epsilon value for a type
//
template < class Real >
Real getEps()
{
  return std::numeric_limits<Real>::epsilon();
}

//
// get pi for a type
//
template < class Real >
Real getPI()
{
  Real a(1.0);
  return 4.0*atan(a);
}

//
// Constructor
//
template < class Real, class RealA >
BoundDelaunayTriangulation<Real,RealA>
::BoundDelaunayTriangulation( const std::vector<double> & px,
                              const std::vector<double> & py,
                              const int nb,
                              const double xscale,
                              const double yscale )
  : numBP(nb), tris(0), points(px.size()), debug(0), ebdebug(0),
    areaCutoff(0.0), angleCutoff(0.0), PI_RealA(0.0)
{
  //
  // error checks for inputs
  //
  if (px.size() < 3)
    throw std::runtime_error("BoundDelaunayTriangulation: not enough points");
  if (px.size() != py.size())
    throw std::runtime_error("BoundDelaunayTriangulation: x,y sizes different");
  if (numBP < 0) numBP = 0;
  if (numBP > int(px.size()))
    throw std::runtime_error("BoundDelaunayTriangulation: too many boundary points specified");
  if (numBP > 0 && numBP < 3)
    throw std::runtime_error("BoundDelaunayTriangulation: must specify none or at least three boundary points");
  if (xscale <= 0.0 || yscale <= 0.0)
    throw std::runtime_error("BoundDelaunayTriangulation: scaling factors must be positive");

  //
  // test a valid (unique) set of points
  //
  std::vector<Point<double> > p(px.size());
  for (std::size_t i=0;i<p.size();i++) {
    p[i].x = px[i];
    p[i].y = py[i];
  }

  sort(p.begin(),p.end(),PointCompare<double>);

  for (std::size_t i=1;i<points.size();i++) {
    if (fabs(p[i].x-p[i-1].x) < 1.e-14*fabs(p[i].x+p[i-1].x) &&
        fabs(p[i].y-p[i-1].y) < 1.e-14*fabs(p[i].y+p[i-1].y)) {
      throw std::runtime_error("BoundDelaunayTriangulation: found non-unique set of points");
    }
  }

  //
  // fill points for triangulation
  //
  for (std::size_t i=0;i<points.size();i++) {
    points[i].id = i;
    points[i].x = px[i]/xscale;
    points[i].y = py[i]/yscale;
  }

  //
  // set tolerances
  //
  Real eps = getEps<Real>();
  Real epsa = convertRealAtoReal<Real,RealA>(getEps<RealA>());
  angleCutoff = 1.e5*epsa;
  areaCutoff = 2.0*eps;

  //
  // set pi
  //
  PI_RealA = getPI<RealA>();
}

template < class Real, class RealA >
void BoundDelaunayTriangulation<Real,RealA>
::sortCenter( const Point<Real> & c,
              std::vector<Point<Real> > & v )
{
  //
  // compute new radii^2 to the center point c
  //
  for (typename std::vector<Point<Real> >::iterator p=v.begin();p<v.end();p++) {
    p->rc2 = (p->x-c.x)*(p->x-c.x)+(p->y-c.y)*(p->y-c.y);
  }

  std::sort(v.begin(),v.end(),PointRadiusCompare<Real>);

}

template < class Real, class RealA >
int BoundDelaunayTriangulation<Real,RealA>
::findCircumcircle( Triangle<Real> & t )
{
  Real d;
  Real bx,by,cx,cy;

  bx = t.p[1].x-t.p[0].x;
  by = t.p[1].y-t.p[0].y;
  cx = t.p[2].x-t.p[0].x;
  cy = t.p[2].y-t.p[0].y;

  d = 2.*(bx*cy-by*cx);
  if (fabs(d) < areaCutoff) {
    return -1;
  }

  // center of circumcircle
  t.cc.x = t.p[0].x+(cy*(bx*bx+by*by)-by*(cx*cx+cy*cy))/d;
  t.cc.y = t.p[0].y+(bx*(cx*cx+cy*cy)-cx*(bx*bx+by*by))/d;

  // radius squared of circumcircle
  t.cc.rc2 = (t.cc.x-t.p[0].x)*(t.cc.x-t.p[0].x)+(t.cc.y-t.p[0].y)*(t.cc.y-t.p[0].y);
  //std::cout << " cc x " << t.cc.x << " y  " << t.cc.y << " r2 " << t.cc.rc2 << std::endl;

  return 0;

}

template < class Real, class RealA >
void BoundDelaunayTriangulation<Real,RealA>
::setAntiTri( Triangle<Real> & t )
{
  // if anti-clockwise following determinant is positive
  // ignore degenerate case

  Real x1 = t.p[1].x-t.p[0].x;
  Real x2 = t.p[2].x-t.p[0].x;
  Real y1 = t.p[1].y-t.p[0].y;
  Real y2 = t.p[2].y-t.p[0].y;
  Real det = x1*y2-x2*y1;

  if (det < 0.0) std::swap(t.p[1],t.p[2]);

}

template < class Real, class RealA >
void BoundDelaunayTriangulation<Real,RealA>
::calcCenterAngle( const Point<Real> & c,
                   Point<Real> & p )
{
  RealA x = p.x-c.x;
  RealA y = p.y-c.y;

  RealA theta = atan2(y,x)+PI_RealA;
  p.theta = convertRealAtoReal<Real,RealA>(theta);
}

template < class Real, class RealA >
void BoundDelaunayTriangulation<Real,RealA>
::printTriangles( const std::vector<Triangle<Real> > & t)
{
  for (std::size_t i=0;i<t.size();i++) {
    std::cout << "# t " << i << std::endl;
    std::cout << "n1 " << t[i].p[0].x << " " << t[i].p[0].y << " id " << t[i].p[0].id << " nt " << t[i].n[0] << std::endl
	      << "n1 " << t[i].p[1].x << " " << t[i].p[1].y << " id " << t[i].p[1].id << " nt " << t[i].n[1] << std::endl
	      << "n1 " << t[i].p[2].x << " " << t[i].p[2].y << " id " << t[i].p[2].id << " nt " << t[i].n[2] << std::endl
	      << "n1 " << t[i].p[0].x << " " << t[i].p[0].y << " id " << t[i].p[0].id << " nt " << t[i].n[0] << std::endl
	      << std::endl;
  }
  std::cout << std::endl;
}

template < class Real, class RealA >
void BoundDelaunayTriangulation<Real,RealA>
::triangulate()
{
  //
  // generate boundary point list and boundary info
  //
  std::vector<int> bl(numBP);
  for (std::size_t i=0;i<bl.size();i++) bl[i] = i;

  //
  // hold the resulting triangles here
  //
  std::vector<Triangle<Real> > t;

  //
  // actually do the triangulation
  //
  if (SweepHullDelaunay(points,numBP,bl,t,false) != 0)
    throw std::runtime_error("BoundDelaunayTriangulation::triangulate: SweepHullDelaunay failed");

  //
  // convert triangles to tris vector, putting smallest node first
  //
  for (std::size_t i=0;i<t.size();i++) {
    int j = 0;
    if (t[i].p[1].id < t[i].p[j].id) j = 1;
    if (t[i].p[2].id < t[i].p[j].id) j = 2;
    if (j==0) {
      tris.push_back(t[i].p[0].id);
      tris.push_back(t[i].p[1].id);
      tris.push_back(t[i].p[2].id);
    }
    else if (j==1) {
      tris.push_back(t[i].p[1].id);
      tris.push_back(t[i].p[2].id);
      tris.push_back(t[i].p[0].id);
    }
    else {
      tris.push_back(t[i].p[2].id);
      tris.push_back(t[i].p[0].id);
      tris.push_back(t[i].p[1].id);
    }
  }
}

template < class Real, class RealA >
int BoundDelaunayTriangulation<Real,RealA>
::SweepHullDelaunay( const std::vector<Point<Real> > & points,
                     const int nbp,
                     const std::vector<int> & bl,
                     std::vector<Triangle<Real> > & triangles,
                     const bool fixNeighbors )
{
  if (debug > 0) std::cout << "# Starting SweepHullDelaunay" << std::endl;

  //
  // generate boundary info
  //
  BoundaryInfo binfo(nbp,bl);

  //
  // setup points
  //
  std::vector<Point<Real> > pleft = points;
  std::vector<int> bplflags(nbp,1); // 1 if boundary point still left

  Point<Real> center;
  Triangle<Real> besttri;
  Triangle<Real> curtri;

  //
  // select first point as seed
  //
  center = pleft[0];
  curtri.p[0] = center;
  if (center.id < nbp) bplflags[center.id] = 0;
  pleft.erase(pleft.begin());

  //
  // second point is closest to first
  //
  sortCenter(center,pleft);
  curtri.p[1] = pleft[0];
  if (curtri.p[1].id < nbp) bplflags[curtri.p[1].id] = 0;
  pleft.erase(pleft.begin());
 
  // find smallest circum-circle with these two points
  typename std::vector<Point<Real> >::iterator bestp = pleft.begin();
  besttri = curtri;
  besttri.cc.rc2 = 9.e99;
  for (typename std::vector<Point<Real> >::iterator p=pleft.begin();p<pleft.end();p++) {
    curtri.p[2] = *p;
    if (debug > 0) {
      std::cout << " p1 " << curtri.p[0].x << " " << curtri.p[0].y
		<< " p2 " << curtri.p[1].x << " " << curtri.p[1].y
		<< " p3 " << curtri.p[2].x << " " << curtri.p[2].y << std::endl;
    }
    if (findCircumcircle(curtri) != 0) continue;
    
    // update best
    if (curtri.cc.rc2 < besttri.cc.rc2) {
      besttri = curtri;
      bestp = p;
    }

    // stop if best cc diameter smaller than distance to point
    if (debug > 0) std::cout << " cc x " << curtri.cc.x << " y " << curtri.cc.y << " p->rc2 " << p->rc2 << " 4*tri.cc.rc2 " << 4.*besttri.cc.rc2 << std::endl;
    if (p->rc2 > 4.*besttri.cc.rc2) {
      break;
    }
  }
  if (besttri.cc.rc2 > 8.e99) {
    std::cout << "error in SweepHullDelaunay: failed to find first circumcircle" << std::endl;
    return -1;
  }
  // remove point
  if (bestp->id < nbp) bplflags[bestp->id] = 0;
  pleft.erase(bestp);

  setAntiTri(besttri);
  besttri.n[0] = -1;
  besttri.n[1] = -1;
  besttri.n[2] = -1;

  // set circumcenter
  Point<Real> ccenter = besttri.cc;
  
  // set center to the barycenter of the triangle -- avoids the problem of when 
  // the circumcenter is outside the triangle
  center = besttri.cc;
  center.x = (besttri.p[0].x+besttri.p[1].x+besttri.p[2].x)/3.;
  center.y = (besttri.p[0].y+besttri.p[1].y+besttri.p[2].y)/3.;
  
  // calculate angles to x axis at circumcenter
  calcCenterAngle(center,besttri.p[0]);
  calcCenterAngle(center,besttri.p[1]);
  calcCenterAngle(center,besttri.p[2]);
  for (typename std::vector<Point<Real> >::iterator p=pleft.begin();p<pleft.end();p++) {
    calcCenterAngle(center,*p);
  }

  if (debug > 0) std::cout << "center x " << center.x << " y " << center.y << std::endl;

  // initialize hull
  std::list<Point<Real> > hull;

  besttri.p[0].ti = 0;
  besttri.hi[0] = hull.insert(hull.begin(),besttri.p[0]);
  besttri.h[0] = true;
  besttri.p[1].ti = 0;
  besttri.hi[1] = hull.insert(std::lower_bound(hull.begin(),hull.end(),besttri.p[1],PointAngleCompare<Real>),besttri.p[1]);
  besttri.h[1] = true;
  besttri.p[2].ti = 0;
  besttri.hi[2] = hull.insert(std::lower_bound(hull.begin(),hull.end(),besttri.p[2],PointAngleCompare<Real>),besttri.p[2]);
  besttri.h[2] = true;
  if (debug > 0) {
    std::cout << "# hull points" << std::endl;
    for (typename std::list<Point<Real> >::iterator i=hull.begin();i!=hull.end();i++) {
      std::cout << "x " << i->x << " y " << i->y << " t " << i->ti << std::endl;
    }
    std::cout << std::endl << std::endl << std::endl;
  }
  // put in first triangle
  std::vector<Triangle<Real> > newtris;
  newtris.push_back(besttri);
  
  if (debug > 0) {
    for (std::size_t i=0;i<newtris.size();i++) {
      std::cout << "# i " << i << " n1 " << newtris[i].n[0] << " n2 " << newtris[i].n[1] << " n3 " << newtris[i].n[2] << std::endl
		<< "x " << newtris[i].p[0].x << " y " << newtris[i].p[0].y << std::endl
		<< "x " << newtris[i].p[1].x << " y " << newtris[i].p[1].y << std::endl
		<< "x " << newtris[i].p[2].x << " y " << newtris[i].p[2].y << std::endl
		<< "x " << newtris[i].p[0].x << " y " << newtris[i].p[0].y << std::endl
		<< std::endl ;//<< std::endl;
    }
    std::cout << std::endl;
  }
  // resort remaining points according to circumcenter
  sortCenter(ccenter,pleft);

  // make new triangles
  while(pleft.size() > 0) {
    Point<Real> nextp = pleft[0];
    
    if (nextp.id < nbp) bplflags[nextp.id] = 0;
    pleft.erase(pleft.begin());

    if (debug > 0) {
      std::cout << "# hull points a" << std::endl;
      for (typename std::list<Point<Real> >::iterator i=hull.begin();i!=hull.end();i++) {
	std::cout << "x " << i->x << " y " << i->y << " t " << i->ti << " th " << i->theta << std::endl;
      }
      std::cout << std::endl << std::endl << std::endl;
    }
    // put new point into hull
    typename std::list<Point<Real> >::iterator ihullc = hull.insert(std::lower_bound(hull.begin(),hull.end(),nextp,PointAngleCompare<Real>),nextp);
    if (debug > 0) {
      std::cout << "# hull points b" << std::endl;
      for (typename std::list<Point<Real> >::iterator i=hull.begin();i!=hull.end();i++) {
	std::cout << "x " << i->x << " y " << i->y << " t " << i->ti << " th " << i->theta << std::endl;
      }
      std::cout << std::endl << std::endl << std::endl;
    }
    typename std::list<Point<Real> >::iterator ihull1,ihull2;

    // make required triangle
    curtri.p[0] = *ihullc;
    ihull2 = ihullc;
    ringIncrementHullIterator(hull,ihull2,true);
    curtri.p[1] = *ihull2;
    ihull1 = ihullc;
    ringIncrementHullIterator(hull,ihull1,false);
    curtri.p[2] = *ihull1;
    // don't add triangle if outside boundary
    if (debug > 0) {
      std::cout << "# a p1 " << curtri.p[0].id << " p2 " << curtri.p[1].id << " p3 " << curtri.p[2].id << std::endl;
      std::cout << "# a added" << std::endl;
    }
    curtri.n[0] = -1;
    curtri.n[1] = ihull1->ti;
    curtri.n[2] = -1;
    ihullc->ti = newtris.size();
    for (int i=0;i<3;i++) if (newtris[ihull1->ti].hi[i] == ihull1) newtris[ihull1->ti].h[i] = false;
    ihull1->ti = newtris.size();
    curtri.h[0] = true;
    curtri.h[1] = false;
    curtri.h[2] = true;
    curtri.hi[0] = ihullc;
    curtri.hi[2] = ihull1;
    updateTriangleNeighbor(newtris[curtri.n[1]],curtri.p[2].id,newtris.size());
    if (findCircumcircle(curtri) != 0) {
      std::cout << "error in sweep_hull_delaunay: failed to find required tri circumcircle -- 1" << std::endl;
      std::cout << "# n1 " << curtri.n[0] << " n2 " << curtri.n[1] << " n3 " << curtri.n[2] << std::endl
		<< "x " << curtri.p[0].x << " y " << curtri.p[0].y << std::endl
		<< "x " << curtri.p[1].x << " y " << curtri.p[1].y << std::endl
		<< "x " << curtri.p[2].x << " y " << curtri.p[2].y << std::endl
		<< "x " << curtri.p[0].x << " y " << curtri.p[0].y << std::endl
		<< std::endl << std::endl;
      return -1;
    }
    if (debug > 0) std::cout << " curtri cc x " << curtri.cc.x << " y " << curtri.cc.y << std::endl;
    newtris.push_back(curtri);
    if (DelaunayFlip(newtris,nbp,binfo) != 0) {
      std::cout << "error in sweep_hull_delaunay: failed to flip triangles" << std::endl;
      return -1;
    }

    // make new triangles, removing points from hull
    curtri.p[0] = *ihullc;
    ihull2 = ihullc;
    ringIncrementHullIterator(hull,ihull2,true);
    // increasing angle direction
    while (1) {
      ihull1 = ihull2;
      ringIncrementHullIterator(hull,ihull2,true);

      curtri.p[1] = *ihull2;
      curtri.p[2] = *ihull1;

      // add triangle if new angle is positive
      Real cp1 = (curtri.p[1].x-curtri.p[0].x)*(curtri.p[2].y-curtri.p[0].y);
      Real cp2 = (curtri.p[2].x-curtri.p[0].x)*(curtri.p[1].y-curtri.p[0].y);
      Real norm1 = (curtri.p[1].x-curtri.p[0].x)*(curtri.p[1].x-curtri.p[0].x)+(curtri.p[1].y-curtri.p[0].y)*(curtri.p[1].y-curtri.p[0].y);
      Real norm2 = (curtri.p[2].x-curtri.p[0].x)*(curtri.p[2].x-curtri.p[0].x)+(curtri.p[2].y-curtri.p[0].y)*(curtri.p[2].y-curtri.p[0].y);
      if (debug > 0) {
        std::cout << "# increasing angle" << std::endl;
        std::cout << "# p1 x " << curtri.p[0].x << " y " << curtri.p[0].y << std::endl;
        std::cout << "# p2 x " << curtri.p[1].x << " y " << curtri.p[1].y << std::endl;
        std::cout << "# p3 x " << curtri.p[2].x << " y " << curtri.p[2].y << std::endl;
        std::cout << "# sin(o)^2 = " << (cp1-cp2)*(cp1-cp2)/norm1/norm2 << std::endl;
        std::cout << "# cp1 " << cp1 << " cp2 " << cp2 << " fabs(cp1-cp2) " << fabs(cp1-cp2) << " 1.e-14*fabs(cp1+cp2) " << 1.e-14*fabs(cp1+cp2) << std::endl;
      }
      if (cp1 > cp2 && fabs(cp1-cp2) > 1.e-14*fabs(cp1+cp2) && findCircumcircle(curtri) == 0 && (cp1-cp2)*(cp1-cp2)/norm1/norm2 > angleCutoff) {
	if (debug > 0) std::cout << "# b added" << std::endl;
	curtri.n[0] = -1;
	curtri.n[1] = ihull1->ti;
	curtri.n[2] = ihullc->ti;

	for (int i=0;i<3;i++) if (newtris[ihull1->ti].hi[i] == ihull1) newtris[ihull1->ti].h[i] = false;
	for (int i=0;i<3;i++) if (newtris[ihullc->ti].hi[i] == ihullc) newtris[ihullc->ti].h[i] = false;
	ihullc->ti = newtris.size();
	curtri.h[0] = true;
	curtri.h[1] = false;
	curtri.h[2] = false;
	curtri.hi[0] = ihullc;
	updateTriangleNeighbor(newtris[curtri.n[1]],curtri.p[2].id,newtris.size());
	updateTriangleNeighbor(newtris[curtri.n[2]],curtri.p[0].id,newtris.size());
	if (findCircumcircle(curtri) != 0) {
	  std::cout << "error in sweep_hull_delaunay: failed to find required tri circumcircle -- 2" << std::endl;
	  std::cout << "# n1 " << curtri.n[0] << " n2 " << curtri.n[1] << " n3 " << curtri.n[2] << std::endl
		    << "x " << curtri.p[0].x << " y " << curtri.p[0].y << std::endl
		    << "x " << curtri.p[1].x << " y " << curtri.p[1].y << std::endl
		    << "x " << curtri.p[2].x << " y " << curtri.p[2].y << std::endl
		    << "x " << curtri.p[0].x << " y " << curtri.p[0].y << std::endl
		    << std::endl << std::endl;
	  return -1;
	}
	newtris.push_back(curtri);
	if (DelaunayFlip(newtris,nbp,binfo) != 0) {
	  std::cout << "error in sweep_hull_delaunay: failed to flip triangles" << std::endl;
	  return -1;
	}
	
	// remove hull point
	hull.erase(ihull1);
      }
      else break;
    }
    // decreasing angle direction
    ihull2 = ihullc;
    ringIncrementHullIterator(hull,ihull2,false);
    while(1) {
      ihull1 = ihull2;
      ringIncrementHullIterator(hull,ihull2,false);

      curtri.p[1] = *ihull1;
      curtri.p[2] = *ihull2;

      // add triangle if new angle is positive
      Real cp1 = (curtri.p[1].x-curtri.p[0].x)*(curtri.p[2].y-curtri.p[0].y);
      Real cp2 = (curtri.p[2].x-curtri.p[0].x)*(curtri.p[1].y-curtri.p[0].y);
      Real norm1 = (curtri.p[1].x-curtri.p[0].x)*(curtri.p[1].x-curtri.p[0].x)+(curtri.p[1].y-curtri.p[0].y)*(curtri.p[1].y-curtri.p[0].y);
      Real norm2 = (curtri.p[2].x-curtri.p[0].x)*(curtri.p[2].x-curtri.p[0].x)+(curtri.p[2].y-curtri.p[0].y)*(curtri.p[2].y-curtri.p[0].y);
      if (debug > 0) {
        std::cout << "# decreasing angle" << std::endl;
        std::cout << "# p1 x " << curtri.p[0].x << " y " << curtri.p[0].y << std::endl;
        std::cout << "# p2 x " << curtri.p[1].x << " y " << curtri.p[1].y << std::endl;
        std::cout << "# p3 x " << curtri.p[2].x << " y " << curtri.p[2].y << std::endl;
        std::cout << "# sin(o)^2 = " << (cp1-cp2)*(cp1-cp2)/norm1/norm2 << std::endl;
        std::cout << "# cp1 " << cp1 << " cp2 " << cp2 << " fabs(cp1-cp2) " << fabs(cp1-cp2) << " 1.e-14*fabs(cp1+cp2) " << 1.e-14*fabs(cp1+cp2) << std::endl;
      }
      if (cp1 > cp2 && fabs(cp1-cp2) > 1.e-14*fabs(cp1+cp2) && findCircumcircle(curtri) == 0 && (cp1-cp2)*(cp1-cp2)/norm1/norm2 > angleCutoff) {
	if (debug > 0) std::cout << "# c added" << std::endl;
	curtri.n[0] = ihull1->ti;
	curtri.n[1] = ihull2->ti;
	curtri.n[2] = -1;
	for (int i=0;i<3;i++) if (newtris[ihull1->ti].hi[i] == ihull1) newtris[ihull1->ti].h[i] = false;
	for (int i=0;i<3;i++) if (newtris[ihull2->ti].hi[i] == ihull2) newtris[ihull2->ti].h[i] = false;
	ihull2->ti = newtris.size();
	curtri.h[0] = false;
	curtri.h[1] = false;
	curtri.h[2] = true;
	curtri.hi[2] = ihull2;
	updateTriangleNeighbor(newtris[curtri.n[0]],curtri.p[1].id,newtris.size());
	updateTriangleNeighbor(newtris[curtri.n[1]],curtri.p[2].id,newtris.size());
	if (findCircumcircle(curtri) != 0) {
	  std::cout << "error in sweep_hull_delaunay: failed to find required tri circumcircle -- 3" << std::endl;
	  std::cout << "# n1 " << curtri.n[0] << " n2 " << curtri.n[1] << " n3 " << curtri.n[2] << std::endl
		    << "x " << curtri.p[0].x << " y " << curtri.p[0].y << std::endl
		    << "x " << curtri.p[1].x << " y " << curtri.p[1].y << std::endl
		    << "x " << curtri.p[2].x << " y " << curtri.p[2].y << std::endl
		    << "x " << curtri.p[0].x << " y " << curtri.p[0].y << std::endl
		    << std::endl << std::endl;
	  return -1;
	}
	newtris.push_back(curtri);
	if (DelaunayFlip(newtris,nbp,binfo) != 0) {
	  std::cout << "error in sweep_hull_delaunay: failed to flip triangles" << std::endl;
	  return -1;
	}
	
	// remove hull point
	hull.erase(ihull1);
      }
      else break;
    }
    if (debug > 0) {
      std::cout << "# hull points" << std::endl;
      for (typename std::list<Point<Real> >::iterator i=hull.begin();i!=hull.end();i++) {
	std::cout << "x " << i->x << " y " << i->y << " t " << i->ti << std::endl;
      }
      std::cout << std::endl << std::endl << std::endl;
    }

    //
    // check that boundaries properly connected
    //
    if (nextp.id < nbp) {
      bool done = false;
      while (!done) {
        if (enforceBoundary(nextp.id,newtris,bplflags,nbp,binfo) < 1) done = true;
      }
    }

    if (debug > 0) {
      for (std::size_t i=0;i<newtris.size();i++) {
	std::cout << "# i " << i << " n1 " << newtris[i].n[0] << " n2 " << newtris[i].n[1] << " n3 " << newtris[i].n[2] << std::endl
		  << "x " << newtris[i].p[0].x << " y " << newtris[i].p[0].y << std::endl
		  << "x " << newtris[i].p[1].x << " y " << newtris[i].p[1].y << std::endl
		  << "x " << newtris[i].p[2].x << " y " << newtris[i].p[2].y << std::endl
		  << "x " << newtris[i].p[0].x << " y " << newtris[i].p[0].y << std::endl
		  << std::endl ;//<< std::endl;
      }
      std::cout << std::endl;
    }
  }

  //
  // only keep triangles inside the boundary
  //
  std::vector<int> ti;
  for (std::size_t i=0;i<newtris.size();i++) {
    //if (tri_inside(nbp,newtris[i])) {
    if (binfo.isInside(newtris[i].p[0].id,newtris[i].p[1].id,newtris[i].p[2].id)) {
      triangles.push_back(newtris[i]);
      ti.push_back(i);
    }
  }

  if (fixNeighbors) {
    //
    // fix the neighbors
    //
    int last = -1;
    int size = triangles.size();
    for (int i=0;i<size;i++) {
      for (int j=0;j<size;j++) {
        for (int k=0;k<3;k++) {
          if (triangles[j].n[k] == ti[i]) triangles[j].n[k] = i;
          else if (triangles[j].n[k] < ti[i] && triangles[j].n[k] > last) triangles[j].n[k] = -1;
        }
      }
      last = ti[i];
    }
    for (std::size_t j=0;j<triangles.size();j++) {
      for (int k=0;k<3;k++) {
        if (triangles[j].n[k] > ti.back()) triangles[j].n[k] = -1;
      }
    }
  }

  return 0;

}

template < class Real >
struct enforceBoundaryData
{
  enforceBoundaryData( const int pointid,
                       const int neighbor,
                       const bool onhull,
                       const typename std::list<Point<Real> >::iterator hulliterator )
    : p(pointid), n(neighbor), h(onhull), hi(hulliterator) {}

  enforceBoundaryData() : p(-1), n(-1), h(false), hi(0) {}

  int p;
  int n;
  bool h;
  typename std::list<Point<Real> >::iterator hi;
};

template < class Real >
bool enforceBoundaryDataCompare( const enforceBoundaryData<Real> & a,
                                 const enforceBoundaryData<Real> & b )
{
  return a.p < b.p;
}

template < class Real, class RealA >
int BoundDelaunayTriangulation<Real,RealA>
::enforceBoundary( const int id,
                   std::vector<Triangle<Real> > & t,
                   std::vector<int> & bplf,
                   const int nb,
                   const BoundaryInfo & binfo )
{
  //std::cout.flush();

  // get the ids for neighboring boundary points
  std::vector<int> neighbors;
  binfo.getNeighbors(id,neighbors);
  std::vector<int> ids;
  ids.push_back(id);
  for (std::size_t i=0;i<neighbors.size();i++)
    if (bplf[neighbors[i]] == 0) ids.push_back(neighbors[i]);

  // neighboring points not yet in mesh
  if (ids.size() < 2) return 0;

  // get the triangles for each of the ids
  std::vector<std::vector<int> > itris(ids.size());
  for (std::size_t i=0;i<t.size();i++) {
    for (int j=0;j<3;j++) {
      for (std::size_t k=0;k<ids.size();k++) {
	if (t[i].p[j].id == ids[k]) itris[k].push_back(i);
      }
    }
  }
  
  // sort the tri list for easy checking
  for (std::size_t i=0;i<ids.size();i++) std::sort(itris[i].begin(),itris[i].end());

  // check if neighboring boundary points fall in a triangle
  std::vector<int> connected(ids.size()-1,0);
  for (std::size_t i=1;i<=connected.size();i++) {
    std::size_t k=0;
    for (std::size_t j=0;j<itris[i].size();j++) {
      // compare this points tris with its neighbor
      while (itris[i][j] > itris[0][k]) {
	k++;
	if (k == itris[0].size()) break;
      }
      if (k == itris[0].size()) break;
      else if (itris[i][j] == itris[0][k]) {
	connected[i-1] = 1;
	break;
      }
    }
  }

  // done if all connected
  bool done = true;
  for (std::size_t i=0;i<connected.size();i++) if (connected[i] == 0) done = false;
  if (done) return 0;

  if (ebdebug > 0) {
    std::cout << "# startt" << std::endl;
    printTriangles(t);
  }

  // vector for start and target boundary points
  std::vector<Point<Real> > pst;
  // starting point (new boundary node)
  pst.push_back(t[itris[0][0]].p[0]);
  for (int i=1;i<3;i++) if (t[itris[0][0]].p[i].id == ids[0]) pst[0] = t[itris[0][0]].p[i];
  // end point(s)
  for (std::size_t i=0;i<connected.size();i++) {
    if (connected[i] == 0) {
      pst.push_back(t[itris[i+1][0]].p[0]);
      for (int j=1;j<3;j++) if (t[itris[i+1][0]].p[j].id == ids[i+1]) pst.back() = t[itris[i+1][0]].p[j];
    }
  }
  if (ebdebug > 0) {
    std::cout << "# pst";
    for (std::size_t i=0;i<pst.size();i++) std::cout << " " << pst[i].id;
    std::cout << std::endl;
  }

  // convert target points to location relative to start
  for (std::size_t i=1;i<pst.size();i++) {
    pst[i].x -= pst[0].x;
    pst[i].y -= pst[0].y;
  }

  // find triangle and node paths from new boundary node to others
  for (std::size_t p=1;p<2/*pst.size()*/;p++) {
    // start by finding triangle at starting node
    std::vector<int> tripath;
    std::vector<int> prevneig;
    std::vector<Point<Real> > n1path,n2path;
    std::vector<enforceBoundaryData<Real> > bd;
    n1path.push_back(pst[0]);
    n2path.push_back(pst[0]);
    for (std::size_t i=0;i<itris[0].size();i++) {
      Triangle<Real> & ct = t[itris[0][i]];
      int ti = 0;
      for (int j=1;j<3;j++) if (ct.p[j].id == pst[0].id) ti = j;
      int ti1 = (ti+1)%3;
      int ti2 = (ti+2)%3;
      Real angle1 = (ct.p[ti1].x-pst[0].x)*pst[p].y-(ct.p[ti1].y-pst[0].y)*pst[p].x;
      Real angle2 = (ct.p[ti2].x-pst[0].x)*pst[p].y-(ct.p[ti2].y-pst[0].y)*pst[p].x;
      if (angle1*angle2 < 0.) {
        tripath.push_back(itris[0][i]);
        tripath.push_back(ct.n[ti1]);
	n1path.push_back(ct.p[ti1]);
	n2path.push_back(ct.p[ti2]);
	bd.push_back(enforceBoundaryData<Real>(ct.p[ti].id,ct.n[ti],ct.h[ti],ct.hi[ti]));
	bd.push_back(enforceBoundaryData<Real>(ct.p[ti2].id,ct.n[ti2],ct.h[ti2],ct.hi[ti2]));
        break;
      }
    }
    if (tripath.size() == 0) throw std::runtime_error("enforce_boundary: cannot enforce constraint -- no triangles between boundary points");
    if (ebdebug > 0) std::cout << "# pst " << p << " starttri " << tripath[0] << " nexttri " << tripath[1] << std::endl;
    while (1) {
      if (tripath.back() < 0) throw std::runtime_error("enforce_boundary: could not find triangle in path");
      if (tripath.size() > t.size()) throw std::runtime_error("enforce_boundary: more triangles in path than in mesh");

      // check if target node is in the last triangle
      Triangle<Real> & ct = t[tripath.back()];
      bool done = false;
      for (int i=0;i<3;i++) if (ct.p[i].id == pst[p].id) {
        if (ebdebug > 0) std::cout << "# pst end tri " << tripath.back() << std::endl;
	done = true;
	prevneig.push_back((i+1)%3);
	n1path.push_back(ct.p[i]);
	n2path.push_back(ct.p[i]);
	int j = (i+2)%3;
	bd.push_back(enforceBoundaryData<Real>(ct.p[j].id,ct.n[j],ct.h[j],ct.hi[j]));
	bd.push_back(enforceBoundaryData<Real>(ct.p[i].id,ct.n[i],ct.h[i],ct.hi[i]));
      }
      if (done) break;

      // find the next triangle
      Real angles[3];
      for (int i=0;i<3;i++) angles[i] = (ct.p[i].x-pst[0].x)*pst[p].y-(ct.p[i].y-pst[0].y)*pst[p].x;

      int oddnode = 0;
      for (int i=1;i<3;i++) if (angles[(i+1)%3]*angles[(i+2)%3] > 0.) oddnode = i;
      if (ebdebug > 0) std::cout << "# pst " << tripath.back() << " oddnode " << oddnode << std::endl;
      if (ct.n[oddnode] == *(tripath.rbegin()+1)) {
        if (ebdebug > 0) std::cout << "# pst " << tripath.back() << " 1 nexttri " << ct.n[(oddnode+2)%3] << std::endl;
	tripath.push_back(ct.n[(oddnode+2)%3]);
	prevneig.push_back(oddnode);
	n1path.push_back(ct.p[(oddnode+2)%3]);
	int j=(oddnode+1)%3;
	bd.push_back(enforceBoundaryData<Real>(ct.p[j].id,ct.n[j],ct.h[j],ct.hi[j]));
      }
      else {
        if (ebdebug > 0) std::cout << "# pst " << tripath.back() << " 2 nexttri " << ct.n[oddnode] << std::endl;
	tripath.push_back(ct.n[oddnode]);
	prevneig.push_back((oddnode+2)%3);
	n2path.push_back(ct.p[(oddnode+1)%3]);
	int j=(oddnode+1)%3;
	bd.push_back(enforceBoundaryData<Real>(ct.p[j].id,ct.n[j],ct.h[j],ct.hi[j]));
      }
    }
    std::sort(bd.begin(),bd.end(),enforceBoundaryDataCompare<Real>);
    if (ebdebug > 0) {
      std::cout << "# tripath";
      for (std::size_t i=0;i<tripath.size();i++) std::cout << " " << tripath[i];
      std::cout << std::endl;
      std::cout << "# pn";
      for (std::size_t i=0;i<prevneig.size();i++) std::cout << " " << prevneig[i] << "(" << t[tripath[i+1]].n[prevneig[i]] << ")";
      std::cout << std::endl;
      std::cout << "# n1p";
      for (std::size_t i=0;i<n1path.size();i++) std::cout << " " << n1path[i].id;
      std::cout << std::endl;
      std::cout << "# n2p";
      for (std::size_t i=0;i<n2path.size();i++) std::cout << " " << n2path[i].id;
      std::cout << std::endl;
      std::cout << "# bdp";
      for (std::size_t i=0;i<bd.size();i++) std::cout << " " << bd[i].p;
      std::cout << std::endl;
    }

    //
    // triangulate each region
    //
    std::vector<Triangle<Real> > nt1,nt2;

    // reverse path 2 to ensure anti-clockwise ordering
    std::reverse(n2path.begin(),n2path.end());
    //std::cout << "# n2p";
    //for (std::size_t i=0;i<n2path.size();i++) std::cout << " " << n2path[i].id;
    //std::cout << std::endl;

    //
    // remove any duplicate points, keeping them in the boundary
    //
    std::vector<Point<Real> > r1path,r2path;
    std::vector<int> r1bp,r2bp;
    int lastp = 0;
    for (std::size_t i=0;i<n1path.size();i++) {
      int newid = n1path[i].id;
      int found = -1;
      for (std::size_t j=0;j<r1path.size();j++) {
        if (r1path[j].id == newid) found = j;
      }
      if (found < 0) {
        r1path.push_back(n1path[i]);
        r1bp.push_back(lastp++);
      }
      else {
        r1bp.push_back(found);
      }
    }
    lastp = 0;
    for (std::size_t i=0;i<n2path.size();i++) {
      int newid = n2path[i].id;
      int found = -1;
      for (std::size_t j=0;j<r2path.size();j++) {
        if (r2path[j].id == newid) found = j;
      }
      if (found < 0) {
        r2path.push_back(n2path[i]);
        r2bp.push_back(lastp++);
      }
      else {
        r2bp.push_back(found);
      }
    }

    // save point ids
    std::vector<int> n1pids,n2pids;
    for (std::size_t i=0;i<r1path.size();i++) {
      n1pids.push_back(r1path[i].id);
      r1path[i].id = i;
    }
    for (std::size_t i=0;i<r2path.size();i++) {
      n2pids.push_back(r2path[i].id);
      r2path[i].id = i;
    }
    if (ebdebug > 0) {
      std::cout << "# r1p";
      for (std::size_t i=0;i<r1path.size();i++) std::cout << " " << r1path[i].id;
      std::cout << std::endl;
      std::cout << "# r2p";
      for (std::size_t i=0;i<r2path.size();i++) std::cout << " " << r2path[i].id;
      std::cout << std::endl;
      std::cout << "# r1bp";
      for (std::size_t i=0;i<r1bp.size();i++) std::cout << " " << r1bp[i];
      std::cout << std::endl;
      std::cout << "# r2bp";
      for (std::size_t i=0;i<r2bp.size();i++) std::cout << " " << r2bp[i];
      std::cout << std::endl;
    }
    // get triangles
    if (SweepHullDelaunay(r1path,r1path.size(),r1bp,nt1) != 0) 
      throw std::runtime_error("enforce_boundary: sweep_hull_delaunay_trim failed on node path 1");
    if (SweepHullDelaunay(r2path,r2path.size(),r2bp,nt2) != 0) 
      throw std::runtime_error("enforce_boundary: sweep_hull_delaunay_trim failed on node path 2");

    if (ebdebug > 0) {
      std::cout << "# nt1" << std::endl;
      printTriangles(nt1);
      std::cout << "# nt2" << std::endl;
      printTriangles(nt2);
    }

    // restore point ids
    for (std::size_t i=0;i<nt1.size();i++) {
      for (int j=0;j<3;j++) nt1[i].p[j].id = n1pids[nt1[i].p[j].id];
    }
    for (std::size_t i=0;i<nt2.size();i++) {
      for (int j=0;j<3;j++) nt2[i].p[j].id = n2pids[nt2[i].p[j].id];
    }

    // find triangle with shared side
    int st1 = -1;
    int st1s = 0;
    for (std::size_t i=0;i<nt1.size();i++) {
      for (int j=0;j<3;j++) {
	if (pst[0].id == nt1[i].p[j].id && pst[p].id == nt1[i].p[(j+2)%3].id) {
	  st1 = i;
	  st1s = (j+2)%3;
	}
      }
      if (st1 >= 0) break;
    }
    int st2 = -1;
    int st2s = 0;
    for (std::size_t i=0;i<nt2.size();i++) {
      for (int j=0;j<3;j++) {
	if (pst[0].id == nt2[i].p[j].id && pst[p].id == nt2[i].p[(j+1)%3].id) {
	  st2 = i;
	  st2s = j;
	}
      }
      if (st2 >= 0) break;
    }
    if (ebdebug > 0) {
      std::cout << "# st1 " << st1 << " st2 " << st2 << std::endl;
      std::cout << "# st1s " << st1s << " st2s " << st2s << std::endl;
    }

    //
    // put side with most triangles in nt1
    //
    if (nt1.size() < nt2.size()) {
      nt1.swap(nt2);
      std::swap(st1,st2);
      std::swap(st1s,st2s);
    }

    //
    // shift neighbors of nt2
    //
    int numnt1 = nt1.size();
    for (std::size_t i=0;i<nt2.size();i++) {
      for (int j=0;j<3;j++) if (nt2[i].n[j] >= 0) nt2[i].n[j] += numnt1;
    }
    
    //
    // set up neighbors
    //
    nt1[st1].n[st1s] = st2+numnt1;
    nt2[st2].n[st2s] = st1;

    //
    // merge nt2 to nt1
    //
    nt1.insert(nt1.end(),nt2.begin(),nt2.end());
    if (ebdebug > 0) {
      std::cout << "# nt1m" << std::endl;
      printTriangles(nt1);
    }
    
    if (tripath.size() != nt1.size()) throw std::runtime_error("sweep_hull_delaunay enforce_boundary: triangle path and retriangulation have different tri counts");

    //
    // clear all hull information from new triangles
    //
    for (std::size_t i=0;i<nt1.size();i++) {
      for (int j=0;j<3;j++) {
	nt1[i].h[j] = false;
      }
    }
    //
    // merge with mesh
    //
    std::vector<int> trimap(tripath);
    std::sort(trimap.begin(),trimap.end());
    if (ebdebug > 0) {
      std::cout << "# trimap:";
      for (std::size_t i=0;i<trimap.size();i++) std::cout << " " << trimap[i];
      std::cout << std::endl;
    }

    //
    // modify new triangles with updated neighbors
    //
    for (std::size_t i=0;i<nt1.size();i++) {
      for (int j=0;j<3;j++)
	if (nt1[i].n[j] >= 0) nt1[i].n[j] = trimap[nt1[i].n[j]];
    }
    if (ebdebug > 0) {
      std::cout << "# nt1tm" << std::endl;
      printTriangles(nt1);
    }

    //
    // update neighbors/hull info along edge of transplant
    //
    std::vector<int> flipt;
    enforceBoundaryData<Real> d;
    for (std::size_t i=0;i<nt1.size();i++) {
      for (int j=0;j<3;j++) {
	if (nt1[i].n[j] < 0) {
	  // find boundary info
	  d.p = nt1[i].p[j].id;
	  typename std::vector<enforceBoundaryData<Real> >::iterator f = std::lower_bound(bd.begin(),bd.end(),d,enforceBoundaryDataCompare<Real>);
	  if (f == bd.end() || f->p != d.p) throw std::runtime_error("sweep_hull_delaunay enforce_boundary: failed to find saved boundary information");
	  nt1[i].n[j] = f->n;
	  nt1[i].h[j] = f->h;
	  nt1[i].hi[j] = f->hi;
	  if (nt1[i].h[j]) nt1[i].hi[j]->ti = trimap[i];
	  // update neighbor's neighbor list
	  if (nt1[i].n[j] >= 0) {
	    Triangle<Real> & ntn = t[nt1[i].n[j]];
	    for (int k=0;k<3;k++) {
	      if (ntn.p[k].id == d.p) {
		ntn.n[(k+2)%3] = trimap[i];
		break;
	      }
	    }
	  }
	  if (flipt.empty()) flipt.push_back(trimap[i]);
	  else if (flipt.back() != trimap[i]) flipt.push_back(trimap[i]);
	}
      }
    }

    //
    // insert new triangles in the mesh
    //
    for (std::size_t i=0;i<nt1.size();i++) t[trimap[i]] = nt1[i];
    if (ebdebug > 0) {
      std::cout << "# tbf" << std::endl;
      printTriangles(t);
    }

    //
    // fix up triangles along the transplant boundary
    //
    for (std::size_t i=0;i<flipt.size();i++) DelaunayFlip(t,nb,binfo,flipt[i]);
    if (ebdebug > 0) {
      std::cout << "# taf" << std::endl;
      printTriangles(t);
    }

  }

  return pst.size()-2;
}

template < class Real, class RealA >
int BoundDelaunayTriangulation<Real,RealA>
::DelaunayFlip( std::vector<Triangle<Real> > & t,
                const int nb,
                const BoundaryInfo & binfo,
                const int tri )
{
  std::stack<std::pair<int,int> > flipq;

  if (tri < 0 || tri >= int(t.size())) {
    // put last tri on the queue
    flipq.push(std::pair<int,int>(t.size()-1,t.back().n[0]));
    flipq.push(std::pair<int,int>(t.size()-1,t.back().n[1]));
    flipq.push(std::pair<int,int>(t.size()-1,t.back().n[2]));
  }
  else {
    // put desired
    flipq.push(std::pair<int,int>(tri,t[tri].n[0]));
    flipq.push(std::pair<int,int>(tri,t[tri].n[1]));
    flipq.push(std::pair<int,int>(tri,t[tri].n[2]));
  }

  int i=0;
  while (flipq.size() > 0) {
    i++;
    std::pair<int,int> cp = flipq.top();
    flipq.pop();
    
    // check if pair of triangles is not Delaunay

    // skip boundary pair
    if (cp.first == -1 || cp.second == -1) continue;

    // find non-shared point location in second tri
    Point<Real> p;
    int t2n;
    if (t[cp.second].n[0] == cp.first) {
      p = t[cp.second].p[2];
      t2n = 2;
    }
    else if (t[cp.second].n[1] == cp.first) {
      p = t[cp.second].p[0];
      t2n = 0;
    }
    else {
      p = t[cp.second].p[1];
      t2n = 1;
    }
    
    // find non-shared point in first tri
    int t1n;
    if (t[cp.first].n[0] == cp.second) t1n = 2;
    else if (t[cp.first].n[1] == cp.second) t1n = 0;
    else t1n = 1;

    if (nb > 0) {
      // don't flip if shared points align with a boundary
      int id1 = t[cp.first].p[(t1n+1)%3].id;
      int id2 = t[cp.first].p[(t1n+2)%3].id;
      //if (((id1+1)%nb == id2 || (id2+1)%nb == id1) && id1 < nb && id2 < nb) continue;
      if (binfo.onBoundary(id1,id2)) continue;
    }

    // if cc radius is greater than radius to second tri point -> flip
    Real rad2 = (p.x-t[cp.first].cc.x)*(p.x-t[cp.first].cc.x)
      + (p.y-t[cp.first].cc.y)*(p.y-t[cp.first].cc.y);

    // also flip if non-shared points align with a boundary and make a proper triangle pair
    bool dflip = false;
    if (nb > 0) {
      int id1 = t[cp.first].p[t1n].id;
      int id2 = t[cp.second].p[t2n].id;
      //if (((id1+1)%nb == id2 || (id2+1)%nb == id1) && id1 < nb && id2 < nb) {
      if (binfo.onBoundary(id1,id2)) {
	Real x1 = t[cp.first].p[(t1n+1)%3].x-t[cp.first].p[t1n].x;
	Real x2 = t[cp.second].p[t2n].x-t[cp.first].p[t1n].x;
	Real y1 = t[cp.first].p[(t1n+1)%3].y-t[cp.first].p[t1n].y;
	Real y2 = t[cp.second].p[t2n].y-t[cp.first].p[t1n].y;
	Real det1 = x1*y2-x2*y1;
	x1 = t[cp.second].p[t2n].x-t[cp.first].p[t1n].x;
	x2 = t[cp.first].p[(t1n+2)%3].x-t[cp.first].p[t1n].x;
	y1 = t[cp.second].p[t2n].y-t[cp.first].p[t1n].y;
	y2 = t[cp.first].p[(t1n+2)%3].y-t[cp.first].p[t1n].y;
	Real det2 = x1*y2-x2*y1;
	if (det1 > 0.0 && det2 > 0.0) dflip = true;
      }
    }
    if ((t[cp.first].cc.rc2 > rad2 && fabs(t[cp.first].cc.rc2-rad2) > 1.e-12*fabs(t[cp.first].cc.rc2+rad2)) || dflip) {
      // flip triangles
      Triangle<Real> newt1,newt2;

      newt1.p[0] = t[cp.first].p[t1n];
      newt1.p[1] = t[cp.first].p[(t1n+1)%3];
      newt1.p[2] = t[cp.second].p[t2n];
      newt1.n[0] = t[cp.first].n[t1n];
      newt1.n[1] = t[cp.second].n[(t2n+2)%3];
      newt1.n[2] = cp.second;
      newt1.h[0] = t[cp.first].h[t1n];
      newt1.h[1] = t[cp.second].h[(t2n+2)%3];
      newt1.h[2] = false;
      newt1.hi[0] = t[cp.first].hi[t1n];
      newt1.hi[1] = t[cp.second].hi[(t2n+2)%3];
      if (newt1.n[1] >= 0) updateTriangleNeighbor(t[newt1.n[1]],newt1.p[2].id,cp.first);
      if (newt1.h[1]) newt1.hi[1]->ti = cp.first;
      newt2.p[0] = t[cp.second].p[t2n];
      newt2.p[1] = t[cp.second].p[(t2n+1)%3];
      newt2.p[2] = t[cp.first].p[t1n];
      newt2.n[0] = t[cp.second].n[t2n];
      newt2.n[1] = t[cp.first].n[(t1n+2)%3];
      newt2.n[2] = cp.first;
      newt2.h[0] = t[cp.second].h[t2n];
      newt2.h[1] = t[cp.first].h[(t1n+2)%3];
      newt2.h[2] = false;
      newt2.hi[0] = t[cp.second].hi[t2n];
      newt2.hi[1] = t[cp.first].hi[(t1n+2)%3];
      if (newt2.n[1] >= 0) updateTriangleNeighbor(t[newt2.n[1]],newt2.p[2].id,cp.second);
      if (newt2.h[1]) newt2.hi[1]->ti = cp.second;
      if (findCircumcircle(newt1) != 0 || findCircumcircle(newt2) != 0) {
	std::cout << "# t1n " << t1n << " t2n " << t2n << std::endl;
	std::cout << "# nt1 n1 " << newt1.n[0] << " n2 " << newt1.n[1] << " n3 " << newt1.n[2] << std::endl
		  << "x " << newt1.p[0].x << " y " << newt1.p[0].y << std::endl
		  << "x " << newt1.p[1].x << " y " << newt1.p[1].y << std::endl
		  << "x " << newt1.p[2].x << " y " << newt1.p[2].y << std::endl
		  << "x " << newt1.p[0].x << " y " << newt1.p[0].y << std::endl
		  << std::endl << std::endl;
	std::cout << "# nt2 n1 " << newt2.n[0] << " n2 " << newt2.n[1] << " n3 " << newt2.n[2] << std::endl
		  << "x " << newt2.p[0].x << " y " << newt2.p[0].y << std::endl
		  << "x " << newt2.p[1].x << " y " << newt2.p[1].y << std::endl
		  << "x " << newt2.p[2].x << " y " << newt2.p[2].y << std::endl
		  << "x " << newt2.p[0].x << " y " << newt2.p[0].y << std::endl
		  << std::endl << std::endl;
      	std::cout << "error in delaunay_flip: failed to find required tri circumcircle" << std::endl;
	return -1;
      }
      t[cp.first] = newt1;
      t[cp.second] = newt2;

      // add additional pairs to test for further flipping
      flipq.push(std::pair<int,int>(cp.first,t[cp.first].n[0]));
      flipq.push(std::pair<int,int>(cp.first,t[cp.first].n[1]));
      flipq.push(std::pair<int,int>(cp.second,t[cp.second].n[0]));
      flipq.push(std::pair<int,int>(cp.second,t[cp.second].n[1]));

    }
  }

  return 0;
}

#endif
