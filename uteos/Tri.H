/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SUTTG, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#ifndef TRI_H
#define TRI_H

#include <vector>
#include <cmath>

#include "EOSData.H"
#include "Node.H"

class Tri {
public:
  Tri( int n1,
       int n2,
       int n3 )
    : maxerr(-1.0) { 
    n[0] = n1;
    n[1] = n2;
    n[2] = n3;
    maxerrloc[0] = -1.;
    maxerrloc[1] = -1.;
  }

  ~Tri() { }

  // find barycentric coords
  void getBcoords( const double x,
                   const double y,
                   std::vector<double> & b )
  {
    b[0] = a[1]*(a[5]-y)-(a[2]-x)*a[4];
    b[1] = a[3]*(a[2]-x)-(a[5]-y)*a[0];
    b[2] = 1.-b[0]-b[1];
    if (fabs(b[0]) < 2.e-14) b[0] = 0.;
    if (fabs(b[1]) < 2.e-14) b[1] = 0.;
    if (fabs(b[2]) < 2.e-14) b[2] = 0.;
  }
  
  void computeBcoords( const std::vector<Node> & nodes,
                       const int logvars )
  {
    double px[3],py[3];

    // compute barycentric coordinate helpers
    if (logvars) {
      for (int j=0;j<3;j++) {
        px[j] = log(nodes[n[j]].d.inputs[0]);
        py[j] = log(nodes[n[j]].d.inputs[1]);
      }
    }
    else {
      for (int j=0;j<3;j++) {
        px[j] = nodes[n[j]].d.inputs[0];
        py[j] = nodes[n[j]].d.inputs[1];
      }
    }
    a[0] = px[0]-px[2];
    a[1] = px[1]-px[2];
    a[2] = px[2];
    a[3] = py[0]-py[2];
    a[4] = py[1]-py[2];
    a[5] = py[2];
    double invdet = 1./(a[0]*a[4]-a[1]*a[3]);
    a[0] *= invdet;
    a[1] *= invdet;
    a[3] *= invdet;
    a[4] *= invdet;
  }

  // node indices in anti-clockwise direction
  int n[3];

  // maximum error on triangle
  double maxerr;

  // location of maximum error
  double maxerrloc[2];

  // vector for faster calculation of barycentric coordinates
  double a[6];
};

inline bool operator < ( const Tri & lhs,
                         const Tri & rhs )
{
  if (lhs.n[0] < rhs.n[0]) return true;
  else if (lhs.n[0] == rhs.n[0]) {
    if (lhs.n[1] < rhs.n[1]) return true;
    else if (lhs.n[1] == rhs.n[1]) {
      if (lhs.n[2] < rhs.n[2]) return true;
      else return false;
    }
    else return false;
  }
  else return false;
}

inline bool operator == ( const Tri & lhs,
                          const Tri & rhs )
{
  if (lhs.n[0] == rhs.n[0] && lhs.n[1] == rhs.n[1] && lhs.n[2] == rhs.n[2]) return true;
  else return false;
}

inline bool operator != ( const Tri & lhs, const Tri & rhs ) { return !operator==(lhs,rhs); }
inline bool operator >  ( const Tri & lhs, const Tri & rhs ) { return  operator< (rhs,lhs); }
inline bool operator <= ( const Tri & lhs, const Tri & rhs ) { return !operator> (lhs,rhs); }
inline bool operator >= ( const Tri & lhs, const Tri & rhs ) { return !operator< (lhs,rhs); }

#endif
