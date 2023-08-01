/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SUTTG, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#ifndef BOUNDARY_INFO_H
#define BOUNDARY_INFO_H

//
// BoundaryInfo.H:
//
// A class for holding information about points that form a
// boundary. Points are numbered along the boundary. The boundary is
// assumed to be defined in an anti-clockwise manner.
//

#include <map>       // for std::multimap
#include <vector>    // for std::vector
#include <algorithm> // for std::sort

class BoundaryInfo
{
public:
  //
  // A triangle with nodes A,B,C oriented anti-clockwise
  //
  struct btri {
    btri( const int A,
          const int B,
          const int C )
      : a(A),b(B),c(C) {}

    int a,b,c;
  };

  //
  // Constructor
  //
  BoundaryInfo( const int num_points,
                const std::vector<int> & boundary_point_list )
    : nbp(num_points), bl(boundary_point_list), bpm()
  {
    //
    // make boundary point map
    //
    for (std::size_t i=0;i<bl.size();i++) bpm.insert(std::pair<int,int>(bl[i],i));

  }

  //
  // Destructor
  //
  ~BoundaryInfo() {}

  //
  // return true if the points are next to each other on the boundary
  //
  bool onBoundary( const int p1,
                   const int p2 ) const
  {
    // first check that both points are individually on the boundary
    if (p1 >= nbp || p2 >= nbp) return false;

    // get locations in boundary for p1
    std::pair<std::multimap<int,int>::const_iterator, std::multimap<int,int>::const_iterator> bi
      = bpm.equal_range(p1);

    for (std::multimap<int,int>::const_iterator it=bi.first; it!=bi.second; ++it) {
      if (bl[(it->second+1)%nbp] == p2 || bl[(it->second-1+nbp)%nbp] == p2)
        return true;
    }

    return false;
  }

  //
  // return the neighbors of a certain boundary point
  //
  void getNeighbors( const int p,
                     std::vector<int> & n ) const
  {
    // get locations in boundary for p
    std::pair<std::multimap<int,int>::const_iterator, std::multimap<int,int>::const_iterator> bi
      = bpm.equal_range(p);

    // get all neighbors
    std::vector<int> neighbors;
    for (std::multimap<int,int>::const_iterator it=bi.first; it!=bi.second; ++it) {
      neighbors.push_back(bl[(it->second+1)%nbp]);
      neighbors.push_back(bl[(it->second-1+nbp)%nbp]);
    }

    // only return unique neighbors
    std::sort(neighbors.begin(),neighbors.end());
    n.clear();
    n.push_back(neighbors[0]);
    for (std::size_t i=1;i<neighbors.size();i++) {
      if (neighbors[i] != n.back()) n.push_back(neighbors[i]);
    }
  }

  //
  // return true if the point indices define a triangle in
  // anti-clockwise orientation that lies inside boundary
  //
  bool isInside( const int a,
                 const int b,
                 const int c ) const
  {
    // if any are not on the boundary we are inside
    if (a >= nbp || b >= nbp || c >= nbp) return true;

    // generate all possible boundary index combinations
    std::pair<std::multimap<int,int>::const_iterator, std::multimap<int,int>::const_iterator> bi;
    bi = bpm.equal_range(a);
    std::vector<btri> bic1;
    for (std::multimap<int,int>::const_iterator it=bi.first; it!=bi.second; ++it) {
      bic1.push_back(btri(it->second,-1,-1));
    }
    bi = bpm.equal_range(b);
    std::vector<btri> bic2;
    for (std::size_t i=0;i<bic1.size();i++) {
      for (std::multimap<int,int>::const_iterator it=bi.first; it!=bi.second; ++it) {
        bic2.push_back(btri(bic1[i].a,it->second,-1));
      }
    }
    bi = bpm.equal_range(c);
    std::vector<btri> bic3;
    for (std::size_t i=0;i<bic2.size();i++) {
      for (std::multimap<int,int>::const_iterator it=bi.first; it!=bi.second; ++it) {
        bic3.push_back(btri(bic2[i].a,bic2[i].b,it->second));
      }
    }

    for (std::size_t i=0;i<bic3.size();i++) {
      // inside if any combination is true
      if (bic3[i].b > bic3[i].a) {
        if (bic3[i].c > bic3[i].b || bic3[i].a > bic3[i].c) return true;
      }
      else {
        if (bic3[i].c > bic3[i].b && bic3[i].c < bic3[i].a) return true;
      }
    }

    return false;

  }

private:

  // number of boundary points
  const int nbp;

  // reference to anti-clockwise list of points along boundary
  const std::vector<int> & bl;

  // map of boundary point number to list index
  std::multimap<int,int> bpm;

};

#endif
