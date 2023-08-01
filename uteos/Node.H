/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SUTTG, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#ifndef NODE_H
#define NODE_H

#include <vector>

#include "EOSData.H"

class Node {
public:
  Node( int ni,
        int no ) :
    d(ni,no) { }

  ~Node() { }

  // node model data
  EOSData d;
  // node location in boundary
  double x;
  // list of nodes which share locations
  std::vector<int> shared;

  // index to triangles this node belongs to
  std::vector<int> tris;

  // min and max errors of triangles associated with this node
  double minerr;
  double maxerr;

  bool operator == ( const Node & rhs ) const {
    if (d.inputs == rhs.d.inputs && d.outputs == rhs.d.outputs) return true;
    else return false;
  }

};

#endif
