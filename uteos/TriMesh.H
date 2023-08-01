/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SUTTG, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#ifndef TRI_MESH_H
#define TRI_MESH_H

#include <vector>

#include "PhaseBoundaryInfo.H"
#include "EOSModel.H"
#include "Node.H"
#include "Tri.H"

class Boundary {
public:
  Boundary() {}

  Boundary( const int t,
	    const std::vector<int> & en1,
	    const std::vector<int> & en2,
	    const double x,
	    const int p,
	    const int p2 = 0,
            const int tpf = 0,
            const int tpf2 = 0 ) :
    type(t),pbmnum(-1),pbmside(-1),pbmnum2(-1),pbmside2(-1),
    phase(p),phase2(p2),tpflag(tpf),tpflag2(tpf2),T(0.),rho(0.),P(0.),
    endnodes1(en1),endnodes2(en2)
  {
    if (type == 1) T = x;
    else if (type == 2) rho = x;
    else if (type == 3) P = x;
    else throw std::logic_error("Boundary: invalid constant type");
  }

  Boundary( const int t,
	    const std::vector<int> & en1,
	    const std::vector<int> & en2,
	    const int pn,
	    const int ps,
	    const int pn2 = -1,
	    const int ps2 = -1 ) :
    type(t),pbmnum(pn),pbmside(ps),pbmnum2(pn2),pbmside2(ps2),
    phase(0),phase2(0),tpflag(0),tpflag2(0),T(0.),rho(0.),P(0.),
    endnodes1(en1),endnodes2(en2)
  {
    if (type != 0) throw std::logic_error("Boundary: invalid pbm type");
  }

  // kind of boundary:
  //  0 -- pbm boundary defined by pbm num and pbm side
  //  1 -- constant T boundary defined by phase and T
  //  2 -- constant rho boundary defined by phase and rho
  //  3 -- constant P boundary defined by phase and P
  int type;
  int pbmnum;   // boundary index in pbm
  int pbmside;  // phase index in pbm
  int pbmnum2;  // boundary index in pbm
  int pbmside2; // phase index in pbm
  int phase;    // phase number
  int phase2;   // second phase number
  int tpflag;   // flag for which side of the const T triple point the boundary lies
  int tpflag2;  //  -- nonzero means on the low energy side
  double T;     // fixed boundary temperature
  double rho;   // fixed boundary density
  double P;     // fixed boundary pressure
  // node ids at ends of boundary
  std::vector<int> endnodes1;
  std::vector<int> endnodes2;
};

class Region {
public:
  //
  // index to region boundaries in anti-clockwise direction
  //
  std::vector<int> bounds;

  //
  // phase the region represents
  //
  int phase;
};

class TriMesh {
public:
  enum MeshTypes {
    RT_MESH,
    RE_MESH,
    PH_MESH,
    PS_MESH,
    PE_MESH,
    TOT_MESH_TYPES
  };

  TriMesh( std::vector<PhaseBoundaryInfo *> * phaseboundaries,
           EOSModel * eosmodel ) :
    pbis(phaseboundaries), model(eosmodel) { }

  ~TriMesh() {}

  void evaluate_boundary( int bnum,
                          int side,
                          MeshTypes meshtype,
                          EOSData & d,
                          double x );

  void evaluate_region( int phase,
                        MeshTypes meshtype,
                        EOSData & d );

  void nodes_swap_RTtoRE( );
  void nodes_swap_PTtoPX( const MeshTypes mt );
  void nodes_swap_PHtoPX( const MeshTypes mt );

  std::vector<PhaseBoundaryInfo *> * pbis;
  std::vector<Node> nodes;
  std::vector<Tri> tris;
  std::vector<Boundary> boundaries;
  std::vector<Region> regions;
  EOSModel * model;
  
};

#endif
