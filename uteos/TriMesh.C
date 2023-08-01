/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SUTTG, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#include <stdexcept>
#include <iostream>

#include "TriMesh.H"

void TriMesh::evaluate_boundary( int bnum,
                                 int side,
                                 MeshTypes meshtype,
                                 EOSData & d,
                                 double x )
{
  Boundary & b = boundaries[bnum];
  int bphase;
  double r1,r2;
  double t1,t2;
  double targetT;

  switch (b.type) {
  case 0:
    //
    // evaluate via phase boundary info
    //
    int p1,pm,p2,pl,ph;
    pbis[0][b.pbmnum]->getPhases(p1,pm,p2,pl,ph);
    if (side == 1) {
      if (b.pbmside == 0) bphase = p1;
      else bphase = p2;
    }
    else bphase = pm;

    //
    // phase boundaries currently lookup in temperature so find value given x
    //
    // grab the temperature or pressure from the endnodes to get the target
    //
    if (meshtype == RT_MESH) {
      t1 = nodes[b.endnodes1[0]].d.inputs[1];
      t2 = nodes[b.endnodes2[0]].d.inputs[1];
    }
    else if (meshtype == RE_MESH) {
      t1 = nodes[b.endnodes1[0]].d.outputs[1];
      t2 = nodes[b.endnodes2[0]].d.outputs[1];
    }
    else if (meshtype == PE_MESH || meshtype == PS_MESH || meshtype == PH_MESH) {
      t1 = nodes[b.endnodes1[0]].d.inputs[0];
      t2 = nodes[b.endnodes2[0]].d.inputs[0];
    }
    else {
      throw std::runtime_error("TriMesh::evaluate_boundary: invalid mesh type");
    }
    targetT = t1+(t2-t1)*x;

    if (meshtype == PE_MESH || meshtype == PS_MESH || meshtype == PH_MESH) {
      //
      // Assume that PX mesh models are using exact boundary calculations and always return enthalpy
      // targetT is actually pressure here
      //
      double val1,val2;
      pbis[0][b.pbmnum]->getExactVars(targetT,val1,val2);
      d.inputs[0] = targetT;
      if (b.pbmside == 0) d.inputs[1] = val1;
      else d.inputs[1] = val2;
      model->getThermalEOS_PH_D(d,bphase);
      //
      // swap to appropriate PX space
      //
      if (meshtype == PE_MESH) std::swap(d.inputs[1],d.outputs[2]);
      else if (meshtype == PS_MESH) std::swap(d.inputs[1],d.outputs[3]);
    }
    else throw std::runtime_error("TriMesh::evaluate_boundary: invalid mesh type for pbm boundary");

    break;
  case 1:
    //
    // constant T boundary -> set rho or P from x, then set T and evaluate
    //
    if (side == 1) {
      bphase = b.phase;
      // check special flag for triple point
      if (b.tpflag) bphase = -bphase;
      r1 = nodes[b.endnodes1[0]].d.inputs[0];
      r2 = nodes[b.endnodes2[0]].d.inputs[0];
    }
    else {
      bphase = b.phase2;
      // check special flag for triple point
      if (b.tpflag2) bphase = -bphase;
      r1 = nodes[b.endnodes1[1]].d.inputs[0];
      r2 = nodes[b.endnodes2[1]].d.inputs[0];
    }
    
    d.inputs[0] = r1+(r2-r1)*x;
    d.inputs[1] = b.T;
    if (meshtype == RT_MESH || meshtype == RE_MESH) model->getThermalEOS_F_D(d,bphase);
    //
    // note that the PT method will throw if the constant T boundary
    // crosses a phase boundary, as this would cause a discontinuity
    // in the boundary (multivalued function values)
    //
    else if (meshtype == PE_MESH || meshtype == PH_MESH || meshtype == PS_MESH) model->getThermalEOS_PT_D(d,bphase);
    else throw std::runtime_error("TriMesh::evaluate_boundary: invalid mesh type for constant T boundary");

    //
    // for RE mesh we need to swap T and E
    //
    if (meshtype == RE_MESH) {
      std::swap(d.inputs[1],d.outputs[1]);
    }

    //
    // for PX meshes we need to swap T and X
    //
    if (meshtype == PE_MESH) {
      std::swap(d.inputs[1],d.outputs[1]);
      std::swap(d.inputs[1],d.outputs[2]);
    }
    if (meshtype == PS_MESH) {
      std::swap(d.inputs[1],d.outputs[1]);
      std::swap(d.inputs[1],d.outputs[3]);
    }
    if (meshtype == PH_MESH) {
      std::swap(d.inputs[1],d.outputs[1]);
    }

    break;
  case 2:
    //
    // constant rho boundary -> set T or E from x, then set rho and evaluate
    //
    if (side == 1) {
      bphase = b.phase;
      t1 = nodes[b.endnodes1[0]].d.inputs[1];
      t2 = nodes[b.endnodes2[0]].d.inputs[1];
    }
    else {
      bphase = b.phase2;
      t1 = nodes[b.endnodes1[1]].d.inputs[1];
      t2 = nodes[b.endnodes2[1]].d.inputs[1];
    }

    d.inputs[0] = b.rho;
    d.inputs[1] = t1+(t2-t1)*x;
    if (meshtype == RT_MESH) model->getThermalEOS_F_D(d,bphase);
    else if (meshtype == RE_MESH) model->getThermalEOS_S_D(d,bphase);
    else throw std::runtime_error("TriMesh::evaluate_boundary: invalid mesh type for constant rho boundary");

    break;
  case 3:
    //
    // constant P boundary -> set E, S, or H from x, then set P and evaluate
    //
    if (side == 1) {
      bphase = b.phase;
      t1 = nodes[b.endnodes1[0]].d.inputs[1];
      t2 = nodes[b.endnodes2[0]].d.inputs[1];
    }
    else {
      bphase = b.phase2;
      t1 = nodes[b.endnodes1[1]].d.inputs[1];
      t2 = nodes[b.endnodes2[1]].d.inputs[1];
    }

    d.inputs[0] = b.P;
    d.inputs[1] = t1+(t2-t1)*x;
    if (meshtype == PE_MESH) model->getThermalEOS_PE_D(d,bphase);
    else if (meshtype == PS_MESH) model->getThermalEOS_PS_D(d,bphase);
    else if (meshtype == PH_MESH) model->getThermalEOS_PH_D(d,bphase);
    else throw std::runtime_error("TriMesh::evaluate_boundary: invalid mesh type for constant P boundary");

    break;
  default:
    throw std::runtime_error("evaluate_boundary: invalid type in boundary");
    break;
  }
}

void TriMesh::evaluate_region( int phase,
                               MeshTypes meshtype,
                               EOSData & d )
{
  if (meshtype == RT_MESH) model->getThermalEOS_F_D(d,phase);
  else if (meshtype == RE_MESH) model->getThermalEOS_S_D(d,phase);
  else if (meshtype == PH_MESH) model->getThermalEOS_PH_D(d,phase);
  else if (meshtype == PS_MESH) model->getThermalEOS_PS_D(d,phase);
  else model->getThermalEOS_PE_D(d,phase); // PE_MESH
}

void TriMesh::nodes_swap_RTtoRE()
{
  for (std::size_t i=0;i<nodes.size();i++) {
    std::swap(nodes[i].d.inputs[1],nodes[i].d.outputs[1]);
  }
}

void TriMesh::nodes_swap_PTtoPX( const MeshTypes mt )
{
  if (mt == PE_MESH) {
    for (std::size_t i=0;i<nodes.size();i++) {
      std::swap(nodes[i].d.inputs[1],nodes[i].d.outputs[1]);
      std::swap(nodes[i].d.inputs[1],nodes[i].d.outputs[2]);
    }
  }
  else if (mt == PS_MESH) {
    for (std::size_t i=0;i<nodes.size();i++) {
      std::swap(nodes[i].d.inputs[1],nodes[i].d.outputs[1]);
      std::swap(nodes[i].d.inputs[1],nodes[i].d.outputs[3]);
    }
  }
  else if (mt == PH_MESH) {
    for (std::size_t i=0;i<nodes.size();i++) {
      std::swap(nodes[i].d.inputs[1],nodes[i].d.outputs[1]);
    }
  }
  else throw std::runtime_error("TriMesh::nodes_swap_PTtoPX: invalid meshtype");
}

void TriMesh::nodes_swap_PHtoPX( const MeshTypes mt )
{
  if (mt == PE_MESH) {
    for (std::size_t i=0;i<nodes.size();i++) {
      std::swap(nodes[i].d.inputs[1],nodes[i].d.outputs[2]);
    }
  }
  else if (mt == PS_MESH) {
    for (std::size_t i=0;i<nodes.size();i++) {
      std::swap(nodes[i].d.inputs[1],nodes[i].d.outputs[3]);
    }
  }
  else {
    if (mt != PH_MESH) {
      throw std::runtime_error("TriMesh::nodes_swap_PHtoPX: invalid meshtype");
    }
  }
}
