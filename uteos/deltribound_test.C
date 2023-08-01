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
#include <vector>
#include <fstream>
#include <cstdlib>

#include "BoundDelaunayTriangulation.H"
#include "qd/qd_real.h"

void read_file( const char * filename,
                std::vector<double> & x,
                std::vector<double> & y,
                double & xmin,
                double & xmax,
                double & ymin,
                double & ymax,
                int & nb )
{
  std::ifstream in(filename, std::ios::in);

  if (in.good()) {
    std::string junk;
    in >> junk;
    in >> nb;
    in >> junk;
    in >> xmin;
    in >> xmax;
    in >> junk;
    in >> ymin;
    in >> ymax;
  }
  else throw std::runtime_error("read_file: failed to open file");

  while (in.good()) {
    std::string junk;
    double p;
    in >> junk;
    if (!in.good()) break;
    in >> p;
    x.push_back(p);
    in >> junk;
    in >> p;
    y.push_back(p);
  }
}

template<>
double convertRealAtoReal<double,qd_real>( const qd_real & var )
{
  return to_double(var);
}


int main(int argc,char *argv[])
{
  if (argc < 2) throw std::runtime_error("need input file name");

  std::vector<double> px,py;
  double xmin,xmax,ymin,ymax;
  int nb;

  read_file(argv[1],px,py,xmin,xmax,ymin,ymax,nb);

  if (px.size() != py.size()) throw std::runtime_error("px != py");

  unsigned int old_cw;
  fpu_fix_start(&old_cw);

  BoundDelaunayTriangulation<double,qd_real> bdt(px,py,nb,xmax-xmin,ymax-ymin);
  bdt.triangulate();
  fpu_fix_end(&old_cw);
  const std::vector<int> & t = bdt.getTriangles();
  for (std::size_t i=0;i<px.size();i++) {
    std::cout << "n " << i << " x " << px[i] << " y " << py[i] << std::endl;
  }
  std::cout << std::endl << std::endl;

  std::cout << "# nt " << t.size()/3 << std::endl;
  for (std::size_t i=0;i<t.size()/3;i++) {
    std::cout << "# t " << i << " n " << t[3*i+0] << " " << t[3*i+1] << " " << t[3*i+2] << std::endl;
    std::cout << "x1 " << px[t[3*i+0]] << " y1 " << py[t[3*i+0]] << std::endl;
    std::cout << "x2 " << px[t[3*i+1]] << " y2 " << py[t[3*i+1]] << std::endl;
    std::cout << "x3 " << px[t[3*i+2]] << " y3 " << py[t[3*i+2]] << std::endl;
    std::cout << "x1 " << px[t[3*i+0]] << " y1 " << py[t[3*i+0]] << std::endl;
    std::cout << std::endl;
  }


  return 0;
}
