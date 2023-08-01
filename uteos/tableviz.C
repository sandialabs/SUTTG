/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SUTTG, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

//
// tableviz.C:
//
// This is a simple program to take a utri table file and output a
// sample grid along with the triangles stored in the file for later
// visualization.

#include <stdexcept>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>

#include "utri_eos_support.h"
#include "utri_eos_mig.h"

void print_triangles(eos_table_utri_t & etu,int type)
{
  std::cout << "# index 1-" << etu.nt << " : triangle output" << std::endl;
  std::cout << "# x y r t p e s f dpdr dpdT dEdr dEdT cs dcsdr phase visc thcon" << std::endl;
  for (int i=0;i<etu.nt;i++) {
    eos_table_tri_t & t = etu.tris[i];
    std::cout << "# t " << i << " n1 " << t.p[0] << " n2 " << t.p[1] << " n3 " << t.p[2] << std::endl;
    for (int j=0;j<4;j++) {
      if (type == 1) {
	std::cout << " " << etu.modes[0]->rho[t.p[j%3]] 
		  << " " << etu.modes[0]->T[t.p[j%3]];
      }
      else if (type == 2) {
	std::cout << " " << etu.modes[0]->rho[t.p[j%3]] 
		  << " " << etu.modes[0]->E[t.p[j%3]];
      }
      else if (type == 3) {
	std::cout << " " << etu.modes[0]->p[t.p[j%3]] 
		  << " " << etu.modes[0]->H[t.p[j%3]];
      }
      else if (type == 4) {
	std::cout << " " << etu.modes[0]->p[t.p[j%3]] 
		  << " " << etu.modes[0]->S[t.p[j%3]];
      }
      else {
	std::cout << " " << etu.modes[0]->p[t.p[j%3]] 
		  << " " << etu.modes[0]->E[t.p[j%3]];
      }
      std::cout << " " << etu.modes[0]->rho[t.p[j%3]] 
		<< " " << etu.modes[0]->T[t.p[j%3]]
		<< " " << etu.modes[0]->p[t.p[j%3]]
		<< " " << etu.modes[0]->E[t.p[j%3]]
		<< " " << etu.modes[0]->S[t.p[j%3]]
		<< " " << etu.modes[0]->F[t.p[j%3]]
		<< " " << etu.modes[0]->dpdr[t.p[j%3]]
		<< " " << etu.modes[0]->dpdT[t.p[j%3]]
		<< " " << etu.modes[0]->dEdr[t.p[j%3]]
		<< " " << etu.modes[0]->dEdT[t.p[j%3]]
		<< " " << etu.modes[0]->cs[t.p[j%3]]
		<< " " << etu.modes[0]->dcsdr[t.p[j%3]]
		<< " " << etu.modes[0]->phase[t.p[j%3]]
		<< " " << etu.modes[0]->visc[t.p[j%3]]
		<< " " << etu.modes[0]->thcon[t.p[j%3]]
		<< std::endl;
    }
    std::cout << std::endl << std::endl;
  }

}

int main(int argc,char * argv[])
{
  //
  // get the table from the command line
  //
  if (argc < 3) {
    std::cerr << "Usage: utriTableDump <tablefile> <tabletype> <mode0> <mode1> ... <modeN>" << std::endl;
    return 1;
  }

  //
  // read in the table file
  //
  int utriHandle = utri_eos_handle_open();
  eos_table_utri_t etu;
  if (eostableutrireadnc(argv[1],&etu,utriHandle) != 0) {
    throw std::runtime_error("Failed to open input table file");
  }

  int type = atoi(argv[2]);

  //
  // mix in other modes as specified on command line
  //
  std::vector<double> modes(etu.nmodes,0.);
  modes[0] = 1.;
  for (int i=0;i<argc-3 && i<etu.nmodes;i++) {
    modes[i] = atof(argv[i+3]);
  }

  std::cout << "# scaling mean mode by " << modes[0] << std::endl;
  for (int j=0;j<etu.np;j++) {
    // don't scale the coordinates so we can still output on a good grid
    if (type == 2) etu.modes[0]->T[j]     *= modes[0];
    etu.modes[0]->p[j]     *= modes[0];
    if (type == 1) etu.modes[0]->E[j]     *= modes[0];
    etu.modes[0]->S[j]     *= modes[0];
    etu.modes[0]->F[j]     *= modes[0];
    etu.modes[0]->dpdr[j]  *= modes[0];
    etu.modes[0]->dpdT[j]  *= modes[0];
    etu.modes[0]->dEdr[j]  *= modes[0];
    etu.modes[0]->dEdT[j]  *= modes[0];
    etu.modes[0]->cs[j]    *= modes[0];
    etu.modes[0]->dcsdr[j] *= modes[0];
    etu.modes[0]->visc[j]    *= modes[0];
    etu.modes[0]->thcon[j]    *= modes[0];
  }

  std::cout << "# mixing in " << etu.nmodes-1 << " modes with strengths: ";
  for (std::size_t i=1;i<modes.size();i++) std::cout << " m" << i << " = " << modes[i];
  std::cout << std::endl << std::endl;

  for (std::size_t i=1;i<modes.size();i++) {
    for (int j=0;j<etu.np;j++) {
      etu.modes[0]->rho[j]   += modes[i]*etu.modes[i]->rho[j];
      etu.modes[0]->T[j]     += modes[i]*etu.modes[i]->T[j];
      etu.modes[0]->p[j]     += modes[i]*etu.modes[i]->p[j];
      etu.modes[0]->E[j]     += modes[i]*etu.modes[i]->E[j];
      etu.modes[0]->S[j]     += modes[i]*etu.modes[i]->S[j];
      etu.modes[0]->F[j]     += modes[i]*etu.modes[i]->F[j];
      etu.modes[0]->dpdr[j]  += modes[i]*etu.modes[i]->dpdr[j];
      etu.modes[0]->dpdT[j]  += modes[i]*etu.modes[i]->dpdT[j];
      etu.modes[0]->dEdr[j]  += modes[i]*etu.modes[i]->dEdr[j];
      etu.modes[0]->dEdT[j]  += modes[i]*etu.modes[i]->dEdT[j];
      etu.modes[0]->cs[j]    += modes[i]*etu.modes[i]->cs[j];
      etu.modes[0]->dcsdr[j] += modes[i]*etu.modes[i]->dcsdr[j];
      etu.modes[0]->visc[j]  += modes[i]*etu.modes[i]->visc[j];
      etu.modes[0]->thcon[j] += modes[i]*etu.modes[i]->thcon[j];
    }
  }
    
  // convert independent variables to log values if needed
  // this must be done before rptree regeneration
  //
  if (etu.flags[0] == 1) {
    if (type == 1) {
      for (int i=0;i<etu.np;i++) etu.modes[0]->rho[i] = log(etu.modes[0]->rho[i]);
      for (int i=0;i<etu.np;i++) etu.modes[0]->T[i] = log(etu.modes[0]->T[i]);
    }
    else if (type == 2) {
      for (int i=0;i<etu.np;i++) etu.modes[0]->rho[i] = log(etu.modes[0]->rho[i]);
      for (int i=0;i<etu.np;i++) etu.modes[0]->E[i] = log(etu.modes[0]->E[i]);
    }      
    else {
      for (int i=0;i<etu.np;i++) etu.modes[0]->p[i] = log(etu.modes[0]->p[i]);
      if (type == 3) for (int i=0;i<etu.np;i++) etu.modes[0]->H[i] = log(etu.modes[0]->H[i]);
      else if (type == 4) for (int i=0;i<etu.np;i++) etu.modes[0]->S[i] = log(etu.modes[0]->S[i]);
      else if (type == 5) for (int i=0;i<etu.np;i++) etu.modes[0]->E[i] = log(etu.modes[0]->E[i]);
    }
  }

  //
  // Fix the boundary
  //
  if (type == 2 || type == 5) {
    if (construct_boundary_tris(&etu,type) != 0) {
      throw std::runtime_error("Failed to construct boundaries");
    }
  }

  //
  // regenerate the rptree
  //
  if (calculate_rptree(&etu,type-1) != 0) {
    throw std::runtime_error("Failed to regenerate rptree");
  }

  //
  // allocate space for calculations
  //
  eos_data_points_t *f = NULL;
  if ((f = (eos_data_points_t *)eospointsrealloc(f,1)) == NULL) {
    throw std::runtime_error("Failed to allocate data points");
  }

  //
  // shut up uninitialized value warning in printouts
  //
  f->rho[0] = 0.;  f->T[0] = 0.;  f->p[0] = 0.;  f->E[0] = 0.;
  f->S[0] = 0.;  f->F[0] = 0.;
  f->dpdr[0] = 0.;  f->dpdT[0] = 0.;  f->dEdr[0] = 0.;  f->dEdT[0] = 0.;
  f->cs[0] = 0.;  f->dcsdr[0] = 0.;  f->visc[0] = 0.;  f->thcon[0] = 0.;
  f->phase[0] = 0;

  //
  // evaluate on a grid
  //
  int numxgrid = 100;
  int numygrid = 100;
  int clip[4];
  double scratch[5];
  int iscratch[1];
  std::cout << "# index 0 : grid evaluation" << std::endl;
  std::cout << "# x y r t p e s f dpdr dpdT dEdr dEdT cs dcsdr phase visc thcon" << std::endl;
  for (int i=0;i<numygrid;i++) {
    double y = etu.rptree.ybounds[0]+(etu.rptree.ybounds[1]-etu.rptree.ybounds[0])*i/(numygrid-1.);
    for (int j=0;j<numxgrid;j++) {
      double x = etu.rptree.xbounds[0]+(etu.rptree.xbounds[1]-etu.rptree.xbounds[0])*j/(numxgrid-1.);
      eostableutrievaluatec0(&etu,1,1,&x,&y,f,type,clip,scratch,iscratch);
      std::cout << " " << x << " " << y << " " << f->rho[0] << " " << f->T[0]
		<< " " << f->p[0] << " " << f->E[0] << " " << f->S[0] << " " << f->F[0]
		<< " " << f->dpdr[0] << " " << f->dpdT[0]
		<< " " << f->dEdr[0] << " " << f->dEdT[0] << " " << f->cs[0] << " " << f->dcsdr[0]
		<< " " << f->phase[0] << " " << f->visc[0] << " " << f->thcon[0]
		<< std::endl;
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  //
  // print out the triangles
  //
  print_triangles(etu,type);

  //
  // cleanup
  //
  eos_table_utri_delete(&etu);
  eospointsfree(f);
  free(f);

  return 0;
}
