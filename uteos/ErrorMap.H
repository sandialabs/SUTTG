/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SUTTG, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#ifndef ERROR_MAP_H
#define ERROR_MAP_H

#include <stdexcept>
#include <cmath>
#include <vector>
#include <iostream>

#include "EOSData.H"

class ErrorMap
{
public:

  //
  // Constructor
  //
  ErrorMap( const double tolerance,
	    const double ftolerance,
	    const std::vector<int> variables )
    : xlocation(1,0.0), ylocation(1,0.0), xrange(1,0.0), yrange(1,0.0),
      tolerances(1,tolerance), ftolerances(1,ftolerance), diffvars(1,variables)
  {
    if (tolerances[0] <= 0.0)
      throw std::runtime_error("ErrorMap: tolerance must be positive");
    if (ftolerances[0] < 0.0)
      throw std::runtime_error("ErrorMap: ftolerance must be non-negative");
    for (std::size_t i=0;i<variables.size();i++) {
      if (diffvars[0][i] < 0)
	throw std::runtime_error("ErrorMap: variable indices must be non-negative");
    }
  }

  //
  // Destructor
  //
  ~ErrorMap() {};

  //
  // get the global tolerance
  //
  double getTolerance() const { return tolerances[0]; }

  //
  // Compute the error for an EOSData structure. This is the maximum
  // over all inputs and specified outputs.
  //
  double compute_error( EOSData & p0,
			EOSData & p,
			const int debug = 0 ) const
  {
    double toterr(0.);

    if (debug > 0) {
      std::cout << " compute_error: " << std::endl;
      for (std::size_t i=0;i<xlocation.size();i++) {
        std::cout << " i " << i << " xl " << xlocation[i] << " yl " << ylocation[i]
                  << " xr " << xrange[i] << " yr " << yrange[i]
                  << " tol " << tolerances[i] << " ftol " << ftolerances[i] << " ints";
        for (std::size_t j=0;j<diffvars[i].size();j++) std::cout << " " << diffvars[i][j];
        std::cout << std::endl;
      }
    }

    //
    // find the appropriate error map entry
    //
    int mapi = 0;
    for (std::size_t i=1;i<xlocation.size();i++) {
      if (debug > 0) std::cout << "compute_error diffs: " << fabs(p0.inputs[0]-xlocation[i]) << " " << fabs(p0.inputs[1]-ylocation[i]) << std::endl;
      if (fabs(p0.inputs[0]-xlocation[i]) < xrange[i] &&
	  fabs(p0.inputs[1]-ylocation[i]) < yrange[i])
	mapi = i;
    }

    //
    // Always check the inputs
    //
    for (std::size_t i=0;i<p0.inputs.size();i++) {
      double relerr = calcrelerr(p0.inputs[i],p.inputs[i],tolerances[mapi],ftolerances[mapi]);
      if (relerr > toterr) toterr = relerr;
      if (debug > 0) std::cerr << " ini " << i << " ve " << p0.inputs[i] << " vi " << p.inputs[i] << " relerr " << relerr << std::endl;
    }

    //
    // Check the desired outputs
    //
    for (std::size_t i=0;i<diffvars[mapi].size();i++) {
      int k = diffvars[mapi][i];
      double relerr = calcrelerr(p0.outputs[k],p.outputs[k],tolerances[mapi],ftolerances[mapi]);
      if (relerr > toterr) toterr = relerr;
      if (debug > 0) std::cerr << " outi " << k << " ve " << p0.outputs[k] << " vi " << p.outputs[k] << " relerr " << relerr << std::endl;
    }

    return toterr;
  }

  //
  // Add a new error region to the map
  //
  void add_error_region( const double xcenter,
			 const double ycenter,
			 const double xradius,
			 const double yradius,
			 const double tolerance,
			 const double ftolerance,
			 const std::vector<int> variables )
  {
    if (tolerance <= 0.0)
      throw std::runtime_error("ErrorMap::add_error_region: tolerance must be positive");
    if (ftolerance < 0.0)
      throw std::runtime_error("ErrorMap::add_error_region: ftolerance must be non-negative");
    if (xradius < 0.0 || yradius < 0.0)
      throw std::runtime_error("ErrorMap::add_error_region: radii must be non-negative");
    for (std::size_t i=0;i<variables.size();i++) {
      if (variables[i] < 0)
	throw std::runtime_error("ErrorMap: variable indices must be non-negative");
    }

    xlocation.push_back(xcenter);
    ylocation.push_back(ycenter);
    xrange.push_back(xradius);
    yrange.push_back(yradius);
    tolerances.push_back(tolerance);
    ftolerances.push_back(ftolerance);
    diffvars.push_back(variables);
  }

private:

  //
  // Method to compute the relative error referenced to 1 for a scalar
  //
  double calcrelerr( const double f_exact,
		     const double f_comp,
		     const double tol,
		     const double ftol ) const
  {
    double err = tol*fabs(f_exact);
    if (ftol > err) err = ftol;
    return fabs(f_comp-f_exact)/err;
  }

  std::vector<double> xlocation;
  std::vector<double> ylocation;
  std::vector<double> xrange;
  std::vector<double> yrange;
  std::vector<double> tolerances;
  std::vector<double> ftolerances;
  std::vector<std::vector<int> > diffvars;

};

#endif
