/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SUTTG, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#ifndef MCMIN_H
#define MCMIN_H

#include <vector>
#include <cstdlib>

struct MCMinParams
{
  MCMinParams()
    : costTemp(1.),
      paramTempInit(1.),
      paramTempMult(1.01),
      paramTempMax(1.e10),
      accTarget(0.2),
      bestCostStop(-9.e99),
      maxsteps(1000),
      accsize(100)
  {}

  double costTemp;
  double paramTempInit;
  double paramTempMult;
  double paramTempMax;
  double accTarget;
  double bestCostStop;
  int maxsteps;
  int accsize;
};

class MCMin
{
public:

  //
  // Constructor
  //
  MCMin( const bool master );

  //
  // Destructor
  //
  virtual ~MCMin();

  //
  // Minimization method
  //
  double minimize( std::vector<double> & state );

  //
  // Initialize control constants
  //
  void initialize( const MCMinParams & p ) {
    costTemp = p.costTemp;
    paramTemp = p.paramTempInit;
    paramTempMult = p.paramTempMult;
    paramTempMax = p.paramTempMax;
    accTarget = p.accTarget;
    bestCostStop = p.bestCostStop;
    maxsteps = p.maxsteps;
    accsize = p.accsize;
  }

protected:

  //
  // Compute the function value at the given state
  //
  virtual double f( const std::vector<double> & state ) = 0;

  //
  // Update the sample widths
  //
  virtual void updateWidths( std::vector<double> & widths ) {
    widths.assign(widths.size(),1.);
  }

  //
  // Return true if the state is valid ( done before an evaluation )
  //
  virtual bool checkState( const std::vector<double> & /*state*/ ) { return true; }

  //
  // true if the master mpi process
  //
  bool mpimaster;

private:

  //
  // acceptance rate control
  //
  void updateAcceptance( const bool accepted );

  //
  // create a new state vector
  //
  void generateState();

  //
  // random number generator for uniform double in [0,1)
  //
  double mydrand()
  {
    return random()/(RAND_MAX+1.0);
  }

  std::vector<double> bestState;
  std::vector<double> curState;
  std::vector<double> nextState;
  std::vector<double> widths;
  double bestCost;
  double bestCostStop;
  double curCost;
  double costTemp;
  double paramTemp;
  double paramTempMult;
  double paramTempMax;
  double paramControl;
  double accTarget;
  std::vector<int> accCounts;
  int curAccBin;
  int steps;
  int maxsteps;
  int accsize;
};

#endif
