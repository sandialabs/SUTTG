/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SUTTG, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#include <limits>
#include <iostream>
#include <cmath>

#include "mpi.h"
#include "MCMin.H"

//
// Constructor
//
MCMin::MCMin(const bool master)
  : mpimaster(master), bestState(0), curState(0), nextState(0), widths(0),
    bestCost(0.), bestCostStop(-9.e99), curCost(0.), costTemp(1.), paramTemp(1.),
    paramTempMult(1.01), paramTempMax(1.e10), paramControl(1.), accTarget(0.2),
    accCounts(0), curAccBin(0), steps(0), maxsteps(1000), accsize(100)
{

}

//
// Destructor
//
MCMin::~MCMin()
{

}

//
// Minimization method
//
double MCMin::minimize( std::vector<double> & state )
{
  //
  // initial state evaluation
  //
  double initcost = f(state);

  //
  // initialize states
  //
  nextState = bestState = curState = state;
  bestCost = curCost = initcost;
  widths.assign(curState.size(),1.);

  if (mpimaster) {
    //
    // initialize acceptance array
    //
    curAccBin = 0;
    int startAcc = int(accsize*accTarget);
    accCounts.resize(accsize);
    for (int i=0;i<accsize;i++) accCounts[i] = (accsize-1-i)*startAcc/accsize;

    //
    // initialize everything else
    //
    steps = 0;
    paramControl = 1.;

  }

  int done = 0;
  while(1) {
    //
    // for skipping calculation if state is known to be invalid
    //
    int badstates = 0;

    if (mpimaster) {
      steps++;

      //
      // update widths for state generation
      //
      updateWidths(widths);

      //
      // give up on generating a proper state after 10 tries
      //
      while (badstates < 10000) {
        badstates++;
        generateState();
        if (checkState(nextState)) break;
      }
    }

    //
    // communicate new state
    //
    MPI_Bcast(&badstates,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(nextState.data(),nextState.size(),MPI_DOUBLE,0,MPI_COMM_WORLD);

    //
    // calculate the new state value
    //
    double cost = std::numeric_limits<double>::max();
    if (badstates < 10000) {
      cost = f(nextState);
    }

    if (mpimaster) {
      //
      // update best
      //
      if (cost < bestCost) {
        bestCost = cost;
        bestState = nextState;
      }

      //
      // update current state
      //
      bool accepted = false;
      if (cost < curCost-costTemp*log(mydrand()) && cost != curCost) {
        curCost = cost;
        curState = nextState;
        accepted = true;
      }

      //
      // acceptance rate control and annealing
      //
      updateAcceptance(accepted);
      if (paramControl < 1.) paramTemp *= paramTempMult;
      if (paramControl <= 0.01) paramTemp *= paramTempMult;

      //
      // stop
      //
      if (steps >= maxsteps || paramTemp > paramTempMax || bestCost < bestCostStop) done = 1;
      else done = 0;
    }

    //
    // everyone exit loop when done
    //
    MPI_Bcast(&done,1,MPI_INT,0,MPI_COMM_WORLD);
    if (done != 0) break;
  }

  //
  // communicate the final results
  //
  MPI_Bcast(&bestCost,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(bestState.data(),bestState.size(),MPI_DOUBLE,0,MPI_COMM_WORLD);

  state = bestState;
  return bestCost;
}

void MCMin::updateAcceptance( const bool accepted )
{
  //
  // update acceptance count array
  //
  if (accepted) for (std::size_t i=0;i<accCounts.size();i++) accCounts[i]++;

  //
  // modify control variable if needed
  //
  int tot = accCounts[curAccBin];
  int target = int(accTarget*accsize);
  if (accepted) {
    if (tot > target*2.0) paramControl *= 2.;
    else if (tot > target*1.4) paramControl *= 1.1;
    else if (tot > target*1.1) paramControl *= 1.01;
    if (paramControl > 1.) paramControl = 1.;
  }
  else {
    if (tot < target/2.0) paramControl /= 2.;
    else if (tot < target/1.4) paramControl /= 1.1;
    else if (tot < target/1.1) paramControl /= 1.01;
    if (paramControl < 0.01) paramControl = 0.01;
  }

  accCounts[curAccBin] = 0;
  curAccBin = (curAccBin+1)%accsize;

}

void MCMin::generateState()
{
  for (std::size_t i=0;i<curState.size();i++) {
    // ASA generator
    double x = mydrand();
    double y = -1.;
    if (x < 0.5) y = 1.;
    nextState[i] = curState[i]+widths[i]*y/paramTemp*(pow(1.+paramTemp,fabs(2.*x-1.))-1.);
  }
}
