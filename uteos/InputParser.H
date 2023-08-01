/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SUTTG, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#ifndef INPUT_PARSER_H
#define INPUT_PARSER_H

#include "XMLElement.h"
#include "EOSFactory.H"
#include "EOSParam.H"

class SharedParamData
{
public:

  //
  // Constructor
  //
  SharedParamData( const double param_value = 0.0 )
    : value(param_value) {}

  ~SharedParamData() {}

  double value;

  // vectors denoting name of the model, param location, and scale
  // factor to which this shared param is copied
  std::vector< std::string > names;
  std::vector< int >         plocs;
  std::vector< double >      scales;

};

class InputParser
{
public:
  //
  // Constructor
  //
  InputParser( const std::string XMLinput,
	       const int verbose = 0 );

  //
  // Destructor
  //
  ~InputParser() {}

  void getTableStrings( std::string & type,
			std::string & outname ) {
    type = tabtype;
    outname = outputfilebase;
  }

  const std::string & getMeshVars() const { return meshvars; }

  int getLogVars() const { return logvars; }

  int getNumThreads() const { return numthreads; }

  Parameters & getParameters() { return p; }

  const std::string & getModel() const { return runmodel; }

  void getBounds( std::vector<double> & Xb,
                  std::vector<double> & Yb ) {
    Xb = Xbounds;
    Yb = Ybounds;
  }

  void getTolerances( double & tol,
                      double & ftol ) {
    tol = tolerance;
    ftol = ftolerance;
  }

  void getSamples( int & boundary,
                   int & region ) {
    boundary = bsamples;
    region = rsamples;
  }

  void getNodeOptimizer( std::string & alg,
                         int & anodes,
                         double & errlev,
                         double & mult ) {
    alg = nodealgorithm;
    anodes = addnodes;
    errlev = steperrorlevel;
    mult = stepmultiplier;
  }

  void getDefaults( std::vector<double> & Defaults ) {
    Defaults[0]= atomic_number;
    Defaults[1]= atomic_weight;
    Defaults[2]= rho_ref;
    Defaults[3]= t_ref;
    Defaults[4]= 0.; //To be filled later as bulk modulus for some table formats
  }

  void copySharedParam( const std::string pname,
                        const double val,
                        Parameters & params );

  int getVerbosity() const { return verbosity; }

private:

  //
  // the parsed XML input file
  //
  RefPtr<XMLElement> xml_tree;

  //
  // the type of tabulation
  //
  std::string tabtype;

  //
  // base name for output table names
  //
  std::string outputfilebase;

  //
  // the independent variables for the mesh
  //
  std::string meshvars;

  //
  // whether to log scale the independent variables
  //
  int logvars;

  //
  // number of threads for state evaluation in mesh error calculation
  //
  int numthreads;
  
  // Reference values
  double atomic_number;
  double atomic_weight;
  double rho_ref;
  double t_ref;

  //
  // lower and upper bounds for table
  //
  std::vector<double> Xbounds;
  std::vector<double> Ybounds;

  //
  // tolerances for table generation and optimization
  //
  double tolerance;
  double ftolerance;

  //
  // error sample size for table generation and optimization
  //
  int bsamples;
  int rsamples;

  //
  // number of nodes to add each time
  //
  int addnodes;

  //
  // algorithm and options for node optimization
  //
  std::string nodealgorithm;
  double steperrorlevel;
  double stepmultiplier;

  //
  // model parameters
  //
  Parameters p;

  //
  // shared parameter data
  //
  std::map< std::string, SharedParamData > sp;

  //
  // name of the model to use
  //
  std::string runmodel;

  //
  // method to setup a shared parameter
  //
  double getSharedParam( RefPtr<XMLElement> & mi,
			 const std::string name,
			 const std::string param_name,
			 const int param_loc,
			 const double param_scale );

  //
  // method to parse the EOSModel input section, creating the Parameters variable
  //
  void parseEOSModel( const std::string mname,
		      RefPtr<XMLElement> & mi );

  //
  // method to parse the Tabulation input section
  //
  void parseTabulation( RefPtr<XMLElement> & mi );

  //
  // parser helper
  //
  EOSFactory eosfactory;

  //
  // flag for debugging info
  //
  int debug;

  //
  // flag for debugging table generation
  //
  int verbosity;

};

#endif
