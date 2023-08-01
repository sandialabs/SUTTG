/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SUTTG, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

// buildtables.cpp
//
// Given input of an EOS model and associated sets of parameters,
// build a corresponding set of utri tables that represent the given
// model instances. The specified base table controls the default
// topology for all other tables, as determined by meeting the desired
// accuracy tolerance. If other model instances cannot conform to the
// topology, then the process will fail.

#include <vector>
#include <stdexcept>
#include <iostream>
#include <istream>
#include <sstream>
#include <cstdlib>
#include <deque>
#include <fstream>
#include <string>
#include <cerrno>

#include <mpi.h>

//
// EOS library defines
//
#include "EOSParam.H"
#include "EOSModel.H"
#include "EOSFactory.H"
#include "PhaseBoundaryInfo.H"

//
// build table defines
//
#include "Node.H"
#include "TriMesh.H"
#include "BMesh.H"
#include "RMesh.H"
#include "Tables.H"
#include "UTriTables.H"

//
// table output defines
//
#include "utri_eos_support.h"

//
// input parsing define
//
#include "InputParser.H"

int read_file_contents( const char *filename,
                        char *& s,
                        int & slen )
{
  std::ifstream in(filename, std::ios::in | std::ios::binary);
  if (in)
  {
    in.seekg(0, std::ios::end);
    slen = in.tellg();
    s = new char[slen];
    in.seekg(0, std::ios::beg);
    in.read(&s[0], slen);
    in.close();
  }
  return (errno);
}

int main( int argc,
          char * argv[] )
{
  int mpi_thread_level;
  int mpi_error = MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&mpi_thread_level);

  if (mpi_error != MPI_SUCCESS /*|| mpi_thread_level < MPI_THREAD_FUNNELED*/ ) {
    MPI_Finalize();
    std::ostringstream oss;
    oss << "MPI_Init_thread failed: error " << mpi_error << " thread " << mpi_thread_level << std::endl;
    throw std::runtime_error(oss.str());
  }
  //else {
  //  std::cout << MPI_Init_thread requested level " << MPI_THREAD_FUNNELED << " got " << mpi_thread_level << std::endl;
  //}

  int mpirank;
  int mpisize;
  mpi_error = MPI_Comm_rank(MPI_COMM_WORLD,&mpirank);
  mpi_error = MPI_Comm_size(MPI_COMM_WORLD,&mpisize);

  if (argc < 2 || argc > 2) {
    MPI_Finalize();
    if (mpirank == 0) throw std::runtime_error("Usage: utriTableBuilder <XMLInputFileName>");
    else return -1;
  }

  //
  // read in the file and bcast it
  //
  char * xml_infile = NULL;
  int xml_infile_length = 0;
  if (mpirank == 0) {
    //
    // get the input file
    //
    int readerr = read_file_contents(argv[1],xml_infile,xml_infile_length);
    if (readerr) {
      std::string estr("Failed to open input file "+std::string(argv[1]));
      perror(estr.c_str());
      MPI_Abort(MPI_COMM_WORLD,-1);
    }
  }
  MPI_Bcast(&xml_infile_length,1,MPI_INT,0,MPI_COMM_WORLD);
  if (mpirank != 0) xml_infile = new char[xml_infile_length];
  MPI_Bcast(xml_infile,xml_infile_length,MPI_CHAR,0,MPI_COMM_WORLD);

  //
  // make string to pass to parser
  //
  std::string xml_content(xml_infile,xml_infile_length);
  delete [] xml_infile;

  InputParser iparser(xml_content/*,1*/);

  // get desired debug level
  int debug = iparser.getVerbosity();

  std::string tabletype;
  std::string tablename;

  iparser.getTableStrings(tabletype,tablename);

  const std::string meshvars = iparser.getMeshVars();

  //
  // setup the parameters for each instance with a factory helper to determine location
  //
  std::vector<Parameters> model_params(1,iparser.getParameters());

  //
  // get which model to initialize out of parameter list
  //
  const std::string modelname = iparser.getModel();

  //
  // get bounds for the tables
  //
  std::vector<double> Xbounds;
  std::vector<double> Ybounds;
  iparser.getBounds(Xbounds,Ybounds);

  //
  // get table tolerances
  //
  double table_tolerance;
  double table_ftolerance;
  iparser.getTolerances(table_tolerance,table_ftolerance);

  //
  // get sample sizes
  //
  int boundary_samples;
  int region_samples;
  iparser.getSamples(boundary_samples,region_samples);

  //
  // get node optimizer info
  //
  std::string node_algorithm;
  int add_nodes;
  double steperrorlevel;
  double stepmultiplier;
  iparser.getNodeOptimizer(node_algorithm,add_nodes,steperrorlevel,stepmultiplier);

  //
  // get log scaling flag
  //
  const int logvars = iparser.getLogVars();

  //
  // get thread count
  //
  const int numthreads = iparser.getNumThreads();
  
  //
  // defaults vector
  //
  std::vector<double> defaults(5);
  iparser.getDefaults(defaults);

  //
  // seed system RNG after all calls to the parser, since libexpat uses the RNG
  //
  srandom(1);

  //
  // Create the tables
  //
  Tables * tab = 0;
  if(tabletype=="utri") {
    tab = new UTriTables(model_params,modelname,meshvars,Xbounds,Ybounds,
                         logvars,numthreads,boundary_samples,region_samples,
                         node_algorithm,add_nodes,stepmultiplier,steperrorlevel,
                         mpirank,mpisize,debug);
    tab->mesh_tables(table_tolerance);
  }
  else {
    std::cout << "Unknown output table type " << tabletype << " was requested." << std::endl;
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }

  //
  // output the tables
  //
  tab->write_tables(tablename,defaults);

  //
  // cleanup
  //
  if (tab != 0) delete tab;

  MPI_Finalize();

  return 0;

}
