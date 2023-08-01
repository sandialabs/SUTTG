/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SUTTG, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#ifndef UTRI_TABLES_H
#define UTRI_TABLES_H

#include "Tables.H"
#include "MCMin.H"
#include "RMesh.H"

class UTriMMOpt : public MCMin
{
public:

  UTriMMOpt( const bool master,
             std::vector<RMesh *> & regionmeshes,
             RMeshMMOptData & meshdata,
             const int logscale,
             const int mpirank,
             const int mpisize,
             const int region_samples );

  ~UTriMMOpt();

  double f( const std::vector<double> & state );

  void getInitState( std::vector<double> & state );

  void updateWidths( std::vector<double> & widths );

  bool checkState( const std::vector<double> & state );

  void printMeshState( const std::vector<double> & state,
                       std::ostream & ofile = std::cout );

private:

  std::vector<RMesh *> rmeshes;
  RMeshMMOptData & data;
  std::vector<Tri> reftris;
  int logvars;
  int myrank,mpicommsize;
  int rsamples;
  int debug;

  std::map<int,int> boundarymap;
  int nbp;

};

class UTriTables : public Tables
{
public:

  UTriTables( std::vector<Parameters> & model_params,
	      const std::string model_name,
	      const std::string mesh_type,
	      const std::vector<double> & xvar_bounds,
	      const std::vector<double> & yvar_bounds,
	      const int logxy,
	      const int totalthreads,
	      const int boundary_samples,
	      const int region_samples,
              const std::string & node_algorithm,
	      const int add_nodes,
              const double step_multiplier,
              const double step_error_level,
	      const int mpirank = 0,
	      const int mpisize = 1,
	      const int debug_level = 0 );

  ~UTriTables();

  void mesh_tables( const double tolerance ) {
    mesh_tables_PX(tolerance);
  }

private:

  void mesh_tables_PX( const double tolerance );

  double refine_regions( std::vector<std::vector<TMesh *> > & rms,
                         int region = -1 );

  double refine_regions_pt( std::vector<std::vector<TMesh *> > & rms,
                            int region = -1 );

  double refine_boundaries( std::vector<std::vector<BMesh *> > & bms,
                            bool uniformerr = true,
                            int boundary = -1,
                            int debug = 0 );

  double setnodelocs( std::vector<std::vector<BMesh *> > & bms,
                      const int b,
                      double rtol,
                      int nnodes );

  double compute_region_error( std::vector<RMesh *> & rmeshes,
                               std::vector<int> & tnodes,
                               std::vector<std::pair<int,std::vector<double> > > & maxerrtris );

  double compute_region_error_initial( std::vector<RMesh *> & rmeshes,
                                       std::vector<std::pair<int,double> > & maxerrtris );

  double compute_region_error_initial( std::vector<RMesh *> & rmeshes,
                                       std::vector<int> & tnodes,
                                       std::vector<std::pair<int,double> > & maxerrtris );

  double optimize_nodes( std::vector<RMesh *> & rmeshes,
			 RMeshMMOptData & utmmod,
			 std::vector<std::pair<int,double> > & trierrs,
                         const bool accurate,
                         const double lastmaxerr );

  void setup_opt( std::vector<RMesh *> & rmeshes,
                  const int curtri,
                  std::vector<int> & lasttrinodes,
                  RMeshMMOptData & rmmmod );

  std::vector<std::vector<BMesh *> > bmst;
  std::vector<std::vector<BMesh *> > bmse;

  int bsamples;
  int rsamples;

  std::string nodealg;
  int addnodes;
  double stepmult;
  double steperr;

  int debug;

};

#endif
