/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SUTTG, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#ifndef TABLES_H
#define TABLES_H

#include <vector>
#include <string>

#include "TMesh.H"
#include "utri_eos_support.h"
#include "RMesh.H"

class Tables
{
public:

  Tables( std::vector<Parameters> & model_params,
	  const std::string model_name,
	  const std::string mesh_type,
	  const std::vector<double> & xvar_bounds,
	  const std::vector<double> & yvar_bounds,
	  const int logxy,
	  const int totalthreads,
	  const int mpirank = 0,
	  const int mpisize = 1 );

  virtual ~Tables();
  
  virtual void mesh_tables( const double tolerance ) = 0;

  void write_tables( const std::string basename,
                     const std::vector<double> & defaults );

protected:

  void create_regions_PX();

  std::vector<TriMesh> rtmeshes;

  TriMesh::MeshTypes meshtype;

  std::vector<double> Xbounds;
  std::vector<double> Ybounds;

  std::vector<std::vector<TMesh *> > rmst;
  std::vector<std::vector<TMesh *> > rmse;

  int mintab,maxtab;
  bool mpimaster;
  int myrank;
  int mpicommsize;
  int logvars;
  int totthr;

private:

  void write_table( const std::string filename,
                    std::vector<TMesh *> & rms,
                    const TriMesh::MeshTypes mtype,
                    const std::vector<double> & defaults );

  void set_eosdatapoint( eos_data_points_t *d,
                         int p,
                         TriMesh::MeshTypes mtype,
                         Node & n );

};

#endif
