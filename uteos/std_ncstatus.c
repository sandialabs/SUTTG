/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SUTTG, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "netcdf.h"

#include "ncstatus.h"

/*
 * If a bad netcdf status is detected, stdout is flushed and an error
 * message is printed on stderr.
 *
 */
int bad_netcdf_status( int status,
                       char *ncfun,
                       const char *database,
                       char *fun,
                       char *file,
                       int line )
{
  if (status != 0) {
    fflush(stdout);
    fprintf(stderr,"NetCDF error %d \"%s\" calling: %s, database: %s, function: %s, file: %s, line: %d\n",
            status,nc_strerror(status),ncfun,database,fun,file,line);
  }

  return status;
}
