/*
 * File:        lsm_file.h
 * Copyright:   (c) 2005-2007 Masa Prodanovic and Kevin T. Chu
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Routines for checking if file management (zipping etc.)
 */

#ifndef included_lsm_file_h
#define included_lsm_file_h

#include "LSMLIB_config.h"

#define NO_ZIP 0
#define GZIP   1
#define BZIP2  2

#ifdef __cplusplus
extern "C" {
#endif


/*! 
 * checkUnzipFile() checks if file has .gz or .bz2 extention and 
 *   uncompresses the file.
 *
 *   Arguments:
 *    - file_name(in):     name of the file
 *    - pzip_status(out):  integer pointer that points to status of the
 *                         file (NO_ZIP, GZIP, BZIP2)
 *    - pfile_base (out):  pointer to string containing the file name base 
 *                         (file_name of the uncompressed file, same as the 
 *                         file name if there was no compression to begin with).
 *
 *  Return value:   none
 *
 *  Notes: 
 *     - There is no attempt to open the file.
 *     - file_base is allocated to 256 bytes within the function, and needs 
 *       to be freed later on
 *
 */        
void checkUnzipFile(char *file_name,int *pzip_status,char **pfile_base);

/*! zipFile() compresses the file according to its status and frees
 *   file_base pointer.
 *
 *   Arguments:
 *    - file_base(in):     name of the uncompressed file
 *    - zip_status(in):    integer compression status of the file 
 *                        (NO_ZIP,GZIP,BZIP2)
 *
 *  Return value:   none
 * 
 *  Notes: In case the status is NO_ZIP, function doesn't do anything.
 */        
void   zipFile(char *file_base,int zip_status);

#ifdef __cplusplus
}
#endif

#endif
