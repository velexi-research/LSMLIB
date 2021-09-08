/*
 * File:        lsm_file.c
 * Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
 *                  Regents of the University of Texas.  All rights reserved.
 *              (c) 2009 Kevin T. Chu.  All rights reserved.
 * Revision:    $Revision$
 * Modified:    $Date$
 * Description: Implementation file for file management functions (zipping etc.)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lsmlib/lsmlib_config.h"
#include "lsmlib/lsm_file.h"


void checkUnzipFile(char *file_name,int *pzip_status, char **pfile_base)
{
     int     zip_status = (int)NO_ZIP, length_base;
     char    command[256], *file_base, *gptr, *bptr;
     
     file_base = (char *)malloc(256*sizeof(char));
     
     gptr = strstr(file_name,".gz");
     bptr = strstr(file_name,".bz2");
     if( gptr != (char *)NULL )
     {
       zip_status = GZIP;
       sprintf(command,"gunzip -f %s",file_name);
       system(command);
       //printf("\n%s",command);
       
       length_base = (strlen(file_name)-3);       
       file_base = strncpy(file_base,file_name,length_base);
       file_base[length_base] = '\0';
       //printf("\n%s",file_base); fflush(stdout); fflush(stderr);
     }
     else if( bptr!= (char *)NULL )
     {
       zip_status = BZIP2;
       sprintf(command,"bunzip2 -f %s",file_name);
       system(command);
       length_base =strlen(file_name)-4;
       file_base = strncpy(file_base,file_name,length_base);
       file_base[length_base] = '\0';
     }
     else 
       file_base = strcpy(file_base,file_name);
     
    *pzip_status = zip_status;
    *pfile_base = file_base;
}

void  zipFile(char *file_base,int zip_status)
{
     char command[256];
     
     if( zip_status == GZIP ) 
     {
        sprintf(command,"gzip -f %s",file_base);
        system(command);
     }
     else if( zip_status == BZIP2 ) 
     {
        sprintf(command,"bzip2 -f %s",file_base);
        system(command);
     }
}
