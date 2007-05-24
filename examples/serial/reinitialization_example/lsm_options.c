/*
 * File:        lsm_options.c
 * Copyright:   (c) 2006 Masa Prodanovic 
 * Revision:    $Revision: 1.1 $
 * Modified:    $Date: 2006/06/02 03:02:39 $
 * Description: Implementation file for Options structure.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "lsm_options.h"
#include "lsm_grid.h"


#define DSZ sizeof(double)

/*======== Helper Functions for Options structure manipulation ========*/
/*!
 * allocateOptions() allocates memory for an Options structure
 *
 * Arguments: none
 *
 * Return value: a pointer to the allocated memory
 *
 */
Options *allocateOptions(void);




/*================= Options structure manipulation ==================*/

Options *allocateOptions(void)
{
  Options *options = (Options *)calloc(1,sizeof(Options));

  return  options;
}


Options *createOptionsDefault(void)
{
  Options *options;
  
  options = allocateOptions();
  
  /* Main part of the structure */

  sprintf(options->outfile,"out_file");
  sprintf(options->path,"./");

  options->dx = 0.02;
  options->tmax = 20.0;
  sprintf(options->accuracy,"medium");
  options->accuracy_id = MEDIUM;

  options->a = 0.1;
  options->b = 0.05;

  options->save_data = 1;
  options->do_reinit = 1;  
  options->do_mask = 0;

  /* User additions */
  
  /* end User additions */
    
  return options;
}


Options *copyOptions(Options *options_src)
{
  Options *options;
  
  options = allocateOptions();
  
  /* Main part of the structure */

  sprintf(options->outfile,options_src->outfile);
  sprintf(options->path,options_src->path);

  options->dx = options_src->dx;
  options->tmax = options_src->tmax;
  sprintf(options->accuracy,options_src->accuracy);
  options->accuracy_id = options_src->accuracy_id;

  options->a = options_src->a;
  options->b = options_src->b;

  options->save_data = options_src->save_data;
  options->do_reinit = options_src->do_reinit;
  options->do_mask = options_src->do_mask;	

  /* User additions */
    
  /* end User additions */
    
  return options;
}

Options *createOptionsFromInputFile(char *filename)
{
  FILE    *fp;
  char    line[100], c, word[100], *exist, *new;
  Options *options;

  int     i, n, found, len, tmp1, cmp;
  double  tmp;

  int N_accur_menu = sizeof(Accuracy_settings_menu)/
                     sizeof(Accuracy_settings_item);

  options = createOptionsDefault();

  if( (fp = fopen(filename,"r")) == NULL )
  {
     printf("\nCouldnt open file %s",filename);
     return options;
  }

  while( fgets(line,80,fp) != NULL )
  {
    /* skip empty spaces at the beginning of the line */
    n = 0;
    while(isspace(line[n])) n++;

    c = tolower(line[n]);

     /* Main part of the structure */

    if ( c == 'a' )
    { /* could be 'a' or 'accuracy' */
      if ( tolower(line[n+1]) == 'c' )
      {
        sscanf(line+n,"accuracy %s",word);
    
        /* try to identify input with one of the accuracy menu items */
        found = 0;
        for(i = 0; i < N_accur_menu; i++)
        {
          cmp = strcmp(word,Accuracy_settings_menu[i].name);
          if( cmp == 0 )
          {
            found = 1;
            options->accuracy_id = i;
            sprintf(options->accuracy,"%s",Accuracy_settings_menu[i].name); 
          }   
        }

        if( !found )
        {
          printf("\nAccuracy type %s not found, set to default.",word);
        }
      }
      else
      {
        sscanf(line+n,"%*s %lf ",&tmp);
        options->a = tmp;
      }
    }
    else if( c == 'b' )
    {  /* 'b' */
      sscanf(line+n,"%*s %lf ",&tmp);
      options->b = tmp; 
    }
    else if( c == 't' )
    { /* tmax */
      sscanf(line+n,"%*s %lf ",&tmp);
      options->tmax = tmp;
    }
    else if( c == 'd' )
    { /* could be 'dx' or 'do_reinit' */
      if ( tolower(line[n+1]) == 'x' )
      {
        sscanf(line+n,"%*s %lf ",&tmp);
        options->dx = tmp;
      }
      else if( line[n+3] == 'r')
      {
        sscanf(line+n,"%*s %d ",&tmp1);
        if ( (tmp1 == 0) || (tmp1 == 1))
          options->do_reinit = tmp1;
        else
        {
          printf("\nIncorrect do_reinit option %d, set to default.\n",tmp1);
        }    
      }
    }
    else if( c == 'o' )
    { /* 'outfile' */
      sscanf(line+n,"outfile %s",word);
      sprintf(options->outfile,"%s",word);

      /* test filename 
      fp1 = fopen(options->outfile,"w");
      if( fp1 == NULL )
      {
        printf("\nCannot open file %s, outfile set to default.",
               options->outfile);
        sprintf(options->outfile,"out");
      }
      else fclose(fp1); */


      /* set path to the one given in outfile name */
      exist = strrchr(options->outfile,'/');
      if( exist )
      {
        len = strlen(options->outfile);
        for(i = len-1; i >= 0; i--)
        {
          if( (options->outfile)[i] == '/' ) break;
        }
        new  = malloc(i+1);
        strncpy(new,options->outfile,i+1);
        sprintf(options->path,"%s",new);
      }
      else sprintf(options->path,"./");
    }
    else if( c == 's' && (tolower(line[n+1]) == 'a') )
    { /* 'save_data' */  
      sscanf(line+n,"%*s %d ",&tmp1);
      if ( (tmp1 == 0) || (tmp1 == 1))
        options->save_data = tmp1;
      else
      {
        printf("\nIncorrect save_data option %d, set to default.\n",tmp1);
      }
    }


    /* User additions */
    /* Watch for overlaps and words starting with the same letter */
    /* end User additions */
  }

  fclose(fp);
  return options;
}


void  printOptions(Options *options,FILE *fp)
{
  /* Main part of the structure */

  fprintf(fp,"\nOptions:\n");
  fprintf(fp,"  outfile   %8s  [ output file name ]\n",options->outfile);
  /* no need to print path - used internally */

  fprintf(fp,"  dx        %8g  [ grid spacing ]\n",options->dx);
  fprintf(fp,"  tmax      %8g  [ max running time allowed ]\n",options->tmax);
  fprintf(fp,"  accuracy  %8s  [ accuracy options: low, medium, high, very_high ]\n",options->accuracy);

  fprintf(fp,"  a         %8g  [ constant; motion in normal direction ]\n",
                                                                    options->a);
  fprintf(fp,"  b         %8g  [ constant; motion by mean curvature   ]\n",
                                                                    options->b);                       

  fprintf(fp,"  save_data %8d  [ save data (1) or not (0)       ]\n",
                                                            options->save_data);
  fprintf(fp,"  do_reinit %8d  [ reinitilize periodically (1) or not (0)]\n",
                                                            options->do_reinit);
  fprintf(fp,"  do_mask   %8d  [ impose mask (1) or not (0)]\n",
                                                              options->do_mask);							    

  /* User additions */
  /* end User additions */
}


void destroyOptions(Options *options)
{
  if(options) free(options);
}
