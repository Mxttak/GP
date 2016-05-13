/************************************************************************
*
*  lap.h
   version 1.0 - 21 june 1996
   author  Roy Jonker, MagicLogic Optimization Inc.
   
   header file for LAP
*
**************************************************************************/

#include "stdio.h"

/*************** CONSTANTS  *******************/

  #define BIG 100000

/*************** TYPES      *******************/

  typedef int row;
  typedef int col;
  typedef double cost;

/*************** FUNCTIONS  *******************/

//extern cost lap(int dim, cost **assigncost,
//             row *rowsol, col *colsol, cost *u, cost *v);

cost lap(unsigned int dim, 
        cost **assigncost,
        col *rowsol, 
        row *colsol, 
        cost *u, 
        cost *v);

extern void checklap(int dim, int **assigncost,
                     int *rowsol, int *colsol, int *u, int *v);

