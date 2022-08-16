/*
**********************************************************
impacta_headers

version 1.8

agrt
11/2/07 big change - now put in ib heating into f2 equation

30/8/07 build file/makefile and header changed so that compiling can be with or without xnbody

16/8/07 this version with visit added 2011 - removed, obsolete

17/7/07 
this version fixed for namespace ambiguity which showed up when compiling on
linux machine

15/4/07
this version fixed for memory leaks

8/3/07

important notes:
inputting the equations is unfortunately going to be a nightmare:

for efficiency reasons each equation x must have a number of functions 
written for them. these are:
__________________________________________________________
xequation_inner, xequation_outer

these are  codes which will be almost identical for putting in the 
appropriate matrix elements. inner code has no boundary checking to make
it quicker - outer code is the same but checks the boundary elements
__________________________________________________________
xequation_count 

this code looks for non-zeros only in the elements which should have
non-zeros - i.e. those given by the above equations.
__________________________________________________________
xequation_clean_inner,xequation_clean_outer

these zero the appropriate elements in the vectorfor next use.

should point out that f0 and f1 equations will need to be 
updated when f2 and f3 added.
__________________________________________________________
xequation_pack_inner,xequation_pack_outer

these pack the sparse matrix efficiently.
should point out that f0 and f1 equations will need to be 
updated when f2 and f3 added.
**********************************************************


**********************************************************
21/3/07

mpi message codes - 

200 - gathering the whole vector for proc 0
300 exchanging ghost cells - lower boundary.
400          "             - upper boundary
700 - gathering a single k string from one processor

28/8/07 - have to be careful with mpi_barrier - 
i was using it rather loosely before.
*/
// to include xnbody or not
#undef impacta_xnbody_on
// whether to use ansi colour codes or not
#undef impactawithcolour

#ifndef __cplusplus
#error a c++ compiler is required!
#endif 

//libraries
#include <mpi.h>
#include <iostream>
#include <fstream>

#include <petscksp.h>
#include <cmath>
#include <ctime>
#include <sstream>
#include <cstring>
#include <cstdlib>

//setup
#include "impacta_definitions.h"
#include "impacta_environment.h"

using namespace globalconsts;
bool IMPACTA_Version_HPC=true;

const bool if_time_messages = false;

//objects
#include "impacta_objects/impacta_dims.h"
//using namespace with_vec_err_checking;
using namespace no_vec_err_checking;
#include "impacta_objects/impacta_vector.h"
#include "impacta_objects/impacta_config.h"
#include "impacta_objects/impacta_parallel.h"
#include "impacta_objects/impacta_matrix.h"
#include "impacta_objects/impacta_parallel_sparse.h"


//operators
#include "impacta_operators/impacta_int_v.h"
#include "impacta_operators/impacta_stenops.h"
using namespace differentials_x_c_y_c; // center difference x and y
//using namespace differentials_x_t_y_t; // test difference x and y
#include "impacta_operators/impacta_variables.h"
//#include "impacta_operators/impacta_variables_debug.h"

// new ! function pointers for boundary conditions
void (*IMPACT_f0_bound)(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus);
void (*IMPACT_f1_E_bound)(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1);
void (*IMPACT_B_bound)(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1);
void (*IMPACT_f2_bound)(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1, IMPACT_Dim *x2);
void (*IMPACT_ni_bound)(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus);
void (*IMPACT_Ci_bound)(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1);
char boundaries[2][11];
int Boundaries_open_x;
int Boundaries_open_y;

#include "impacta_operators/impacta_boundaries.h"
//boundary conditions entered here:
/*
  these can be:
  p - periodic
  r - reflecting, bz flips across boundary
  o - open boundaries
  b - reflecting, bz constant over boundary
 */
//using namespace BC_x_r_y_p; -> now in main program


#include "impacta_operators/impacta_operators.h"

//communication between processors..
#include "impacta_functions/impacta_swapghosts.h"

//exit
#include "impacta_functions/impacta_exit.h"

//sparse matrix functions.
#include "impacta_functions/impacta_sparse_count.h"
#include "impacta_functions/impacta_cleanvector.h"
using namespace no_sparserow_checking;
//using namespace with_sparserow_checking;
#include "impacta_functions/impacta_fast_sparse_pack.h"
using namespace no_zero_vector_checking; //check v-zeroing is correct
//using namespace with_zero_vector_checking;
//the checks really slow it down!

//initialization and input/output
#include "impacta_io/impacta_init.h"
#include "impacta_io/impacta_moments.h"
#include "impacta_functions/impacta_moment_funcs.h"
#include "impacta_io/impacta_dists.h"
#include "impacta_io/impacta_rw_moments.h"
#include "impacta_io/impacta_filehandling.h"

#include "impacta_io/impacta_input_maths.h"
#include "impacta_io/impacta_input_functions.h"
#include "impacta_io/impacta_input_functions_2d.h"
#include "impacta_io/impacta_input.h"

//equations
#include "impacta_equations/impacta_collisions.h"
#include "impacta_equations/impacta_f0equation.h"
#include "impacta_equations/impacta_f1equation.h"
#include "impacta_equations/impacta_f2equation.h"
#include "impacta_equations/impacta_eequation.h"
#include "impacta_equations/impacta_bequation.h"
#include "impacta_equations/impacta_ionmotion.h"
#include "impacta_equations/impacta_ionization.h"

//functions for forming and checking sparse matrix
#include "impacta_functions/impacta_form_sparse.h"
#include "impacta_functions/impacta_switch_elements.h"
#include "impacta_functions/impacta_explicit.h" 
//forms explicit sparse matrix - reintroduced (25/10/07)

//#include "impacta_functions/impacta_debugging_code.h"
#include "impacta_functions/impacta_lag_check.h"

//matrix solver
//#include "impacta_roughmatrixsolver.h"
#include "impacta_extpks/impacta_matrix_solver.h"
//#include "impacta_extpks/xnbody_funcs.h"

// messages and main functions
#include "impacta_io/impacta_startmessages.h"
#include "impacta_help.h"
#include "impacta_functions/impacta_main_funcs.h"

// Ray Tracing package
#include "impacta_raytrace/matrix_ops.h"
#include "impacta_raytrace/vector_ops.h"
//#include <boost/math/special_functions/erf.hpp>

#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/math/special_functions/airy.hpp>
#include <boost/date_time/posix_time/posix_time.hpp> 
#include <boost/date_time/gregorian/gregorian.hpp>
#include "impacta_raytrace/tracer.h"

