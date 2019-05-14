/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#ifndef DEOMCONST_H_
#define DEOMCONST_H_

#include "armadillo"

using namespace arma;

static const double    deom_pi = datum::pi;
static const cx_double deom_ci = cx_double(0.0,1.0); 
static const cx_double deom_c1 = cx_double(1.0,0.0); 

// static const double    deom_cm2unit = 4.5554927e-6;
// static const double    deom_kt2unit = 3.1662e-6;
// static const double    deom_fs2unit = 41.34902;

static const double    deom_cm2unit = 1.0;
static const double    deom_kt2unit = 1.0;
static const double    deom_fs2unit = 1.0;

#endif
