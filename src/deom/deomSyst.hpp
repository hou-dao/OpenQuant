/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#ifndef DEOMSYST_H_
#define DEOMSYST_H_

#include <string>
#include "armadillo"
#include "JsonParser.hpp"

using namespace std;
using namespace arma;
using namespace json11;

class syst {
    
    public:
    
    int nsys;
    int nmod;
    cx_mat  ham1;
    cx_cube qmd1;
    
    syst (const Json& json) {
        printf ("$InitSyst\n");
        ham1 = json2zmat(json["ham1"]);
        qmd1 = json2zcube(json["qmd1"]);
        ham1.print("ham1");
        qmd1.print("qmd1");
        nsys = qmd1.n_rows;
        nmod = qmd1.n_slices;
        printf ("$InitSyst\n\n");
    }
    
    syst (const cx_mat& _h, const cx_cube& _q): nsys(_q.n_rows), nmod(_q.n_slices), ham1(_h), qmd1(_q) {}

    syst (const syst& _s): nsys(_s.nsys), nmod(_s.nmod), ham1(_s.ham1), qmd1(_s.qmd1) {}
        
   ~syst () {}
};

#endif
