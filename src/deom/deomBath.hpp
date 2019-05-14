/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#ifndef DEOMBATH_H_
#define DEOMBATH_H_

#include <map>
#include <string>
#include "armadillo"
#include "JsonParser.hpp"
#include "deomConst.hpp"

using namespace std;
using namespace arma;
using namespace json11;

class bath {

    public:

    vec temperature;
    ivec   modLabel;
    cx_vec coef_lft;
    cx_vec coef_rht;
    vec    coef_abs;
    cx_vec expn_gam;
    vec    delt_res;

    bath (const Json& json) 
    {
        printf ("$InitBath\n");
        temperature = json2vec(json["temp"]);
        modLabel = json2ivec(json["mode"]);
        coef_lft = json2zvec(json["etal"]);
        coef_rht = json2zvec(json["etar"]);
        coef_abs = json2vec(json["etaa"]);
        expn_gam = json2zvec(json["expn"]);
        delt_res = json2vec(json["delr"]);
        temperature.print("temperature");
        modLabel.print("modLabel");
        coef_lft.print("coef_lft");
        coef_rht.print("coef_rht");
        coef_abs.print("coef_abs");
        expn_gam.print("expn_gam");
        delt_res.print("delt_res");
        printf ("$InitBath\n\n");
    }

    bath (const bath& rhs): 
          temperature(rhs.temperature),
          modLabel(rhs.modLabel),
          coef_lft(rhs.coef_lft), 
          coef_rht(rhs.coef_rht), 
          coef_abs(rhs.coef_abs), 
          expn_gam(rhs.expn_gam), 
          delt_res(rhs.delt_res) {}

    bath (const vec& temp, const ivec& mlbl, const cx_vec& etal, const cx_vec& etar, const vec& etaa, 
          const cx_vec& expn, const vec& delr): 
          temperature(temp),
          modLabel(mlbl),
          coef_lft(etal), 
          coef_rht(etar), 
          coef_abs(etaa),
          expn_gam(expn), 
          delt_res(delr) {}

    bath& operator= (const bath& rhs) {
        if (this != &rhs) {
            temperature = rhs.temperature;
            modLabel = rhs.modLabel;
            coef_lft = rhs.coef_lft;
            coef_rht = rhs.coef_rht;
            coef_abs = rhs.coef_abs;
            expn_gam = rhs.expn_gam;
            delt_res = rhs.delt_res;
        }
        return *this;
    }

   ~bath () {};

};

#endif
