/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#ifndef DEOMPULSE_H_
#define DEOMPULSE_H_

#include <map>
#include <string>
#include <cstdio>
#include "armadillo"
#include "JsonParser.hpp"
#include "deomConst.hpp"

using namespace std;
using namespace arma;
using namespace json11;

/**
 * Chirped Gaussian pulse
 *
 * E(t) = \frac{ampl*\pi}{\sqrt{2\pi}sigm}
 *        \exp(-\frac{t^2}{2\sigm^2})
 *        \cos(freq*t+varp*t^2/2)
 *
 **/

class pulse {

    public:

    bool   on;
    double ampl;
    double sigm;
    double freq;
    double varp;

    pulse (const Json& json) {

        ampl = json["ampl"].number_value();
        sigm = json["sigm"].number_value();
        freq = json["freq"].number_value();
        varp = json["varp"].number_value();

        printf ("pulse data:\n");
        printf ("ampl : %12.4f\n", ampl);
        printf ("sigm : %12.4f\n", sigm);
        printf ("freq : %12.4f\n", freq);
        printf ("varp : %12.4f\n", varp);
        printf ("\n\n");
    }

    pulse (double _ampl, double _sigm, double _freq, double _varp):
           on(false), ampl(_ampl), sigm(_sigm), freq(_freq), varp(_varp) {}

    pulse (const pulse& rhs): on (false), ampl(rhs.ampl), sigm(rhs.sigm), freq(rhs.freq), varp(rhs.varp) {}

   ~pulse () {on=false, ampl=sigm=freq=varp=0;}

    pulse& operator=(const pulse& rhs) {
        if (this != &rhs) {
            on   = rhs.on;
            ampl = rhs.ampl;
            sigm = rhs.sigm;
            freq = rhs.freq;
            varp = rhs.varp;
        }
        return *this;
    }

    double et (double t) const {
        double _et = 0;
        if (on) {
            _et = (ampl*deom_pi/(sqrt(2.0*deom_pi)*sigm))*exp(-0.5*t*t/(sigm*sigm))*cos(freq*t+varp*t*t);
        } 
        return _et;
    }

    void turnoff () {on = false;}

    void turnon  () {on = true;}
};

#endif
