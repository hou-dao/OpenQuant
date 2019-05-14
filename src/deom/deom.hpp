/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#ifndef DEOM_H_
#define DEOM_H_
   
#include "armadillo"
#include "Trie.hpp"
#include "JsonParser.hpp"

#include "deomConst.hpp"
#include "deomSyst.hpp"
#include "deomBath.hpp"
#include "deomHidx.hpp"
#include "deomPulse.hpp"

using namespace std;
using namespace arma;
using namespace json11;

class deom: public syst, public bath, public hidx {

    public:

    cx_cube ddos1;
    cx_cube ddos2;

    deom (const Json& json): syst (json["syst"]), bath (json["bath"]), hidx (json["hidx"]) {
        ddos1.set_size(nsys,nsys,nmax);
        ddos2.set_size(nsys,nsys,nmax);
    }

    deom (const syst& s, const bath& b, const hidx& h): syst (s), bath (b), hidx (h) {
        ddos1.set_size(nsys,nsys,nmax);
        ddos2.set_size(nsys,nsys,nmax);
    }

    deom (const deom& rhs): syst(rhs.ham1,rhs.qmd1), 
        bath(rhs.temperature, rhs.modLabel, rhs.coef_lft, rhs.coef_rht, rhs.coef_abs, rhs.expn_gam, rhs.delt_res), 
        hidx(rhs.nind, rhs.lmax, rhs.nmax, rhs.lddo, rhs.nddo, rhs.ferr, rhs.tree, rhs.keys, rhs.expn) {
        ddos1.set_size(nsys,nsys,nmax);
        ddos2.set_size(nsys,nsys,nmax);
    }

    ~deom () {}


    void oprAct (cx_cube& d_ddos, const cx_mat& sdip, const cx_cube& ddos, const char lrc='l');

    void iniHei (cx_cube& d_ddos, const cx_mat& sdip);

    void remSch (cx_cube& d_ddos, const cx_cube& ddos, const double t);

    void rem (cx_cube& d_ddos, const cx_cube& ddos, const double t, const char sch_hei = 's') {
        if (sch_hei == 's') {
            remSch (d_ddos, ddos, t);
        } else {
            printf("sch_hei is invalid!\n");
        }
    }

    void propagation (cx_cube& ddos, const double dt=0.005, const int nt=1000, const int nk=10);

    void equilibrium (cx_cube& ddos, const double dt=0.005, const double err=2.0e-8, const int nk=10, const string& method=string("Scit") ) {
        if (method == "Scit") {
            EqSolverScit (ddos, dt, err, nk);
        } else if (method == "Prop") {
            EqSolverProp (ddos, dt, err, nk);
        } else {
            printf ("Wrong method!\n");
        }
    }

    void EqSolverScit (cx_cube& ddos, const double dt=0.005, const double err=2.0e-8, const int nk=10);

    void EqSolverProp (cx_cube& ddos, const double dt=0.005, const double err=2.0e-8, const int nk=10);

    inline bool is_valid (const cx_mat& ddo) const {
        // return any(abs(vectorise(ddo))>ferr);
        for (unsigned int i=0; i<ddo.n_rows; ++i) {
            for (unsigned int j=0; j<ddo.n_cols; ++j) {
                if (abs(ddo(j,i))>ferr) return true;
            }
        }
        return false;
    }

    void filter (cx_cube& ddos) {
        int n = 1;
        int l = 0;
        for (int iddo=1; iddo<nddo; ++iddo) {
            shared_ptr<TrieNode> p = keys[iddo];
            if (is_valid(ddos.slice(iddo))) {
                if (n != iddo) {
                    p->rank = n;
                    keys[n] = p;
                    ddos.slice(n) = ddos.slice(iddo);
                }
                l = l>(p->tier)?l:(p->tier);
                ++n;
            } else {
                p->rank = -9527;
                keys[iddo] = nullptr;
            }
        }
        lddo = l;
        nddo = n;
    }

    bool notconverged (const int& nddo_backup, const hkey& keys_backup, 
                    const cx_cube& ddos_backup, const cx_cube& ddos, const double& tol) {
        for (int iddo=0; iddo<nddo_backup; ++iddo) {
            shared_ptr<TrieNode> nod = keys_backup[iddo];
            double maxDiff = 0.0;
            if (nod != nullptr && nod->rank>=0)
                maxDiff = max(abs(vectorise(ddos_backup.slice(iddo)-ddos.slice(nod->rank))));
            else
                maxDiff = max(abs(vectorise(ddos_backup.slice(iddo))));
            if (maxDiff > tol) {
                printf("Error bigger than %16.6e, while tol=%16.6e\n", maxDiff, tol);
                return true;
            }
        }
        return false;
    }

    template<typename... Tc>
    void rk4 (cx_cube& ddos, const double t, const double dt, const Tc&... args) {
        cx_cube* rhot = &ddos;
        for (int k=0; k<4; ++k) {
            const double dtp = dt/(4-k);
            const int ntmp = nddo;
            rem (ddos1,*rhot,t,args...);
            rhot = (k<3)?(&ddos2):(&ddos);
            rhot->head_slices(ntmp) = ddos.head_slices(ntmp)+ddos1.head_slices(ntmp)*dtp;
            if (nddo > ntmp)
                rhot->slices(ntmp,nddo-1) = ddos1.slices(ntmp,nddo-1)*dtp;
        }
        filter (ddos);
    }

    cx_double Trace (const cx_mat& sdip, const cx_cube& ddos) const {
        return trace(sdip*ddos.slice(0));
    }

    double entropy (const cx_cube& ddos, const string& type=string("vn")) const {
        vec popt = zeros<vec>(nsys);
        if (type == "vn") {
            eig_sym(popt,ddos.slice(0));
        } else if (type == "sh") {
            vec eval = zeros<vec>(nsys);
            cx_mat evec = zeros<cx_mat>(nsys,nsys);
            eig_sym(eval, evec, ham1);
            popt = real(diagvec(evec.st()*ddos.slice(0)*evec));
        }
        double s = 0.0;
        for (int i=0; i<nsys; ++i) {
            if (popt(i) < 0)
                return -9527;
            else if (popt(i) > 1.e-20)
                s -= popt(i)*log(popt(i));
        }
        return s;
    }

    double Umicro (const cx_cube& ddos) const {
        return real(trace(ham1*ddos.slice(0)));
    }
};

void copyright();

#endif
