/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include <cstdio>
#include <cstdlib>
#include <memory>
#include "deom.hpp"

static const int ntMax = 100000000;

void deom::EqSolverScit (cx_cube& ddos, const double dt, const double err, const int nk) {

    printf("$SCI session!\n");

    vec eigh = eig_sym(ham1);
    const double OMG = 2*max(max(abs(eigh)),sqrt(lmax*max(coef_abs)));
    printf ("OMG=%16.6e\n",OMG);

    if (ddos.is_empty() || abs(trace(ddos.slice(0))-1.0)>1.e-15) 
    {
        ddos = zeros<cx_cube>(nsys,nsys,nmax);
        mat eham = expmat_sym(-real(ham1)/mean(temperature));
        ddos.slice(0) = eham/trace(eham)*deom_c1;
    } 

    cx_mat rho_old(nsys,nsys);
    cx_mat htmp(deom_ci*ham1);
    cx_mat hlft(htmp);
    cx_mat hrht(htmp);
    double max_diff = 0.0;
    int iter = 0;
    do 
    { 
        max_diff = 0.0;
        for (int iddo=0; iddo<nddo; ++iddo) 
        {
            shared_ptr<TrieNode> nod = keys[iddo];
            const bool flag = is_valid(ddos.slice(iddo));

            cx_mat rhob = zeros<cx_mat>(nsys,nsys);
            for (int m=0; m<nmod; ++m) 
            {
                rho_old = qmd1.slice(m)*ddos.slice(iddo)-ddos.slice(iddo)*qmd1.slice(m);
                if (abs(delt_res(m))>1.0e-15) 
                {
                    rhob -= delt_res(m)*(qmd1.slice(m)*rho_old-rho_old*qmd1.slice(m));
                }
            }
            if (nod->tier < lmax) 
            {
                for (int mp=0; mp<nind; ++mp) 
                {
                    ivec key1 = gen_key(nod->nvec, mp, 1);
                    const int m = modLabel(mp);
                    const int n = key1(mp);
                    const cx_double sn = -deom_ci*sqrt(n*coef_abs(m));
                    shared_ptr<TrieNode> ptr = tree.try_insert(key1,expn,-1);
                    const int jddo = ptr->rank;
                    if (jddo>=0 && jddo!=nddo)
                    {
                        rhob += sn*(qmd1.slice(m)*ddos.slice(jddo)-ddos.slice(jddo)*qmd1.slice(m));
                    } 
                    else if (flag) 
                    {
                        ddos.slice(nddo).zeros();
                        ptr->rank = nddo;
                        keys[nddo] = ptr;
                        nddo += 1;
                    }
                }
            }

            for (int mp=0; mp<nind; ++mp) 
            {
                ivec key1 = gen_key(nod->nvec, mp, -1);
                if (!key1.is_empty()) 
                {
                    const int m = modLabel(mp);
                    const int n = nod->nvec(mp);
                    const cx_double sn = -deom_ci*sqrt(n/coef_abs(mp));
                    const cx_double cl = sn*coef_lft(mp);
                    const cx_double cr = sn*coef_rht(mp);
                    shared_ptr<TrieNode> ptr = tree.try_insert(key1,expn,-1);
                    const int jddo = ptr->rank;
                    if (jddo>=0 && jddo!=nddo)
                    {
                        rhob += cl*qmd1.slice(m)*ddos.slice(jddo)-cr*ddos.slice(jddo)*qmd1.slice(m); 
                    } 
                    else if (flag) 
                    {
                        ddos.slice(nddo).zeros();
                        ptr->rank = nddo;
                        keys[nddo] = ptr;
                        nddo += 1;
                    }
                }
            }

            rhob += OMG*ddos.slice(iddo);

            rho_old = ddos.slice(iddo);
            hlft = htmp;
            hlft.diag() += OMG+nod->gams;
            hrht =-htmp;
            if (iddo != 0)
            {
                ddos.slice(iddo) = syl(hlft,hrht,-rhob);
            } 
            else 
            {
                cx_mat rhox = syl(hlft,hrht,-rhob);
                ddos.slice(0) = 0.5*(rhox+rhox.t());
            }
            double diff = max(vectorise(abs(rho_old-ddos.slice(iddo))));
            max_diff = max_diff>diff?max_diff:diff;
        }
        
        iter += 1;

        filter(ddos);

        printf("SCI step %d, nddo=%d, lddo=%d, max_diff=%16.6e\n", iter, nddo, lddo, max_diff);
        vec pop = real(ddos.slice(0).diag());
        pop.print("population");

        if (max_diff < err) break;

    } while (iter < ntMax);

    ddos.slice(0).save("rho0.real",raw_ascii);
    printf("$SCI session!\n");
}
