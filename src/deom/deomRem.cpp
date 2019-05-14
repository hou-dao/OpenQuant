/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include <cmath>
#include <cstdio>
#include <memory>
#include "deom.hpp"

void deom::remSch (cx_cube& dtotal, const cx_cube& total, const double t) {
    const int nsav = nddo;
    dtotal.head_slices(nddo).zeros();

    cx_cube qddo(size(qmd1));
    cx_cube ddoq(size(qmd1));

    for (int iado=0; iado<nsav; ++iado) {
        const cx_mat& ado = total.slice(iado);
        if (iado==0 || is_valid (ado)) {
            shared_ptr<TrieNode> nod = keys[iado];

            dtotal.slice(iado) += -deom_ci*(ham1*ado-ado*ham1)-nod->gams*ado;
            for (int m=0; m<nmod; ++m) {
                qddo.slice(m) = qmd1.slice(m)*ado;
                ddoq.slice(m) = ado*qmd1.slice(m);
                if (abs(delt_res(m)) > 1.e-15) {
                    dtotal.slice(iado) -= delt_res(m)*(
                        qmd1.slice(m)*(qddo.slice(m)-ddoq.slice(m))
                       -(qddo.slice(m)-ddoq.slice(m))*qmd1.slice(m));
                }
            }

            if (nod->tier < lmax) {
                for (int mp=0; mp<nind; ++mp) {
                    ivec key1 = gen_key(nod->nvec, mp, 1);
                    const int m = modLabel(mp);
                    const int n = key1(mp);
                    const cx_double sn = -deom_ci*sqrt(n/coef_abs(mp));
                    const cx_double cl = sn*coef_lft(mp);
                    const cx_double cr = sn*coef_rht(mp);
                    shared_ptr<TrieNode> ptr = tree.try_insert(key1,expn,nddo);
                    const int jddo = ptr->rank;
                    if (jddo != nddo) {
                        dtotal.slice(jddo)+= cl*qddo.slice(m)-cr*ddoq.slice(m);
                        if (jddo < nsav) {
                            const cx_double sm = -deom_ci*sqrt(n*coef_abs(m));
                            dtotal.slice(iado)+= sm*(qmd1.slice(m)*total.slice(jddo)-total.slice(jddo)*qmd1.slice(m)); 
                        }
                    } else {
                        dtotal.slice(nddo) = cl*qddo.slice(m)-cr*ddoq.slice(m);
                        keys[nddo] = ptr;
                        nddo += 1;
                    }
                }
            }
        }
    }
}
