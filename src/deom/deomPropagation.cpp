/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include <cstdio>
#include <cstdlib>
#include "deom.hpp"

void deom::propagation (cx_cube& ddos, const double dt, const int nt, const int nk) {

    FILE *frho = fopen("propagation.log","w");
    for (int it=0; it<nt; ++it) {
        const double t = it*dt;
        printf ("Propagation %5.1f%%: nddo=%6d, lddo=%3d\n",
                100*it/static_cast<double>(nt), nddo, lddo);
        if (it%nk == 0) {
            fprintf (frho, "%16.6e", t/deom_fs2unit);
            for (int i=0; i<nsys; ++i) {
                fprintf (frho, "%16.6e", real(ddos(i,i,0)));
            }
            fprintf (frho, "\n");
        }
        rk4 (ddos,t,dt);
    }
    fclose (frho);
}
