/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include "deom.hpp"

int main (int argc, char *argv[])
{
    Json json;

    if (JsonParser(argc, argv, json))
    {
        printf("Error in reading input file!\n");
    }
    else
    {
        copyright();

        deom d(json["deom"]);

        const int    nt = json["nt"].int_value();
        const int    nk = json["nk"].int_value();
        const double dt = json["dt"].number_value();
        const int inistate = json["inistate"].int_value();

        cx_cube ddos = zeros<cx_cube>(d.nsys,d.nsys,d.nmax);
        ddos(inistate,inistate,0) = 1.0;

        FILE *frho = fopen("prop-rho.dat","w");
        FILE *fent = fopen("entropy.dat","w");
        for (int it=0; it<nt; ++it) {
            const double t = it*dt;
            if (it%nk == 0) {
                printf ("Propagation %5.1f%%: nddo=%6d, lddo=%3d\n", 100*it/static_cast<double>(nt), d.nddo, d.lddo);
                fprintf (frho, "%16.6e", t);
                for (int i=0; i<d.nsys; ++i)
                    fprintf (frho, "%20.10e", real(ddos(i,i,0)));
                fprintf (frho, "\n");
                double s1 = d.entropy(ddos,"vn");
                double s2 = d.entropy(ddos,"sh");
                fprintf (fent, "%16.6e%16.6e%16.6e\n", t, s1, s2);
            }
            d.rk4 (ddos,t,dt);
        }
        fclose (fent);
        fclose (frho);
    }

    return 0;
}
