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

static void corr (const int& nt, const double& dt, const double& staticErr, const int& nk,
                  const cx_mat& sdip1, const cx_mat& sdip2,
                  const syst& s, const bath& b, const hidx& h) {

    deom d0(s,b,h);

    cx_cube rho0 = zeros<cx_cube>(d0.nsys,d0.nsys,d0.nmax);
    cx_cube rho1 = zeros<cx_cube>(d0.nsys,d0.nsys,d0.nmax);
    cx_cube rho2 = zeros<cx_cube>(d0.nsys,d0.nsys,d0.nmax);

    d0.equilibrium (rho0,dt,staticErr,nk);
    //
    cx_double avg1 = d0.Trace(sdip1,rho0);
    cx_double avg2 = d0.Trace(sdip2,rho0);
    FILE *fs = fopen("stationary.dat","w");
    fprintf(fs, "%16.6e%16.6e%16.6e%16.6e\n", real(avg1), imag(avg1), real(avg2), imag(avg2));
    fclose(fs);

    printf("Equilibrium pop:");
    for (int i=0; i<d0.nsys; ++i)
        printf("%16.6e", real(rho0(i,i,0)));
    printf("\n");

    // c11, c12, c21, c22
    cx_vec ct11 = zeros<cx_vec>(nt);
    cx_vec ct12 = zeros<cx_vec>(nt);
    cx_vec ct21 = zeros<cx_vec>(nt);
    cx_vec ct22 = zeros<cx_vec>(nt);

    deom d1(d0);
    deom d2(d0);
    // \mu\rho
    d1.oprAct(rho1,sdip1,rho0,'l');
    d2.oprAct(rho2,sdip2,rho0,'l');
    // iTr[\mu G(t)\mu\rho]
    FILE *fc11 = fopen("corr-ll.t","w");
    FILE *fc12 = fopen("corr-lr.t","w");
    FILE *fc21 = fopen("corr-rl.t","w");
    FILE *fc22 = fopen("corr-rr.t","w");
    for (int it=0; it<nt; ++it) {
        double t = it*dt;
        ct11(it) = d1.Trace(sdip1,rho1)-avg1*avg1;
        ct12(it) = d1.Trace(sdip2,rho1)-avg1*avg2;
        ct21(it) = d2.Trace(sdip1,rho2)-avg2*avg1;
        ct22(it) = d2.Trace(sdip2,rho2)-avg2*avg2;
        if (it%nk == 0) {
            printf ("In corr: it=%d, nddo=%d, lddo=%d\n", it, d1.nddo, d1.lddo);
            printf ("In corr: it=%d, nddo=%d, lddo=%d\n", it, d2.nddo, d2.lddo);
        }
        fprintf(fc11, "%16.6e%16.6e%16.6e\n", t/deom_fs2unit, real(ct11(it)), imag(ct11(it)));
        fprintf(fc12, "%16.6e%16.6e%16.6e\n", t/deom_fs2unit, real(ct12(it)), imag(ct12(it)));
        fprintf(fc21, "%16.6e%16.6e%16.6e\n", t/deom_fs2unit, real(ct21(it)), imag(ct21(it)));
        fprintf(fc22, "%16.6e%16.6e%16.6e\n", t/deom_fs2unit, real(ct22(it)), imag(ct22(it)));
        d1.rk4 (rho1,t,dt);
        d2.rk4 (rho2,t,dt);
    }
    fclose(fc22);
    fclose(fc21);
    fclose(fc12);
    fclose(fc11);

    // 1D FFT
    ct11(0) *= 0.5;
    ct12(0) *= 0.5;
    ct21(0) *= 0.5;
    ct22(0) *= 0.5;
    const double dw = 2.0*deom_pi/(nt*dt);
    const cx_vec& cw11 = ifft(ct11)*nt*dt;
    const cx_vec& cw12 = ifft(ct12)*nt*dt;
    const cx_vec& cw21 = ifft(ct21)*nt*dt;
    const cx_vec& cw22 = ifft(ct22)*nt*dt;

    fc11 = fopen("corr-ll.w","w");
    fc12 = fopen("corr-lr.w","w");
    fc21 = fopen("corr-rl.w","w");
    fc22 = fopen("corr-rr.w","w");
    for (int iw=nt/2; iw<nt; ++iw) {
        double w = (iw-nt)*dw/deom_cm2unit;
        fprintf(fc11, "%16.6e%16.6e%16.6e\n", w, real(cw11(iw)), imag(cw11(iw)));
        fprintf(fc12, "%16.6e%16.6e%16.6e\n", w, real(cw12(iw)), imag(cw12(iw)));
        fprintf(fc21, "%16.6e%16.6e%16.6e\n", w, real(cw21(iw)), imag(cw21(iw)));
        fprintf(fc22, "%16.6e%16.6e%16.6e\n", w, real(cw22(iw)), imag(cw22(iw)));
    }
    for (int iw=0; iw<nt/2; ++iw) {
        double w = iw*dw/deom_cm2unit;
        fprintf(fc11, "%16.6e%16.6e%16.6e\n", w, real(cw11(iw)), imag(cw11(iw)));
        fprintf(fc12, "%16.6e%16.6e%16.6e\n", w, real(cw12(iw)), imag(cw12(iw)));
        fprintf(fc21, "%16.6e%16.6e%16.6e\n", w, real(cw21(iw)), imag(cw21(iw)));
        fprintf(fc22, "%16.6e%16.6e%16.6e\n", w, real(cw22(iw)), imag(cw22(iw)));
    }
    fclose(fc22);
    fclose(fc21);
    fclose(fc12);
    fclose(fc11);
}

int main (int argc, char *argv[]) {

    Json json;

    if (JsonParser(argc, argv, json)) {
        printf("Error in reading input file!\n");
    } else {

        copyright();

        syst s(json["deom"]["syst"]);
        bath b(json["deom"]["bath"]);
        hidx h(json["deom"]["hidx"]);

        const int    nt = json["nt"].int_value();
        const int    nk = json["nk"].int_value();
        const double dt = json["dt"].number_value();
        const double staticErr = json["staticErr"].number_value();
        cx_mat sdip1 = json2zmat(json["sdip1"]);
        cx_mat sdip2 = json2zmat(json["sdip2"]);

        corr (nt, dt, staticErr, nk, sdip1, sdip2, s, b, h);
    }

    return 0;
}
