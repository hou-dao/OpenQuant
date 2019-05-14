/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#ifndef DEOMHIDX_H_
#define DEOMHIDX_H_

#include <map>
#include <memory>
#include <cstdio>
#include <string>
#include <vector>
#include "armadillo"
#include "JsonParser.hpp"
#include "Trie.hpp"

using namespace std;
using namespace arma;
using namespace json11;

typedef vector<shared_ptr<TrieNode> > hkey;

class hidx 
{

public:

    int    nind;
    int    lmax;
    int    nmax;
    int    lddo;
    int    nddo;
    double ferr;
    Trie   tree;
    hkey   keys;
    cx_vec expn;

    hidx (const Json& json) 
    {
        expn = json2zvec(json["expn"]);
        nind = expn.n_rows;
        lmax = json["lmax"].int_value();
        nmax = json["nmax"].int_value();
        ferr = json["ferr"].number_value();

        int nmax_kl = static_cast<int>(get_nmax(nind,lmax));
        if (nmax_kl<nmax) {
            nmax = nmax_kl;
        }

        keys.resize(nmax); 
        tree.init(lmax+1);
        keys[0] = tree.try_insert(zeros<ivec>(1),expn,0);
        nddo = 1;
        lddo = 0;

        printf ("$InitHidx\n");
        printf ("nind = %d\n", nind);
        printf ("lddo = %d\n", lddo);
        printf ("nddo = %d\n", nddo);
        printf ("lmax = %d\n", lmax);
        printf ("nmax = %d\n", nmax);
        printf ("ferr = %g\n", ferr);
        printf ("$InitHidx\n\n");
    }

    hidx (const hidx& rhs): nind(rhs.nind), lmax(rhs.lmax),nmax(rhs.nmax),
        lddo (rhs.lddo), nddo(rhs.nddo), ferr(rhs.ferr), 
        tree (rhs.tree), keys(rhs.keys), expn(rhs.expn) {}

    hidx (const int _nind, const int _lmax, const int _nmax, const int _lddo, const int _nddo,
          const double _ferr, const Trie& _tree, const hkey& _keys, const cx_vec& _expn):
          nind(_nind), lmax(_lmax), nmax(_nmax), lddo(_lddo), nddo(_nddo), ferr(_ferr), 
          tree(_tree), keys(_keys), expn(_expn) {}

   ~hidx () {}

    hidx& operator= (const hidx& rhs) {
        if (this != &rhs) {
            nind = rhs.nind;
            lmax = rhs.lmax;
            nmax = rhs.nmax;
            lddo = rhs.lddo;
            nddo = rhs.nddo;
            ferr = rhs.ferr;
            tree = rhs.tree;
            keys = rhs.keys;
            expn = rhs.expn;
        }
        return *this;
    }
    
    /* The maximal number of ddos given K and L
     *
     */
    unsigned long get_nmax (const int K, const int L) const {
        unsigned long ntot = 1;
        for (int k=1; k<=K; ++k) {
            ntot *= L+k;
            ntot /= k;
            if (ntot > 300000) {
                // printf ("Be careful! too many elements!\n");
                return nmax;
            }
        }
        return ntot;
    }

    inline ivec gen_key (const ivec& str, const int pos, const int chg) const {
        int len0 = str.n_rows;
        if (len0>=pos+1 && chg>0) {
            ivec key (str);
            key(pos) += chg;
            return key;
        } else if (len0<pos+1 && chg>0) {
            ivec key = zeros<ivec>(pos+1);
            for (int i=0; i<len0; ++i)
                key(i) = str(i);
            key(pos) += chg;
            return key;
        } else if (len0>=pos+1 && chg<0 && str(pos)+chg>=0) {
            ivec key(str);
            key(pos) += chg;
            int npos = len0;
            while (key(npos-1) == 0 && npos!=1)
                npos -= 1;
            return key.head(npos);
        }
        return ivec();
    }

    void init_hierarchy () {
        int ni = 0;
        int nf = nddo;
        for (int l=0; l<lmax; ++l) {
            for (int iddo=ni; iddo<nf; ++iddo) {
                shared_ptr<TrieNode> nod = keys[iddo];
                for (int mp=0; mp<nind; ++mp) {
                    ivec key1 = gen_key(nod->nvec, mp, 1);
                    shared_ptr<TrieNode> ptr = tree.try_insert(key1,expn,nddo);
                    if (ptr->rank == nddo) {
                        keys[nddo] = ptr;
                        nddo += 1;
                    }
                }
            }
            ni = nf;
            nf = nddo;
        }
        lddo = lmax;
        nmax = nddo;
    }

    ivec gen_nbar (ivec nvec) {
        int k = 0;
        int len = nvec.n_rows;
        ivec nbar = zeros<ivec>(len+1);
        while (k<len) {
            if (abs(imag(expn(k))) < 1.0e-20) {
                nbar(k) = nvec(k);
                k += 1;
            } else {
                if (k+1<len)
                    nbar(k) = nvec(k+1);
                else
                    nbar(k) = 0;
                nbar(k+1) = nvec(k);
                k += 2;
            }
        }
        return (nbar(len)==0)?nbar.head_rows(len):nbar;
    }

    void status() {
        ivec tcount = zeros<ivec>(lddo+1);
        ivec kcount = zeros<ivec>(nind);
        for (int i=0; i<nddo; ++i) {
            int l = keys[i]->tier;
            tcount(l) += 1;
            for (unsigned int j=0; j<(keys[i]->nvec).n_rows; ++j) {
                kcount(j) = kcount(j)>(keys[i]->nvec(j))?kcount(j):(keys[i]->nvec(j));
            }
        }
        printf("DDO index histogram\n");
        for (int l=0; l<=lddo; ++l) {
            int pcen = (int)((double)tcount(l)/(double)nddo*200);
            printf("Tier %3d: ", l);
            for (int i=0; i<pcen; ++i) {
                printf("*");
            }
            printf("(%d)\n",(int)tcount(l));
        }
        kcount.print("Max Levels");
    }
};

#endif
