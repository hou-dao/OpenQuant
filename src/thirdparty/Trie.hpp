#ifndef TRIE_H_
#define TRIE_H_

#include <cstdio>
#include <memory>
#include <vector>
#include "armadillo"

class TrieNode {

    public:

    TrieNode (): gams(0.0), tier(0), rank(-9527) {}

    TrieNode (const int n_child): gams(0.0), tier(0), rank(-9527), child(n_child) {}

    inline TrieNode (const int n_child, const arma::cx_double _gams, const int _tier, const int _rank, const arma::ivec& _nvec): 
        gams(_gams), tier(_tier), rank(_rank), nvec(_nvec), child(n_child) {}

   ~TrieNode () {}

    arma::cx_double   gams;
    int               tier;
    int               rank;
    arma::ivec        nvec;
    std::vector<std::shared_ptr<TrieNode> > child;
};

class Trie {

    public:

    Trie (): n_child(0), n_size(1) {
        root = std::make_shared<TrieNode>();
    }

    void init (const int _n_child) {
        if (n_child==0 && n_size==1 && root->child.empty()) {
            n_child = _n_child;
            root->child = std::vector<std::shared_ptr<TrieNode> >(n_child);
        } else {
            printf("Error in initializing Trie!\n");
        }
    }

    Trie (const int _n_child): n_child(_n_child), n_size(1) {
        root = std::make_shared<TrieNode>(n_child);
    }

    Trie (const Trie& rhs): n_child(rhs.n_child), n_size(rhs.n_size) {
        root = rhs.root;
    }

   ~Trie () {}

    // Try to insert a key-value pair; return ptr if insertion succeed or node exist
    std::shared_ptr<TrieNode> try_insert (const arma::ivec& nvec, const arma::cx_vec& expn, const int rank) {
        std::shared_ptr<TrieNode> p(root);
        for (unsigned int k=0; k<nvec.n_rows; ++k) {
            const int nk = nvec[k];
            if (p->child[nk]) {
                p = p->child[nk];
            } else {
                p = p->child[nk] = std::make_shared<TrieNode>(n_child-p->tier-nk,p->gams+expn(k)*static_cast<double>(nk),p->tier+nk,-9527,nvec.head(k+1));
                ++n_size;
            }
        }
        if (p->rank<0)
            p->rank = rank;
        return p;
    }

    // Return node pointer if key is in the trie.
    inline std::shared_ptr<TrieNode> find (const arma::ivec& nvec) const {
        std::shared_ptr<TrieNode> p(root);
        for (unsigned int k=0; k<nvec.n_rows; ++k) {
            const int nk = nvec[k];
            if (!p->child[nk])
                return nullptr;
            else
                p = p->child[nk];
        }
        return p;
    }

    int size () const {
        return n_size;
    }

private:
    int n_child;
    int n_size;
    std::shared_ptr<TrieNode> root;
};

#endif
