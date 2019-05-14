#ifndef JSONPARSER_H_
#define JSONPARSER_H_

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include "json11.hpp"
#include "armadillo"

inline int JsonParser (int argc, char *argv[], json11::Json& data)
{
    if (argc != 2 || (argc == 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")))
    {
        printf("Usage: path-to-binary InputFile\n");
        return -1;
    } 
    else if (argc == 2)
    {
        std::ifstream jsonFile(argv[1]);
        if (jsonFile)
        {
            std::stringstream strStream;
            strStream << jsonFile.rdbuf();
            std::string jsonStr = strStream.str();
            std::string err;
            data = json11::Json::parse(jsonStr,err);
            if (!err.empty())
            {
                printf ("Error in parsing input file: %s\n", err.c_str());
                return -3;
            }
        }
        else
        {
            printf ("Error in find jsonFile!\n");
            return -2;
        }
    }

    return 0;
}

inline arma::vec json2vec (const json11::Json& js)
{
    const unsigned int nrow = js["nrow"].int_value();
    arma::vec arr(nrow);
    if (nrow != js["real"].array_items().size())
        printf("Error in converting jsonarray to arma!\n");
    for (unsigned int i=0; i<nrow; ++i)
        arr(i) = js["real"].array_items()[i].number_value();
    return arr;
}

inline arma::ivec json2ivec (const json11::Json& js)
{
    const unsigned int nrow = js["nrow"].int_value();
    arma::ivec arr(nrow);
    if (nrow != js["real"].array_items().size())
        printf("Error in converting jsonarray to arma!\n");
    for (unsigned int i=0; i<nrow; ++i)
        arr(i) = js["real"].array_items()[i].int_value();
    return arr;
}

inline arma::cx_vec json2zvec (const json11::Json& js)
{
    const unsigned int nrow = js["nrow"].int_value();
    arma::cx_vec arr(nrow);
    if (nrow != js["real"].array_items().size() || nrow != js["imag"].array_items().size())
        printf("Error in converting jsonarray to arma!\n");
    for (unsigned int i=0; i<nrow; ++i)
        arr(i) = arma::cx_double(js["real"].array_items()[i].number_value(),js["imag"].array_items()[i].number_value());
    return arr;
}

inline arma::mat json2mat (const json11::Json& js)
{
    const unsigned int nrow = js["nrow"].int_value();
    const unsigned int ncol = js["ncol"].int_value();
    arma::mat arr(nrow,ncol);
    if (nrow*ncol != js["real"].array_items().size())
        printf("Error in converting jsonarray to arma!\n");
    for (unsigned int i=0; i<nrow; ++i)
        for (unsigned int j=0; j<ncol; ++j)
            arr(i,j) = js["real"].array_items()[i*ncol+j].number_value();
    return arr;
}

inline arma::imat json2imat (const json11::Json& js)
{
    const unsigned int nrow = js["nrow"].int_value();
    const unsigned int ncol = js["ncol"].int_value();
    arma::imat arr(nrow,ncol);
    if (nrow*ncol != js["real"].array_items().size())
        printf("Error in converting jsonarray to arma!\n");
    for (unsigned int i=0; i<nrow; ++i)
        for (unsigned int j=0; j<ncol; ++j)
            arr(i,j) = js["real"].array_items()[i*ncol+j].int_value();
    return arr;
}

inline arma::cx_mat json2zmat (const json11::Json& js)
{
    const unsigned int nrow = js["nrow"].int_value();
    const unsigned int ncol = js["ncol"].int_value();
    arma::cx_mat arr(nrow,ncol);
    if (nrow*ncol != js["real"].array_items().size() || nrow*ncol != js["imag"].array_items().size())
        printf("Error in converting jsonarray to arma!\n");
    for (unsigned int i=0; i<nrow; ++i)
        for (unsigned int j=0; j<ncol; ++j)
            arr(i,j) = arma::cx_double(js["real"].array_items()[i*ncol+j].number_value(),js["imag"].array_items()[i*ncol+j].number_value());
    return arr;
}

inline arma::cube json2cube (const json11::Json& js)
{
    const unsigned int nrow = js["nrow"].int_value();
    const unsigned int ncol = js["ncol"].int_value();
    const unsigned int nslc = js["nslc"].int_value();
    arma::cube arr(nrow,ncol,nslc);
    if (nrow*ncol*nslc != js["real"].array_items().size())
        printf("Error in converting jsonarray to arma!\n");
    for (unsigned int k=0; k<nslc; ++k)
        for (unsigned int i=0; i<nrow; ++i)
            for (unsigned int j=0; j<ncol; ++j)
                arr(i,j,k) = js["real"].array_items()[k*nrow*ncol+i*ncol+j].number_value();
    return arr;
}

inline arma::icube json2icube (const json11::Json& js)
{
    const unsigned int nrow = js["nrow"].int_value();
    const unsigned int ncol = js["ncol"].int_value();
    const unsigned int nslc = js["nslc"].int_value();
    arma::icube arr(nrow,ncol,nslc);
    if (nrow*ncol*nslc != js["real"].array_items().size())
        printf("Error in converting jsonarray to arma!\n");
    for (unsigned int k=0; k<nslc; ++k)
        for (unsigned int i=0; i<nrow; ++i)
            for (unsigned int j=0; j<ncol; ++j)
                arr(i,j,k) = js["real"].array_items()[k*nrow*ncol+i*ncol+j].int_value();
    return arr;
}

inline arma::cx_cube json2zcube (const json11::Json& js)
{
    const unsigned int nrow = js["nrow"].int_value();
    const unsigned int ncol = js["ncol"].int_value();
    const unsigned int nslc = js["nslc"].int_value();
    arma::cx_cube arr(nrow,ncol,nslc);
    if (nrow*ncol*nslc != js["real"].array_items().size() || nrow*ncol*nslc != js["imag"].array_items().size())
        printf("Error in converting jsonarray to arma!\n");
    for (unsigned int k=0; k<nslc; ++k)
        for (unsigned int i=0; i<nrow; ++i)
            for (unsigned int j=0; j<ncol; ++j)
                arr(i,j,k) = arma::cx_double(js["real"].array_items()[k*nrow*ncol+i*ncol+j].number_value(),
                                       js["imag"].array_items()[k*nrow*ncol+i*ncol+j].number_value());
    return arr;
}

#endif
