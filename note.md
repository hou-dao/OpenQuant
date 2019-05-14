------------------------------------------
Dissipaton Equation Of Motion (DEOM)
------------------------------------------

1. Structure
   deom/
     |
     |-input/ 
     |
     |-src/
     |   |
     |   |-apps/ (single-file projects)
     |   |   |
     |   |   |-rhot.cpp (quantum dynamics)
     |   |   |
     |   |   |-corr.cpp (correlation)
     |   |
     |   |--deom/ (deom core library)
     |   |
     |   |--thirdparty/ (trie tree; handle hierarchy index)
     |
     |--scripts/ (scripts for generating input files)
     |
     |--CMakeLists.txt 
     |
     |--README.md

2. Dependencies

   json11 (for parsing input)
   armadillo (linear algebra)
   blas & lapack (linear algebra)

   If there's mkl, append "-mkl" to CMAKE_CXX__FLAGS and delete "blas lapack" in the target_link_libraries in the CMakeLists;
   If there's no mkl, but there's openblas, use "openblas" instead of "blas".

3. Compile

   Requirement: c++ compiler that support c++11 (g++-4.7 or above); 
                cmake; 
                git (fetch json11 from github.com);

   mkdir deom && cd deom && cp /path_to_sourcecodes/deom.tar.gz .
   tar xvf deom.tar.gz
   mkdir build && cd build
   cmake .. && make

The binary files will be found in deom/bin

4. Usage

   step 1. cd deom/input
   step 2. Edit gen_input.py as you need
   step 3. python gen_input.py
   step 4. ../bin/{app}

