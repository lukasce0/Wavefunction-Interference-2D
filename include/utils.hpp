#ifndef __utils_hpp__  
#define __utils_hpp__

#include <iostream>
#include <armadillo>
#include <vector>
#include <string>
#include <cmath>
#include <math.h>
#include <complex>






// A function that prints the structure of a sparse matrix to screen.
void print_sp_matrix_structure(const arma::sp_cx_mat& A);


void writeToFileReal(std::string file_name, QuantumSlit instance, int save_every); //writes real part of wavefunction into a file, saves every save_every element


void writeToFileImag(std::string file_name, QuantumSlit instance, int save_every); //writes imaginary part of wavefunction into a file, saves every save_every element


void writeToFile(std::string file_name, QuantumSlit instance, int save_every); //writes probability mass function into a file, saves every save_every element




#endif