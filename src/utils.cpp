
#include <iostream>
#include <armadillo>
#include <vector>
#include <string>
#include <cmath>
#include <math.h>
#include <complex>


#include "QuantumSlit.hpp"






// A function that prints the structure of a sparse matrix to screen.
void print_sp_matrix_structure(const arma::sp_cx_mat& A)
{
	using namespace std;
	using namespace arma;

	// Declare a C-style 2D array of strings.
	string S[A.n_rows][A.n_cols];

	// Initialise all the strings to " ".
	for (int i = 0; i < A.n_rows; i++)
	{
		for (int j = 0; j < A.n_cols; j++)
		{
			S[i][j] = " ";
		}
	}

	// Next, we want to set the string to a dot at each non-zero element.
	// To do this we use the special loop iterator from the sp_cx_mat class
	// to help us loop over only the non-zero matrix elements.
	sp_cx_mat::const_iterator it = A.begin();
	sp_cx_mat::const_iterator it_end = A.end();

	int nnz = 0;
	for (it; it != it_end; ++it)
	{
		S[it.row()][it.col()] = ".";
		nnz++;
	}

	// Finally, print the matrix to screen.
	cout << endl;
	for (int i = 0; i < A.n_rows; i++)
	{
		cout << " | ";
		for (int j = 0; j < A.n_cols; j++)
		{
			cout << S[i][j] << " ";
		}
		cout << "|\n";
	}

	cout << endl;
	cout << "matrix size: " << A.n_rows << "x" << A.n_cols << endl;
	cout << "non-zero elements: " << nnz << endl;
	cout << endl;
}





void writeToFileReal(std::string file_name, QuantumSlit instance, int save_every) { //writes real part of wavefunction into a file, saves every save_every element

	std::string file = file_name; //filename
	std::ofstream ofile;
	ofile.open(file);

	for (int n = 0; n < instance.results.size(); n += save_every) { //sum over all times
		for (int k = 0; k < instance.N_elements; k++) { //sum over whole room
			ofile << instance.results.at(n).at(k).real() << " ";
		}
		ofile << std::endl;
	}

	ofile.close();
}




void writeToFileImag(std::string file_name, QuantumSlit instance, int save_every) { //writes imaginary part of wavefunction into a file, saves every save_every element

	std::string file = file_name; //filename
	std::ofstream ofile;
	ofile.open(file);

	for (int n = 0; n < instance.results.size(); n += save_every) { //sum over all times
		for (int k = 0; k < instance.N_elements; k++) { //sum over whole room
			ofile << instance.results.at(n).at(k).imag() << " ";
		}
		ofile << std::endl;
	}

	ofile.close();
}




void writeToFile(std::string file_name, QuantumSlit instance, int save_every) { //writes probability mass function into a file, saves every save_every element

	std::string file = file_name; //filename
	std::ofstream ofile;
	ofile.open(file);

	for (int n = 0; n < instance.results.size(); n += save_every) { //sum over all times
		for (int k = 0; k < instance.N_elements; k++) { //sum over whole room
			ofile << pow(abs(instance.results.at(n).at(k)), 2) << " ";
		}
		ofile << std::endl;
	}

	ofile.close();
}