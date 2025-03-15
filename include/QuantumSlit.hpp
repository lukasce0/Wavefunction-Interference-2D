#ifndef __QuantumSlit_hpp__
#define __QuantumSlit_hpp__

#include <iostream>
#include <armadillo>
#include <vector>
#include <string>
#include <cmath>
#include <math.h>
#include <complex>





class QuantumSlit {

public:
	int dim; //dim = M-2, meaning the number of points in the room along one axis
	int N_elements; //(M-2)^2
	arma::sp_cx_dmat A; //the first metrix for time evolution
	arma::sp_cx_dmat B; //second metrix for time evolution
	arma::cx_vec a; //diagonal elements for metrix A
	arma::cx_vec b; //diagonal elements for metrix B
	double dx; //the in room, meaning dx = h
	double dt; //step in time
	arma::cx_double r; //values outside of the diagonal of A and B
	std::vector<arma::cx_vec> results; //vector containing solution at differant times
	arma::vec V; //potential
	double v0; //vale of potential in the wall



	QuantumSlit(double step_room, double step_time, int slits, double wall_potential);


	void initialU(double xc, double yc, double sx, double sy, double px, double py); //initial wave function


	void solve(double T); //evaluate the wavefunction for time T


private: //these function are meant to only be accessed from within the object



	void constructMatrix(); //finds vectors a and b and assignes values to metrix A and B


	void findPotentialSingle(); //finds potential for single slit


	void findPotentialDouble(); //finds potential for double slit


	void findPotentialTriple(); //finds potential for triple slit


	void findVectors(); //find vector a and b


	arma::cx_vec makeVectorComp(arma::sp_cx_dmat A); //transforms a complex (M-2)x(M-2) metrix into a vector with (M-1)^2 elements


	arma::vec makeVectorRel(arma::mat A); //transforms a real (M-2)x(M-2) metrix into a vector with (M-1)^2 elements

};




#endif