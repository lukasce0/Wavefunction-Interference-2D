#include <iostream>
#include <armadillo>
#include <vector>
#include <string>
#include <cmath>
#include <math.h>
#include <complex>


#include "QuantumSlit.hpp"
#include "utils.hpp"




int main() {

	//
	//Solves the SE without a wall in order to study the normalization over time
	//

	double dx = 0.005; //room step (equivalent to h)
	double dt = 2.5e-5; //step in time
	int slits = 2; //number of slits
	double wall_potential = 0.; //potential in the wall (always zero outside of it)

	QuantumSlit no_wall_prob(dx, dt, slits, wall_potential); //create an instance without barriers
	

	double xc = 0.25; //center of initial wave function along x
	double yc = 0.5; //center of initial wave function along y
	double sx = 0.05; //standard diviation along x
	double sy = 0.05; //standard diviation along y
	double px = 200.; //momentum along x
	double py = 0.; //momentum along y
	no_wall_prob.initialU(xc, yc, sx, sy, px, py); //finds initial wavefunction

	double T = 0.008; //time of the simulation
	no_wall_prob.solve(T); //finds wavefunction for next time T

	writeToFile("results_probability_no_wall.txt", no_wall_prob, 10); //writes every 10th probability function to file
	//print_sp_matrix_structure(no_wall_prob.B); //find the structure of the metrix B


	//
	//Solves the SE with a wall (two slits) in order to study the normalization over time
	//

	sy = 0.05;
	wall_potential = 1e10;

	QuantumSlit wall_prob(dx, dt, slits, wall_potential); //instance with a double slit
	wall_prob.initialU(xc, yc, sx, sy, px, py);
	wall_prob.solve(T);
	writeToFile("results_probability_wall.txt", wall_prob, 10);



	//
	//Solves the SE with a wall (two slits) in order to plot the probability mass function and both imaginary and real parts of wavefunction
	//

	sy = 0.2;
	T = 0.002;

	QuantumSlit contur(dx, dt, slits, wall_potential); //creates double slit instance
	contur.initialU(xc, yc, sx, sy, px, py);
	contur.solve(T);
	writeToFile("results_probability_contur.txt", contur, 2);
	writeToFileReal("results_probability_contur_real.txt", contur, 2);
	writeToFileImag("results_probability_contur_imag.txt", contur, 2);


	//
	//Solves the SE with a wall (single slit) in order to plot the probability mass function
	//

	slits = 1;
	QuantumSlit contur_single(dx, dt, slits, wall_potential);
	contur_single.initialU(xc, yc, sx, sy, px, py);
	contur_single.solve(T);
	writeToFile("results_probability_contur_single.txt", contur_single, 2);



	//
	//Solves the SE with a wall (triple slit) in order to plot the probability mass function
	//

	slits = 3;
	QuantumSlit contur_triple(dx, dt, slits, wall_potential);
	contur_triple.initialU(xc, yc, sx, sy, px, py);
	contur_triple.solve(T);
	writeToFile("results_probability_contur_triple.txt", contur_triple, 2);



	
	return 0;
}