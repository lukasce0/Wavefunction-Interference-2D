
#include <iostream>
#include <armadillo>
#include <vector>
#include <string>
#include <cmath>
#include <math.h>
#include <complex>


#include "QuantumSlit.hpp"









	QuantumSlit::QuantumSlit(double step_room, double step_time, int slits, double wall_potential) {
		dt = step_time;
		dx = step_room;
		dim = round((double)1. / dx - 1.);

		N_elements = pow(dim, 2);
		r = arma::cx_double(0, dt / (2 * dx * dx));
		v0 = wall_potential;


		//Assign values to potential metrix (V):
		if (slits == 1) {
			findPotentialSingle();
		}
		if (slits == 2) {
			findPotentialDouble();
		}
		if (slits == 3) {
			findPotentialTriple();
		}

		//Find vectors a and b, and afterwards metrix A and B
		constructMatrix();


	}





	void QuantumSlit::initialU(double xc, double yc, double sx, double sy, double px, double py) { //initial wave function

		arma::sp_cx_dmat u(dim, dim); //metrix of the wavefunction
		double x; //position in room along x
		double y; //position in room along y

		for (int i = 0; i < dim; i++) { //sum over points along x
			x = i * dx + dx;
			for (int j = 0; j < dim; j++) { //sum over points along y
				y = dx + j * dx;

				u(i, j) = arma::cx_double(exp(-pow(x - xc, 2) / (2 * pow(sx, 2)) - pow(y - yc, 2) / (2 * pow(sy, 2))), 0);
				u(i, j) *= arma::cx_double(cos(px * x), sin(px * x));
				u(i, j) *= arma::cx_double(cos(py * y), sin(py * y));
			}
		}


		double norm; //the normalization constant
		arma::cx_double u_loc; //value of wavefuncion at a point (formal requirement)


		for (int i = 0; i < dim; i++) { //sum over x
			for (int j = 0; j < dim; j++) { //sum over y
				u_loc = u(i, j);
				norm += real(u_loc * conj(u_loc));
			}
		}
		u /= sqrt(norm); //normalizing the initial wavefunction

		results.push_back(makeVectorComp(u)); //appending the initial vector to result list
	}




	void QuantumSlit::solve(double T) { //evaluate the wavefunction for time T

		arma::cx_vec c; //intermediate vector for solving a metrix equation
		arma::cx_vec u_loc; //next wavefunction (after dt)

		double N_time = T / dt; //number of time steps

		for (int n = 0; n <= N_time; n++) {
			c = B * results.at(n);
			u_loc = arma::spsolve(A, c);
			results.push_back(u_loc);
		}
	}








	void QuantumSlit::constructMatrix() { //finds vectors a and b and assignes values to metrix A and B

		A = arma::sp_cx_mat(N_elements, N_elements);
		B = arma::sp_cx_mat(N_elements, N_elements);

		findVectors(); //finds vectors a and b

		for (int i = 0; i < N_elements - dim; i++) //finds all sub matricies ((M-2)x(M-2)), except the one in lower right end
		{
			A(i, i) = a.at(i); //values along diagonal
			B(i, i) = b.at(i);

			A(i, i + dim) = -r; //values forthest away from diagonal
			A(i + dim, i) = -r;
			B(i, i + dim) = r;
			B(i + dim, i) = r;

			if ((i + 1) % dim != 0) { //values next to diagonal
				A(i, i + 1) = -r;
				A(i + 1, i) = -r;
				B(i, i + 1) = r;
				B(i + 1, i) = r;
			}
		}

		for (int i = N_elements - dim; i < N_elements; i++) {

			A(i, i) = a.at(i); //values along diagonal
			B(i, i) = b.at(i);

			if (i != dim - 1 and (i + 1) % dim != 0) { //values next to diagonal
				A(i, i + 1) = -r;
				A(i + 1, i) = -r;
				B(i, i + 1) = r;
				B(i + 1, i) = r;
			}
		}
	}






	void QuantumSlit::findPotentialSingle() { //finds potential for single slit

		double x_loc; //local position x
		double y_loc; //local position y
		bool condition_x;
		bool condition_y;
		arma::mat V_mat = arma::mat(dim, dim); //metrix of potential

		for (int j = 0; j < dim; j++) { //sum over y
			y_loc = dx + j * dx;

			for (int i = 0; i < dim; i++) { //sum over x

				x_loc = dx + i * dx;


				condition_x = (x_loc > 0.49) and (x_loc < 0.51); //if true x_loc is in the wall
				condition_y = not (y_loc > 0.475 and y_loc < 0.525); //if true y_loc is in the wall


				if (condition_x and condition_y == true) { //if (x_loc, y_loc) in wall
					V_mat(i, j) = v0;
				}
				else { //if (x_loc, y_loc) not in wall
					V_mat(i, j) = 0.;
				}
			}
		}
		V = makeVectorRel(V_mat); //transforms the result into a vector, such it can be used directly in finding the matricies A and B
	}






	void QuantumSlit::findPotentialDouble() { //finds potential for double slit

		double x_loc; //local position x
		double y_loc; //local position y
		bool condition_x;
		bool condition_y;
		arma::mat V_mat = arma::mat(dim, dim); //metrix of potential

		for (int j = 0; j < dim; j++) { //sum over y
			y_loc = dx + j * dx;

			for (int i = 0; i < dim; i++) { //sum over x

				x_loc = dx + i * dx;


				condition_x = (x_loc > 0.49) and (x_loc < 0.51); //if true x_loc is in the wall
				condition_y = not (y_loc > 0.425 and y_loc < 0.475) and not (y_loc > 0.525 and y_loc < 0.575); //if true y_loc is in the wall


				if (condition_x and condition_y == true) { //if (x_loc, y_loc) in wall
					V_mat(i, j) = v0;
				}
				else {
					V_mat(i, j) = 0.; //if (x_loc, y_loc) not in wall
				}
			}
		}
		V = makeVectorRel(V_mat); //transforms the result into a vector, such it can be used directly in finding the matricies A and B
	}





	void QuantumSlit::findPotentialTriple() { //finds potential for triple slit

		double x_loc;//local position x
		double y_loc; //local position y
		bool condition_x;
		bool condition_y;
		arma::mat V_mat = arma::mat(dim, dim); //metrix of potential

		for (int j = 0; j < dim; j++) { //sum over y
			y_loc = dx + j * dx;

			for (int i = 0; i < dim; i++) { //sum over x

				x_loc = dx + i * dx;


				condition_x = (x_loc > 0.49) and (x_loc < 0.51); //if true x_loc is in the wall
				condition_y = not (y_loc > 0.475 and y_loc < 0.525) and not (y_loc > 0.575 and y_loc < 0.625) and not (y_loc > 0.375 and y_loc < 0.425); //if true y_loc is in the wall


				if (condition_x and condition_y == true) { //if (x_loc, y_loc) in wall
					V_mat(i, j) = v0;
				}
				else { //if (x_loc, y_loc) not in wall
					V_mat(i, j) = 0.;
				}
			}
		}
		V = makeVectorRel(V_mat); //transforms the result into a vector, such it can be used directly in finding the matricies A and B
	}







	void QuantumSlit::findVectors() { //find vector a and b

		a = arma::cx_vec(N_elements); //vector a
		b = arma::cx_vec(N_elements); //vector b

		for (int i; i < N_elements; i++) { //sum over all indices k

			a.at(i) = 1. + 4. * r + arma::cx_double(0, dt / (double)2. * V.at(i));
			b.at(i) = 1. - 4. * r - arma::cx_double(0, dt / (double)2. * V.at(i));
		}
	}







	arma::cx_vec QuantumSlit::makeVectorComp(arma::sp_cx_dmat A) { //transforms a complex (M-2)x(M-2) metrix into a vector with (M-1)^2 elements

		arma::cx_vec w(N_elements); //the result vector

		for (int j = 0; j < dim; j++) { //sum over y
			for (int i = 0; i < dim; i++) { //sum over x
				w.at(j * dim + i) = A(i, j);
			}
		}

		return w;
	}




	arma::vec QuantumSlit::makeVectorRel(arma::mat A) { //transforms a real (M-2)x(M-2) metrix into a vector with (M-1)^2 elements

		arma::vec w(N_elements); //the result vector

		for (int j = 0; j < dim; j++) { //sum over y
			for (int i = 0; i < dim; i++) { //sum over x
				w.at(j * dim + i) = A(i, j);
			}
		}

		return w;
	}



