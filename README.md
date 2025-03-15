# Wavefunction interference

Code in this repository solves the Schrodinger equation numerically, using Crank Nicolson Scheme. The C++ code solves the equation and writes results, as real and imaginary parts of the wavefunction, or the probability distribution into text files. Python code reads the text files and plots the results.

## This repository contains the following code:

### C++:
- QuantumSlit: This is an object that creates either a single, double, or triple slit, creates an initial wavefunction, and finds the evolution of the system using Nicolson Scheme. There is both a public and private section. All attributes are public, together with some functions that are not called from within the object and are meant to be called from a main file. In the private section, there are only associated utility functions meant to be called only from within the class. The code uses the armadillo library, meaning it should be built with -larmadillo flag.

- utils: These files contain four functions. The first print_sp_matrix_structure, prints the structure of an armadillos sp-metrix. The three remaining are writeToFileReal, writeToFileImag, and writeToFil, all of which create text files and fill them with the results, which are respectively the real part of the wavefunction, the imaginary part of the wavefunction, and the probability distribution. This code does use armadillo library, meaning it should be built with -larmadillo flag.

- main.cpp: The main function which creates five instances of QuantumSlit object. The first solves the system without any barriers, the second and third with double slits, the fourth with a single slit, and the last with three slits. This function needs to be built together with both the QuantumSlit and utils. This code does use armadillo library, meaning it should be built with -larmadillo flag. An example of how to build and run this file in a Linux terminal using one command is:
    g++ main.cpp src/utils.cpp src/QuantumSlit.cpp -I include -larmadillo -O3 -o main.exe ; ./main.exe

### Python:
- plotting.py: This code interprets data collected by the above-mentioned code. Its only purpose is to plot the results.

### LaTeX:
- main.tex: This is the source code for the report. Running this code requires figures that are stored in the directory LaTeX/figures.
