#ifndef MAIN_STABILITY_HPP
#define MAIN_STABILITY_HPP

// Includes
#include <random>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <sys/stat.h>
#include "Tools.hpp"
#include "Plots.hpp"
#include "Matrix.hpp"
using ::std::cin;
using ::std::cout;
using ::std::endl;
using ::std::fixed;

// Simulation types
#define SIMTYPE_SW 1
#define SIMTYPE_SF 2
#define SIMTYPE_LE 3
#define SIMTYPE_RG 4
#define SIMTYPE_GA 5

// Switches
//#define DETAILED_OUTPUT // if defined, information about every trial is printed into a file
//#define MAT_OUTPUT // if defined, every single matrix is printed into a file (DETAILED_OUTPUT has to be defined)
//#define EV_OUTPUT // if defined, eigenvalues of every single matrix are printed into a file (DETAILED_OUTPUT has to be defined)
#define NONZEROS_EXACT // if defined, there is always an exactly defined number of non-zeros in the matrix (depending on the connectance)

#endif
