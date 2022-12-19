#ifndef TOOLS_HPP
#define TOOLS_HPP

// Includes
#include <chrono>
#include <cmath>
#include <cstring>
//#include <time.h> // may be needed on certain platforms
#include <errno.h>
#include <complex>
#include <iostream>
#include <fstream>
using ::std::string;
using ::std::complex;
using ::std::ofstream;
using ::std::chrono::system_clock;
using ::std::ios;
using ::std::cout;
using ::std::endl;
using ::std::ifstream;

// Constants
#define EPSILON 1e-12	// numeric constant, very small number that is considered as zero (e.g., for Newton's method)
#define EPSILON2 1e-9	// numeric constant for residue in long division
#define MAX_ITERATIONS	1000000 // maximal number of iterations for Newton's method

// Switches
#define ENFORCE_CONVERGENCE // if defined, findZeros() will try harder (can take much longer but offers a much higher chance to find all zeros)

/*** node_and_degree ***
 * Class holding a node and its degree */
class node_and_degree
{
public:
	int node;
	int degree;
	node_and_degree(int, int);
};

// Function declarations
string concat(const char*, const char *);
string concat(string, const char*);
string concat(const char*, string);
bool sortvec(const node_and_degree&, const node_and_degree&);
string dtos(double, int, bool = false);
string dateStr(string = "", bool = false);
void writeLineToFile(string, ofstream*);
int sgn (double);
int timeMeasure (bool);
void showChDirErrMessage();
unsigned int getClockSeed();
bool copyFile(string, string);
int circle_dist(int, int, int);
void derivePolynomial(complex<double>*, int);
void derivePolynomial(double*, int);
bool polynomialLongDivision(complex<double>*, complex<double>, int n, ofstream* = NULL);
int findZeros(complex<double>*, double*, int, ofstream* = NULL);

#endif
