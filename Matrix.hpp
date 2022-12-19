#ifndef MATRIX_HPP
#define MATRIX_HPP

// Includes
#include <random>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "Tools.hpp"
using ::std::vector;
using ::std::ofstream;
using ::std::minstd_rand0;
using ::std::uniform_int_distribution;
using ::std::uniform_real_distribution;
using ::std::cout;
using ::std::endl;
using ::std::fixed;
using ::std::normal_distribution;

// Methods to determine eigenvalues
#define EV_METHOD_CLASSICAL 1
#define EV_METHOD_GSL 2

// Switches
//#define PAJEK_FORMAT // if defined, the file format for the program Pajek is used by Matrix::saveMatrix()

// Type definition for matrix elements with subtracted lambda
typedef struct 
{
	double element;
	int lambda;
} twoc;

/*** Matrix ***
 * Class that represents a matrix of coefficients for a linear system */
class Matrix 
{
private:
	int N; // matrix row dimension
	double* M; // the matrix, consisting of double values

	minstd_rand0 rg_conn; // linear congruential random generator for connection probability
	minstd_rand0 rg_strength; // linear congruential random generator for coupling strengths values to establish connections

	uniform_int_distribution<int> u_dist_element; // uniform distribution across matrix elements
	uniform_real_distribution<double> u_dist_conn; // uniform distribution for connection probability

	void chPolynomialRec(double*, int, twoc**);

public:
	/*** cNN (macro) ***
	 * Returns a consecutive number for element (i,j) rather than a pair of numbers like (i,j), be  *
	 * aware that it starts with zero, unlike i and j *
	 * - int i: the row where the element is located *
	 * - int j: the column where the element is located */
	#define cNN(i, j) (((i)-1)*N + ((j)-1))

	/*** cNN2 (macro) ***
	 * Returns a consecutive number for element (i,j) rather than a pair of numbers like (i,j), be  *
	 * aware that it starts with zero, unlike i and j (general case for a matrix with dimension d) *
	 * - int i: the row where the element is located *
	 * - int j: the column where the element is located *
	 * - int d: the dimension of the considered matrix */
	#define cNN2(i, j, d) (((i)-1)*d + ((j)-1))

	/*** row (macro) ***
	 * Returns the row number for element n (in consecutive numbering), be  *
	 * aware that it starts with one, unlike the consecutive number *
	 * - int n: the consecutive element number */
	#define row(n) ((((n) - ((n) % N)) / N) + 1)

	/*** row2 (macro) ***
	 * Returns the row number for element n (in consecutive numbering), be  *
	 * aware that it starts with one, unlike the consecutive number (general case for a matrix with dimension d) *
	 * - int n: the consecutive element number *
	 * - int d: the dimension of the considered matrix */
	#define row2(n, d) ((((n) - ((n) % d)) / d) + 1)

	/*** col (macro) ***
	 * Returns the column number for element n (in consecutive numbering), be  *
	 * aware that it starts with one, unlike the consecutive number *
	 * - int n: the consecutive element number */
	#define col(n) (((n) % N) + 1)

	/*** col2 (macro) ***
	 * Returns the column number for element n (in consecutive numbering), be  *
	 * aware that it starts with one, unlike the consecutive number (general case for a matrix with dimension d) *
	 * - int n: the consecutive element number *
	 * - int d: the dimension of the considered matrix */
	#define col2(n, d) (((n) % d) + 1)

	/*** symm (macro) ***
	 * Returns the number of the symmetric element for an element given  *
	 * by its consecutive number *
	 * - int n: the consecutive element number */
	#define symm(n) (cNN(col(n),row(n)))

	/*** symm2 (macro) ***
	 * Returns the number of the symmetric element for an element given  *
	 * by its consecutive number (general case for a matrix with dimension d) *
	 * - int n: the consecutive element number *
	 * - int d: the dimension of the considered matrix */
	#define symm2(n, d) (cNN2(col2(n, d),row2(n, d), d))

	int nonzeroTest();
	void reset();
	void chPolynomial(double*);
	bool eigenvalues(complex<double>*, int, ofstream* = NULL);
	int isStable(ofstream* = NULL);
	int isStableRH(ofstream* = NULL);
	bool isOnMainDiagonal(int) const;
	bool areConnected(int, int) const;
	int sumDegreeDist(int*) const;
	vector<node_and_degree> getOutHubs(int) const;
	node_and_degree getLargestOutHub() const;
	node_and_degree getLargestInHub() const;
	void kronecker(Matrix A, Matrix B);
	void saveMatrix(ofstream *f) const;

	void generateGaussRandom(const double, const double, const double, const double, const double, ofstream* = NULL);
	void generateUniformRandom(const double, const double, const double, const double, const double, 
		                       int = 2, ofstream* = NULL);
	void generateUniformRandomNonExact(const double, const double, const double, const double, const double, int = 2, ofstream* = NULL);
	void generateSmallWorld(const int, int, const double, const double, const double, const double, bool = false, ofstream* = NULL);
	void generateSmallWorld_fullRewire(const int, const int, const double, const double, const double, const double, ofstream* = NULL);

	void generateScaleFree(const int, const int, const double, const double, const double, const double, int = 0, bool = false, ofstream* = NULL);
	int generateLeskovec(const int, const int, const double, const double, const double, const double, bool = false, ofstream* = NULL);

	double trace() const;
	void add(Matrix);
	void multiply(Matrix);
	void multiply(double);
	void mpow(int);
	void transpose();
	double get(int) const;
	bool isNonZero(int) const;
	bool isZero(int) const;
	void set(double*);
	void set(int, double);
	int dim() const;
	void init(int);

	Matrix(int);
	Matrix();
	Matrix(const Matrix&);
	~Matrix();
	Matrix& operator=(const Matrix&);
	Matrix& operator+=(const Matrix&);
};


#endif
