/***************************************************************
 *** Matrix class to describe networks of different topology ***
 ***             (c) Jannik Luboeinski 2016-2022             ***
 ***************************************************************/

#include "Matrix.hpp"

/*** Matrix::chPolynomialRec ***
 * Returns the coefficients of the characteristic polynomial of a given matrix by argument pv *
 * - pv: the polynomial coefficients (coefficient for x^i with array index i), has to be an allocated array of dim+1 numbers *
 * - dim: the dimension of the matrix *
 * - matrix: the matrix (actually two matrizes for degree zero and degree one coefficients) */
void Matrix::chPolynomialRec(double* pv, int dim, twoc** matrix)
{
	twoc** submatrix;
	double* pv2;

	// Exit condition for recursion
	if (dim == 1)
	{
		pv[0] = matrix[0][0].element;
		pv[1] = matrix[0][0].lambda;

		return;
	}

	pv2 = new double[dim];
	submatrix = new twoc* [dim-1];

	// Clear vector of coefficients and allocate submatrix
	for(int i=0; i <= dim; i++)
	{
		pv[i] = 0.;
		if (i < dim-1)
			submatrix[i] = new twoc [dim-1];
	}

	// Laplace's formula applied to the first row
	for(int i=1; i <= dim; i++) // loop over submatrices
	{
		int a=0, b=0;

		// compute (dim-1)x(dim-1) submatrix (leave out first row and i-th column)
		for(int j=1; j <= dim; j++)
		{
			for(int k=1; k <= dim; k++)
			{
				if ((j != 1) && (k != i))
				{
					submatrix[a][b].element = matrix[j-1][k-1].element;
					submatrix[a][b++].lambda = matrix[j-1][k-1].lambda;
				}
			}
			if (j != 1)
				a++;
			b = 0;
		}

		// Recursion
		chPolynomialRec(pv2, dim-1, submatrix);
		for(int l=0; l<dim; l++)
		{
			pv[l] += pow((-1),(i+1)) * matrix[0][i-1].element * pv2[l];
			pv[l+1] += pow((-1),(i+1)) * matrix[0][i-1].lambda * pv2[l];
		}
		// det += Math.pow((-1),(i+1)) * matrix[0][i-1] * Det(submatrix, n-1);
	}
	delete[] pv2;

	// Free space allocated for submatrix
	for(int i=0; i<dim-1; i++)
	{
		delete[] submatrix[i];
	}
	delete[] submatrix;
}

/*** Matrix::nonzeroTest ***
 * Check if exactly the desired number of non-zeros was generated *
 * - return: the number of non-zero values except those on the main diagonal */
int Matrix::nonzeroTest()
{
	int nonzeros = 0;

	for (int i=1; i<=N; i++)
	{
		for (int j=1; j<=N; j++)
		{
			if (i != j)
			{
				if (abs(M[cNN(i,j)]) > EPSILON)
					nonzeros++;
			}
		}
	}

	return nonzeros;
}

/*** Matrix::reset ***
 * Resets the generation of random numbers and all matrix values */
void Matrix::reset()
{
	for (int i=0; i<N*N; i++) // set all entries to zero
	{
		M[i] = 0.;
	}
}


/*** Matrix::chPolynomial ***
 * Returns the coefficients of the characteristic polynomial of the given matrix by argument pv, using the function chPolynomialRec *
 * - pv: the polynomial coefficients (coefficient for x^i with array index i), has to be an allocated array of N+1 numbers */
void Matrix::chPolynomial(double* pv)
{
	twoc** matrix;

	// Copy matrix "M" into new object "matrix" that is able to contain variables in matrix elements;
	//	subtract lambda from all elements on the main diagonal
	matrix = new twoc* [N];

	for(int i=0; i<N; i++)
	{
		matrix[i] = new twoc [N];
		for(int j=0; j<N; j++)
		{
			matrix[i][j].element = -1. * M[cNN(i+1,j+1)];
			if (i == j)
				matrix[i][j].lambda = 1; // marker for "1*lambda"
			else
				matrix[i][j].lambda = 0;
		}
	}

	// Start computation via recursion
	chPolynomialRec(pv, N, matrix);

	// Free space allocated for "matrix"
	for(int i=0; i<N; i++)
	{
		delete[] matrix[i];
	}
	delete[] matrix;
}

/*** Matrix::eigenvalues ***
 * Computes the eigenvalues of the matrix employing Newton's method to solve *
 * the characteristic polynomial *
 * - ev: array of complex eigenvalues that has been allocated before *
 * - method: specifies the method to be used to determine the eigenvalues (EV_METHOD_CLASSICAL for solving using char. pòlynomial/quite slow, EV_METHOD_GSL for solving using GSL which uses Francis QR double-shift) *
 * - err_output [optional]: pointer to a file for error output *
 * - return: true if all eigenvalues were found, false if not */
bool Matrix::eigenvalues(complex<double>* ev, int method, ofstream* err_output)
{
	// find the characteristic polynomial, then solve for its zeros using Horner's method (mind that finding the char. polynomial is slow !)
	if (method == EV_METHOD_CLASSICAL)
	{
		double* chp_coeff = new double [N+1]; // coefficients of the characteristic polynomial

		chPolynomial(chp_coeff);

		if ( N > 0 && findZeros(ev, chp_coeff, N, err_output) == N ) // if there are as many zeros of the char. polynomial as dimensions
		{
			delete[] chp_coeff;
			return true;
		}

		delete[] chp_coeff;
	}

	// use Francis' QR double-shift algorithm as implemented in the GSL
	else if  (method == EV_METHOD_GSL)
	{
		gsl_matrix* M_buf = gsl_matrix_alloc(N, N); // allocate new GSL matrix for eigenvalue computation
		gsl_matrix_view m = gsl_matrix_view_array(M, N, N);
		gsl_matrix_memcpy(M_buf, &m.matrix); // use buffer M_buf to retain the original matrix
		gsl_vector_complex *eval = gsl_vector_complex_alloc(N);
		gsl_eigen_nonsymm_workspace *w = gsl_eigen_nonsymm_alloc(N);

		gsl_eigen_nonsymm(M_buf, eval, w); // use buffer M_buf instead of view on original matrix &m.matrix

		if (w->n_evals == N) // if there are as many eigenvalues as dimensions
		{
			memcpy(ev, eval->data, N*sizeof(complex<double>));

			gsl_eigen_nonsymm_free(w);
			gsl_vector_complex_free(eval);
			gsl_matrix_free(M_buf);
			return true;
		}

		gsl_eigen_nonsymm_free(w);
		gsl_vector_complex_free(eval);
		gsl_matrix_free(M_buf);
	}

	return false;
}

/*** Matrix::isStable ***
 * Tells by the eigenvalues if the system corresponding to the matrix is stable; *
 * stable in this context means that the real parts of all eigenvalues are negative *
 * (cf. Ashby, Design for a brain, p. 255) *
 * - err_output [optional]: pointer to a file for error output *
 * - return: 1 if stable, -1 if unstable and 0 if stability could not be determined */
int Matrix::isStable(ofstream* err_output)
{
	complex<double>* ev = new complex<double>[N]; // array of eigenvalues

	if (!eigenvalues(ev, EV_METHOD_GSL, err_output))
	{
		delete[] ev;
		return 0;
	}

	for (int i=0; i<N; i++)
	{
		if (ev[i].real() > -EPSILON)
		{
			delete[] ev;
			return -1;
		}
	}

	delete[] ev;
	return 1;
}

/*** Matrix::isStableRH ***
 * Tells by using the Routh-Hurwitz criterion if the system corresponding to the matrix is stable; *
 * stable in this context means that the real parts of all eigenvalues are negative *
 * (cf. Ashby, Design for a brain, p. 255) *
 * - err_output [optional]: pointer to a file for error output *
 * - return: 1 if stable, -1 if unstable and 0 if stability could not be determined */
int Matrix::isStableRH(ofstream* err_output)
{
	double* chp_coeff = new double [N+1]; // coefficients of the characteristic polynomial
	double** hurwitz = new double* [N+1]; // Hurwitz matrix
	for (int i=0; i<N+1; i++)
		hurwitz[i] = new double [N+1];

	chPolynomial(chp_coeff); // obtain char. polynomial (computation is slow...)

	for (int i=0; i<N+1; i++)
	{
		bool zitr = false; // indicates if a zero has occurred in the current row
		bool nzitr = false; // indicates if a non-zero has occurred in the current row

		for (int j=0; j<N+1; j++)
		{
			hurwitz[i][j] = 0.; // by default, set all entries to zero

			if (i == 0 && j <= N/2) // the first row
			{
				hurwitz[0][j] = chp_coeff[N - (2*j)];
			}
			else if (i == 1 && j <= (N-1)/2) // the second row
			{
				hurwitz[1][j] = chp_coeff[N - (2*j+1)];
			}
			else if (i >= 2 && (j+1) <= N) // all further rows
			{
				hurwitz[i][j] = ( hurwitz[i-1][0] * hurwitz[i-2][j+1] - hurwitz[i-2][0] * hurwitz[i-1][j+1] ) / hurwitz[i-1][0];
			}

			if (i > 0 && j == 0) // detect sign changes
			{
				if (sgn(hurwitz[i][0]) != sgn(hurwitz[i-1][0])) // sign change -- instability encountered
				{
					if (err_output)
						*err_output << "Sign change in row " << i << " (" << hurwitz[i-1][0] <<  " -> " << hurwitz[i][0] << ")" << endl;

					for (int i=0; i<N+1; i++)
						delete[] hurwitz[i];
					delete[] hurwitz;
					delete[] chp_coeff;

					return -1;
				}
			}

			if (abs(hurwitz[i][j]) < EPSILON) // detect zeros
			{
				zitr = true;

				if (!nzitr && j == N) // if there were only zeros in this row
				{
					if (err_output)
						*err_output << "Row " << i << " has only zeros!" << endl << "Auxiliary polynomial: ";


					// Use auxiliary polynomial, which is located in the row above
					if (i > 0)
					{
						double* aux_poly = new double [N+1];
						memset(aux_poly, 0., (N+1)*sizeof(double)); // set all values to zero

						for (int k=0; k <= (N-(i-1))/2; k++) // obtain auxiliary polynomial
						{
							int index = N-(i-1) - 2*k;
							aux_poly[index] = hurwitz[i-1][k];
						}
						if (err_output)
						{
							for (int k=N;k>0;k--)
								*err_output << aux_poly[k] << " x^" << k << " + ";
							*err_output << aux_poly[0] << endl << endl;
						}

						derivePolynomial(aux_poly, N-(i-1)); // derive auxiliary polynomial

						for (int k=0; k < N+1; k++) // insert derived auxiliary polynomial into Hurwitz matrix
						{
							if (k <= (N-i)/2)
							{
								int index = N-i - 2*k;
								hurwitz[i][k] = aux_poly[index];
							}
						}

						delete[] aux_poly;
					}
					else // there is no row above -- stability cannot be determined
					{
						if (err_output)
							*err_output << "N/A - now row above!" << endl;

						for (int i=0; i<N+1; i++)
							delete[] hurwitz[i];
						delete[] hurwitz;
						delete[] chp_coeff;

						return 0;
					}
				}
			}
			else if (zitr) // non-zero now but zero before in this row -- instability encountered
			{
				if (err_output)
					*err_output << "Zero in non-zero row " << i << "!" << endl;

				for (int i=0; i<N+1; i++)
					delete[] hurwitz[i];
				delete[] hurwitz;
				delete[] chp_coeff;

				return -1;
			}
			else
				nzitr = true;
		}
	}

	// Output of Hurwitz matrix
	/*if (err_output)
	{
		for(int i=0; i<N+1; i++)
		{
			for(int j=0; j<N+1; j++)
			{
				if (hurwitz[i][j] < 0)
					*err_output << fixed << "\t" << hurwitz[i][j]; // has a minus sign
				else
					*err_output << fixed << "\t " << hurwitz[i][j]; // leave a space for the sign
			}
			*err_output << endl;
		}
		*err_output << endl;
	}*/

	for (int i=0; i<N+1; i++)
		delete[] hurwitz[i];
	delete[] hurwitz;
	delete[] chp_coeff;
	return 1;
}

/*** Matrix::isOnMainDiagonal ***
 * Returns whether or not a given element is on the main diagonal *
 * - int n: the number of the element in consecutive ordering *
 * - return: true if element is on main diagonal, false if not */
bool Matrix::isOnMainDiagonal(int n) const
{
	if (row(n) == col(n))
		return true;
	else
		return false;
}


/*** Matrix::areConnected ***
 * Returns whether or not there is a coupling  between element m and element n *
 * - int m: the number of the first element in consecutive ordering *
 * - int n: the number of the second element in consecutive ordering *
 * - return: true if connection from m to n exists, false if not */
bool Matrix::areConnected(int m, int n) const
{
	if (M[cNN(m+1,n+1)] != 0)
		return true;
	else
		return false;
}

/*** Matrix::sumDegreeDist ***
 * Gets the degree distribution of the matrix and adds it to a provided array (which is NOT reset by this function). *
 * The method works only for undirected or bidirectional graphs. *
 * int* dist: array of N+1 elements containing the number of vertices for each degree *
 * return: the sum of the degrees of all vertices */
int Matrix::sumDegreeDist(int* dist) const
{
	int tot_degree = 0; // the sum of all degrees

	for(int i=1; i<=N; i++) // pass through rows
	{
		int nonzeros = 0; // number of nonzeros in the current row

		for(int j=1; j<=N; j++) // pass through columns
		{
			if (abs(M[cNN(i,j)]) > EPSILON) // nonzero found
				nonzeros++;
		}

		dist[nonzeros]++; // increase counter for found degree
		tot_degree += nonzeros; // increase total degree
	}

	return tot_degree;
}

/*** Matrix::getOutHubs ***
 * Returns the largest out-degree hubs in the graph. *
 * num_hubs: number of hubs to be returned (1 would only be the largest, 2 would also include the second-largest, and so on)
 * return: vector of node_and_degree objects containing the index degree for each hub */
vector<node_and_degree> Matrix::getOutHubs(int num_hubs) const
{
	vector<node_and_degree> hubs;
	
	if (num_hubs <= 0)
	{
		return vector<node_and_degree>();
	}
	
	for(int i=1; i<=N; i++) // pass through rows
	{
		int nonzeros = 0; // number of nonzeros in the current row

		for(int j=1; j<=N; j++) // pass through columns and count nonzeros
		{
			if (abs(M[cNN(i,j)]) > EPSILON) // nonzero found
				nonzeros++;
		}
		hubs.push_back(node_and_degree(i, nonzeros));
	}
	sort(hubs.begin(), hubs.end(), sortvec); // sort the vector of hubs using the sortvec() function
	
	hubs.erase(hubs.begin() + num_hubs + 1, hubs.end() + 1); // keep the 'num_hubs' largest hubs and remove all smaller ones
	
	return hubs;
}

/*** Matrix::getLargestOutHub ***
 * Returns the largest hub in the graph (the node with the largest out-degree). *
 * return: the index of the largest out-degree hub, and its degree, as a node_and_degree object */
node_and_degree Matrix::getLargestOutHub() const
{
	int k_max_out = 0; // highest out-degree value that occurs (maximally possible: N)
	int hub_out = -1; // the largest out-degree hub in the graph
	
	for(int i=1; i<=N; i++) // pass through rows
	{
		int nonzeros = 0; // number of nonzeros in the current row

		for(int j=1; j<=N; j++) // pass through columns
		{
			if (abs(M[cNN(i,j)]) > EPSILON) // nonzero found
				nonzeros++;
		}

		if (nonzeros > k_max_out)
		{
			hub_out = i; // the currently considered node is now the out-degree hub
			k_max_out = nonzeros;
		}
	}

	
	return node_and_degree(hub_out, k_max_out);
}

/*** Matrix::getLargestInHub ***
 * Returns the largest hub in the graph (the node with the largest in-degree). *
 * return: the index of the largest in-degree hub, and its degree, as a node_and_degree object */
node_and_degree Matrix::getLargestInHub() const
{
	int k_max_in = 0; // highest out-degree value that occurs (maximally possible: N)
	int hub_in = -1; // the largest out-degree hub in the graph 
	
	for(int i=1; i<=N; i++) // pass through columns
	{
		int nonzeros = 0; // number of nonzeros in the current column

		for(int j=1; j<=N; j++) // pass through rows
		{
			if (abs(M[cNN(j,i)]) > EPSILON) // nonzero found
				nonzeros++;
		}

		if (nonzeros > k_max_in)
		{
			hub_in = i; // the currently considered node is now the in-degree hub
			k_max_in = nonzeros;
		}
	}

	return node_and_degree(hub_in, k_max_in);
}


/*** Matrix::kronecker ***
 * Makes this matrix the Kronecker product of a matrix A with a matrix B
 * - A: first matrix *
 * - B: second matrix */
void Matrix::kronecker(Matrix A, Matrix B)
{
	int dimA = A.dim();
	int dimB = B.dim();
	int dimAB = dimA*dimB;
	Matrix ret(dimAB);

	for (int i=1; i<=dimA; i++) // go through rows of A
	{
		int rowblock = (i-1)*dimB;
		for (int k=1; k<=dimB; k++) // go through rows of B
		{
			for (int j=1; j<=dimA; j++) // go through columns of A
			{
				int colblock = (j-1)*dimB;
				for (int l=1; l<=dimB; l++) // go through columns of B
				{
					double Ag = A.get(cNN2(i,j,dimA));
					double Bg = B.get(cNN2(k,l,dimB));
					ret.set(cNN2(rowblock+k, colblock+l, dimAB), Ag*Bg);
				}
			}
		}
	}

	*this = ret;
}

/*** Matrix::saveMatrix ***
 * Saves all the entries of the matrix to a given file *
 * - ofstream* f: pointer to a file for information output */
void Matrix::saveMatrix(ofstream *f) const
{
	if (f == NULL)
		return;

	*f << "Matrix (N = " << N << "):" << endl << endl;
	f->precision(3);

	for(int i=1; i<=N; i++)
	{
		int nonzeros = N;

		for(int j=1; j<=N; j++)
		{
			if (M[cNN(i,j)] < 0)
			{
				*f << fixed
#ifndef PAJEK_FORMAT
				   << "\t"
#else
				   << " "
#endif
				   << M[cNN(i,j)]; // has a minus sign
			}
			else
			{
				if (abs(M[cNN(i,j)]) < EPSILON)
					nonzeros--;

				*f << fixed
#ifndef PAJEK_FORMAT
				   << "\t " // leave a space for the (not written) plus sign
#else
				   << " "
#endif
				   << M[cNN(i,j)]; // has a minus sign
			}
		}

		*f
#ifndef PAJEK_FORMAT
		   << "\t(" << nonzeros << ")" // write number of non-zeros at the end of the line
#endif
		   << endl;
	}

	*f << endl;
}

/*** Matrix::generateGaussRandom ***
 * Generates Gaussian random values for the whole matrix without any pattern *
 * - double _C: probability of coupling ("connectance" as in Gardner & Ashby 1970) *
 * - double _avS: average coupling strength *
 * - double _sigmaS: standard deviation of coupling strength *
 * - double _avS_self: average coupling strength for self-coupling *
 * - double _sigmaS_self: standard deviation of coupling strength for self-coupling *
 * - ofstream* f [optional]: pointer to a file for information output */
void Matrix::generateGaussRandom(const double _C, const double _avS, const double _sigmaS, const double _avS_self, const double _sigmaS_self, ofstream* f)
{
	minstd_rand0 rg_nonzeros(getClockSeed()); // linear congruential random generator for non-zero values to establish connections
	minstd_rand0 rg_strength(getClockSeed()); // linear congruential random generator for coupling strengths values to establish connections
	uniform_int_distribution<int> u_dist_nonzeros(0, N*N-1); // uniform distribution for non-zero values
	normal_distribution<double> n_dist_strength(_avS,_sigmaS); // normal distribution for coupling strengths
	normal_distribution<double> n_dist_strength_self(_avS_self,_sigmaS_self); // normal distribution for coupling strengths for self-coupling
	int nonzeros = round((N*N - N) * _C); // number of matrix entries that are non-zero (except self-couplings)

	if (f != NULL)
	{
		*f << "Connectance = " << _C << endl;
		*f << "Coupling strength = " << _avS << " +- " << _sigmaS << endl;
		*f << "Self-coupling strength = " << _avS_self << " +- " << _sigmaS_self << endl << endl;
	}

	reset(); // to obtain an empty matrix

	while (nonzeros > 0)
	{
		int choice = u_dist_nonzeros(rg_nonzeros);

		if (row(choice) != col(choice) // if element is not on the main diagonal (i.e., none of the self-couplings)
			&& M[choice] == 0) // and if element has not been set yet
		{
			double rn = n_dist_strength(rg_strength);
			if (abs(rn) > EPSILON) // do not accept zero values
			{
				M[choice] = rn; // draw coupling strength
				nonzeros--;
			}
		}
	}

	for (int i=1; i<=N; i++) // elements on the main diagonal
	{
		M[cNN(i,i)] = n_dist_strength_self(rg_strength); // draw coupling strength from special distribution
	}

	if (f != NULL)
		*f << "Non-zeros desired / obtained: " << int(round((N*N - N) * _C)) << " / " << nonzeroTest() << endl;
}

/*** Matrix::generateUniformRandom ***
 * Generates uniform random values for the whole matrix without any pattern but for an exact number of non-zero values *
 * - double _C: probability of coupling ("connectance" as in Gardner & Ashby 1970) *
 * - double _minS: minimum coupling strength *
 * - double _maxS: maximum coupling strength *
 * - double _minS_self: minimum coupling strength for self-coupling *
 * - double _maxS_self: standard deviation of coupling strength for self-coupling *
 * - int graph_struct [optional]: indicates if the graph shall be undirected (0), asymmetrically bidirectional (1), or freely directed (2) *
 * - ofstream* f [optional]: pointer to a file for information output */
void Matrix::generateUniformRandom(const double _C, const double _minS, const double _maxS, const double _minS_self, const double _maxS_self, 
                                   int graph_struct, ofstream* f)
{
	uniform_real_distribution<double> u_dist_strength(_minS,_maxS); // uniform distribution for coupling strengths
	uniform_real_distribution<double> u_dist_strength_self(_minS_self,_maxS_self); // uniform distribution for coupling strengths for self-coupling

	int nonzeros = round((N*N - N) * _C); // number of matrix entries that are non-zero (except self-couplings)

	if (f != NULL)
	{
		*f << "Connectance = " << _C << endl;
		*f << "Coupling strength in [" << _minS << ":" << _maxS << "]" << endl;
		*f << "Self-coupling strength in [" << _minS_self << ":" << _maxS_self << "]" << endl << endl;
	}

	reset(); // to obtain an empty matrix

	while (nonzeros > 0)
	{
		int choice = u_dist_element(rg_conn);

		if (row(choice) != col(choice) // if element is not on the main diagonal (i.e., none of the self-couplings)
			&& isZero(choice)) // and if element has not been set yet
		{
			double coupl_strength = u_dist_strength(rg_strength); // draw coupling strength
			double coupl_strength2 = _minS;

			if (graph_struct == 0) // if the graph is undirected, copy the connection strength from the opposite (upper or lower) triangle
				coupl_strength2 = coupl_strength; // set to the same value
			else if (graph_struct == 1) // if the graph is bidirectional but asymmetric, draw a new value
				coupl_strength2 = u_dist_strength(rg_strength); // set to a new value
			
			if (abs(coupl_strength) > EPSILON && abs(coupl_strength2) > EPSILON) // if the drawn coupling strength(s) unequal zero
			{
				M[choice] = coupl_strength; // set to drawn value
				nonzeros--;
				if (graph_struct < 2) // if undirected or bidirectional graph is considered
				{
					M[symm(choice)] = coupl_strength2; // set to selected value
					nonzeros--;
				}
			}
		}
	}

	for (int i=1; i<=N; i++) // elements on the main diagonal
	{
		M[cNN(i,i)] = u_dist_strength_self(rg_strength); // draw coupling strength from special distribution
	}

	if (f != NULL)
		*f << "Non-zeros desired / obtained: " << int(round((N*N - N) * _C)) << " / " << nonzeroTest() << endl;
}


/*** Matrix::generateUniformRandomNonExact ***
 * Generates uniform random values for the whole matrix without any pattern for an average number of non-zero values *
 * - double _C: probability of coupling ("connectance" as in Gardner & Ashby 1970) *
 * - double _minS: minimum coupling strength *
 * - double _maxS: maximum coupling strength *
 * - double _minS_self: minimum coupling strength for self-coupling *
 * - double _maxS_self: standard deviation of coupling strength for self-coupling *
 * - int graph_struct [optional]: indicates if the graph shall be undirected (0), asymmetrically bidirectional (1), or freely directed (2) *
 * - ofstream* f [optional]: pointer to a file for information output */
void Matrix::generateUniformRandomNonExact(const double _C, const double _minS, const double _maxS, const double _minS_self, const double _maxS_self, 
                                           int graph_struct, ofstream* f)
{
	uniform_real_distribution<double> u_dist_strength(_minS,_maxS); // uniform distribution for coupling strengths
	uniform_real_distribution<double> u_dist_strength_self(_minS_self,_maxS_self); // uniform distribution for coupling strengths for self-coupling

	if (f != NULL)
	{
		*f << "Connectance = " << _C << endl;
		*f << "Coupling strength in [" << _minS << ":" << _maxS << "]" << endl;
		*f << "Self-coupling strength in [" << _minS_self << ":" << _maxS_self << "]" << endl << endl;
	}

	reset(); // to obtain an empty matrix

	for (int i=0; i<N*N; i++)
	{
		int irow = row(i), icol = col(i);
		if (irow == icol) // elements on the main diagonal
		{
			M[cNN(irow,icol)] = u_dist_strength_self(rg_strength); // draw coupling strength from special distribution
		}
		else if (graph_struct > 1            // if the graph is directed, and the element is not on the main diagonal (i.e., none of the self-couplings),
		      && u_dist_conn(rg_conn) <= _C) // and there shall be a connection                              
		{
			M[cNN(irow,icol)] = u_dist_strength(rg_strength); // draw coupling strength
			
		}
		else if (graph_struct <= 1           // if the graph is undirected or asymmetrically bidirectional, 
		      && icol > irow                 // and the element is in the upper triangle and not on the main diagonal,
		      && u_dist_conn(rg_conn) <= _C) // and there shall be a connection
		{
			M[cNN(irow,icol)] = u_dist_strength(rg_strength); // draw coupling strength
			if (graph_struct == 0) // if the graph is undirected
				M[cNN(icol,irow)] = M[cNN(irow,icol)]; // copy the connection strength from the opposite (upper or lower) triangle
			else if (graph_struct == 1) // if the graph is bidirectional but asymmetric
				M[cNN(icol,irow)] = u_dist_strength(rg_strength); // draw new coupling strength
		}
	}

	if (f != NULL)
		*f << "Non-zeros desired / obtained: " << int(round((N*N - N) * _C)) << " / " << nonzeroTest() << endl;
}

/*** Matrix::generateSmallWorld ***
 * Generates a small world graph of the Watts-Strogatz (1998) type, optionally with additional self-couplings. *
 * Therefore, at first a regular graph with connections from each vertex to its k nearest neighbors is created. *
 * Next, a specified number of edges is rewired to new, random destinations. *
 * - int k: number of nearest neighbors to which connections shall be established (should be an even number) *
 * - int rew_num: number of edges to be rewired *
 * - double _minS: minimum coupling strength *
 * - double _maxS: maximum coupling strength *
 * - double _minS_self: minimum coupling strength for self-coupling *
 * - double _maxS_self: standard deviation of coupling strength for self-coupling *
 * - bool is_directed [optional]: indicates if the graph shall be directed (i.e., the matrix will not be symmetric) *
 * - ofstream* f [optional]: pointer to a file for information output */
void Matrix::generateSmallWorld(const int k, int rew_num, const double _minS, const double _maxS, const double _minS_self, const double _maxS_self, 
                                bool is_directed, ofstream* f)
{
	uniform_int_distribution<int> u_dist_edges(0, k*N/2-1); // uniform distribution across all edges (except for loops)
	uniform_int_distribution<int> u_dist_end(1,2); // uniform distribution for edges end point (binary)
	uniform_int_distribution<int> u_dist_vertex(1,N); // uniform distribution for vertices in one row/column

	uniform_real_distribution<double> u_dist_strength(_minS,_maxS); // uniform distribution for coupling strengths
	uniform_real_distribution<double> u_dist_strength_self(_minS_self,_maxS_self); // uniform distribution for coupling strengths for self-coupling

	int* edge_pos = new int[k*N/2]; // array that contains the positions of edges in the upper triangle of the adjacency matrix (except for loops)
	int en = 0; // index for writing in edge_pos array

	if (f != NULL)
	{
		*f << "Nearest-neighbor connections = " << k << endl;
		*f << "Rewired connections = " << rew_num << endl;
		*f << "Coupling strength in [" << _minS << ":" << _maxS << "]" << endl;
		*f << "Self-coupling strength in [" << _minS_self << ":" << _maxS_self << "]" << endl << endl;
	}

	reset(); // to obtain an empty matrix

	// Generate regular graph
	for (int i=0; i<N*N; i++)
	{
		int irow = row(i), icol = col(i);
		if (irow == icol) // elements on the main diagonal
		{
			M[cNN(irow,icol)] = u_dist_strength_self(rg_strength); // draw coupling strength from special distribution
		}
		else if (icol > irow && circle_dist(icol,irow,N) <= k/2) // if element is in the upper triangle and not on the main diagonal,
		                                                         // and if it is a connection between near neighbors, then establish connection 
		{
			M[cNN(irow,icol)] = u_dist_strength(rg_strength); // draw coupling strength

			if (!is_directed) 
			{
				M[cNN(icol,irow)] = M[cNN(irow,icol)]; // set the same coupling strength to the corresponding element in the lower triangle
			}
			else
			{
				M[cNN(icol,irow)] = u_dist_strength(rg_strength); // draw coupling strength for the corresponding element in the lower triangle
			}
			if (isZero(cNN(irow,icol)) || isZero(cNN(icol,irow))) // do not accept zero values
			{
				i--;
				continue; // go for another random number(s)
			}

			edge_pos[en++] = cNN(irow,icol); // save position of this edge in this matrix
		}
	}

	// Rewiring (differences to Watts & Strogatz: there is no preference for changing small-distance edges; an exact number of edges is rewired)
	int rewire = rew_num; // number of edges that are to be rewired
	int it = 0;
	while (rewire > 0 && it++ < MAX_ITERATIONS)
	{
		int edge_num = u_dist_edges(rg_conn); // the number of the edge that will be rewired (self-couplings are left untouched)
		int edge = edge_pos[edge_num]; // position of the edge in the matrix

		if (edge < -EPSILON) // this edge has been rewired before
			continue;

		int end = u_dist_end(rg_conn); // the end (of the edge) that will be detached
		int new_target = u_dist_vertex(rg_conn); // the new target that the edge will connect to

		int edge_col_rewired = cNN(row(edge),new_target); // adjacency matrix index of the edge whose column value has changed
		int edge_row_rewired = cNN(new_target,col(edge)); // adjacency matrix index of the edge whose row value has changed

		if ( ( end == 1 || isNonZero(edge_row_rewired) ) // if end = 1, or if the edge for end = 2 already exists,
			&& isZero(edge_col_rewired) ) // AND if the new connection does not exist yet -- leave the row of the edge fixed
		{
			M[edge_col_rewired] = M[edge]; // create new connection with same coupling strength
			M[symm(edge_col_rewired)] = M[symm(edge)]; // create symmetric counterpart of the new connection
			M[edge] = 0.; // delete old connection
			M[symm(edge)] = 0.; // delete symmetric counterpart of the old connection

			rewire--; // one rewiring completed
			edge_pos[edge_num] = -1.; // mark as rewired
		}
		else if ( isZero(edge_row_rewired) ) // else if the new connection does not exist yet -- leave the column of the edge fixed
		{
			M[edge_row_rewired] = M[edge]; // create new connection with same coupling strength
			M[symm(edge_row_rewired)] = M[symm(edge)]; // create symmetric counterpart of the new connection
			M[edge] = 0.; // delete old connection
			M[symm(edge)] = 0.; // delete symmetric counterpart of the old connection

			rewire--; // one rewiring completed
			edge_pos[edge_num] = -1.; // mark as rewired
		}

	}

	if (f != NULL)
	{
		if (it >= MAX_ITERATIONS)
			*f << "Only " << rew_num-rewire << " rewirings were performed." << endl;
		*f << "Non-zeros desired / obtained: " << k*N << " / " << nonzeroTest() << endl;
	}

	delete[] edge_pos;
}


/*** Matrix::generateSmallWorld_fullRewire ***
 * Generates a small world graph of the Watts-Strogatz (1998) type, optionally with additional self-couplings. *
 * Therefore, at first a regular graph with connections from each vertex to its k nearest neighbors is created. *
 * Next, a specified number of edges is rewired to new, random destinations. *
 * This function establishes rewirings between two new vertices, thus it does not leave one vertex in place, *
 *  It creates a symmetric matrix in the end. *
 * - int k: number of nearest neighbors to which connections shall be established (should be an even number) *
 * - int rew_num: number of edges to be rewired *
 * - double _minS: minimum coupling strength *
 * - double _maxS: maximum coupling strength *
 * - double _minS_self: minimum coupling strength for self-coupling *
 * - double _maxS_self: standard deviation of coupling strength for self-coupling *
 * - ofstream* f [optional]: pointer to a file for information output */
void Matrix::generateSmallWorld_fullRewire(const int k, const int rew_num, const double _minS, const double _maxS, const double _minS_self, const double _maxS_self, ofstream* f)
{
	uniform_real_distribution<double> n_dist_strength(_minS,_maxS); // uniform distribution for coupling strengths
	uniform_real_distribution<double> n_dist_strength_self(_minS_self,_maxS_self); // uniform distribution for coupling strengths for self-coupling

	if (f != NULL)
	{
		*f << "Nearest-neighbor connections = " << k << endl;
		*f << "Rewired connections = " << rew_num << endl;
		*f << "Coupling strength in [" << _minS << ":" << _maxS << "]" << endl;
		*f << "Self-coupling strength in [" << _minS_self << ":" << _maxS_self << "]" << endl << endl;
	}

	reset(); // to obtain an empty matrix

	// Generate regular graph
	for (int i=0; i<N*N; i++)
	{
		int irow = row(i), icol = col(i);
		if (irow == icol) // elements on the main diagonal
		{
			M[cNN(irow,icol)] = n_dist_strength_self(rg_strength); // draw coupling strength from special distribution
		}
		else if (icol > irow && circle_dist(icol,irow,N) <= k/2) // if element is in the upper triangle and not on the main diagonal,
						  			  		   // establish connection if it is a connection between near neighbors
		{
			M[cNN(irow,icol)] = n_dist_strength(rg_strength); // draw coupling strength
		}
		else if (icol < irow) // if element is in the lower triangle and not on the main diagonal,
					    // copy connection strength from upper triangle (because it is an undirected graph)
		{
			M[cNN(irow,icol)] = M[cNN(icol,irow)]; // copy coupling strength
		}
	}

	// Rewiring
	int rewire = rew_num; // number of edges that are to be rewired
	int it = 0;
	while (rewire > 0 && it++ < MAX_ITERATIONS)
	{
		int choice = u_dist_element(rg_conn); // choice for connection to be eliminated (same rand. distribution as for non-zeros can be used)
		int choice_new = u_dist_element(rg_conn); // choice for new connection to be established (same rand. distribution as for non-zeros can be used)

		if (row(choice) != col(choice) // if element 'choice' is not on the main diagonal (i.e., leave the self-couplings untouched)
			&& row(choice_new) != col(choice_new) // and if element 'choice_new' is not on the main diagonal
			&& isNonZero(choice) // and if element 'choice' constitutes an existing connection
			&& isZero(choice_new)) // and if element 'choice_new' does not constitute an existing connection

		{
			M[choice_new] = M[choice]; // create new connection with same coupling strength
			M[symm(choice_new)] = M[choice]; // do the same to the symmetric counterpart
			M[choice] = 0.; // eliminate old connection
			M[symm(choice)] = 0.; // do the same to the symmetric counterpart
			rewire--;
		}
	}

	if (f != NULL)
		*f << "Non-zeros desired / obtained: " << k*N << " / " << nonzeroTest() << endl;
}

/*** Matrix::generateScaleFree ***
 * Generates a scale-free graph of the Barabási-Albert (1999) type, optionally with additional self-couplings. *
 * Therefore, new vertices with a fixed number of randomly drawn edges are continuously added to an existing *
 * population. *
 * - int m0: initial number of (disconnected) vertices (has to be <= N) *
 * - int m: number of edges to be established for each new vertex (has to be <= m0) *
 * - double _minS: minimum coupling strength *
 * - double _maxS: maximum coupling strength *
 * - double _minS_self: minimum coupling strength for self-coupling *
 * - double _maxS_self: standard deviation of coupling strength for self-coupling *
 * - int stabilize_hubs [optional]: indicates if only the largest hubs, and how many of them, shall receive self-stabilization (0 to disable) *
 * - bool is_directed [optional]: indicates if the graph shall be directed (i.e., the matrix will not be symmetric) *
 * - ofstream* f [optional]: pointer to a file for information output */
void Matrix::generateScaleFree(const int m0, const int m, const double _minS, const double _maxS, const double _minS_self, const double _maxS_self, 
                               int stabilize_hubs, bool is_directed, ofstream* f)
{
	uniform_real_distribution<double> u_dist_strength(_minS,_maxS); // uniform distribution for coupling strengths
	uniform_real_distribution<double> u_dist_strength_self(_minS_self,_maxS_self); // uniform distribution for coupling strengths for self-coupling

	int* k = new int[N]; // array that contains the degrees of all nodes
	int k_tot = 0; // total degree (sum of all degrees)

	if (f != NULL)
	{
		*f << "Initial number of vertices = " << m0 << endl;
		*f << "Edges per new vertex = " << m << endl;
		*f << "Coupling strength in [" << _minS << ":" << _maxS << "]" << endl;
		*f << "Self-coupling strength in [" << _minS_self << ":" << _maxS_self << "]" << endl << endl;
	}

	reset(); // to obtain an empty matrix

	// Connect vertices via prefential attachment
	for (int new_vertex=1; new_vertex<=N; new_vertex++) // loop over new vertices
	{
		k[new_vertex-1] = 0; // initially, the degree of each node is zero
		
		if (stabilize_hubs == 0)
			M[cNN(new_vertex,new_vertex)] = u_dist_strength_self(rg_strength); // draw coupling strength for self-coupling from special distribution

		if (new_vertex > m0) // for all added vertices
		{
			int conn = m; // number of edges that are to be established
			int it = 0; // number of iterations

			while (conn > 0 && it++ < MAX_ITERATIONS) // loop to establish m edges to already existing vertices
			{
				int ex_vertex; // chosen existing vertex to which an edge shall be established, self-couplings are left untouched
				int index = 0; // index of the new edge in the adjacency matrix
				uniform_int_distribution<int> u_dist(1, k_tot+new_vertex-1); // uniform distribution between 1 and the total number of degrees plus current number of vertices
												    // - in this range, each node has an interval of k_i+1 outcomes
				int r = u_dist(rg_conn); // drawn random number
				int r_sum = 1; // cumulated position

				for (int i=0; i<new_vertex-1; i++) // loop over existing vertices, finds the vertex whose interval contains r
				{
					int r_sum_next = r_sum + k[i] + 1;

					if (r >= r_sum && r < r_sum_next) // the drawn r lies in the probability interval of vertex i+1
					{
						ex_vertex = i+1;
						index = cNN(ex_vertex, new_vertex);
						break;
					}
					else
						r_sum = r_sum_next;
				}

				if (!isOnMainDiagonal(index) // if this edge is not a self-coupling 
				    && isZero(index)) //and if it was not chosen before (during the current while-loop)
				{
					double coupl_strength = u_dist_strength(rg_strength); // draw coupling strength
					double coupl_strength2;

					if (!is_directed)
						coupl_strength2 = coupl_strength; // set to the same value
					else
						coupl_strength2 = u_dist_strength(rg_strength); // draw another value
					
					if (abs(coupl_strength) > EPSILON && abs(coupl_strength2) > EPSILON) // if the drawn coupling strength(s) unequal zero
					{
						M[index] = coupl_strength; // create new connection with drawn coupling strength
						M[symm(index)] = coupl_strength2; // create symmetric counterpart with the same or a different coupling strength
						k[ex_vertex-1]++; // increase degree of existing node
						k[new_vertex-1]++; // increase degree of new node
						k_tot += 2; // increase the total degree by 2
						conn--;
					}
				}
			}
			if (f != NULL && it >= MAX_ITERATIONS)
				*f << "new_vertex = " << new_vertex << ": only " << m-conn << " edges were established" << endl;
		}
	}
	
	// Information about hubs
	node_and_degree largest_hub = getLargestOutHub();
	vector<node_and_degree> hubs = getOutHubs(stabilize_hubs);
	if (f != NULL)
	{
		*f << "Largest hub: node #" << largest_hub.node << " with degree " << largest_hub.degree << endl;
		
		if (stabilize_hubs > 0)
		{
			*f << "Further considered hubs: ";
			for (int i=1; i<hubs.size(); i++)
			{
				*f << "#" << hubs[i].node << " (" << hubs[i].degree << ")";
				if (i < hubs.size()-1)
					*f << ", ";
			}
			*f << endl;
		}
	}
	
	// Dedicated stabilization of hubs
	if (stabilize_hubs > 0)
	{
		for (int i=0; i<hubs.size(); i++)
		{
			M[cNN(hubs[i].node,hubs[i].node)] = u_dist_strength_self(rg_strength); // draw self-coupling strength for hub
		}
	}

	if (f != NULL)
	{
		*f << "Non-zeros desired / obtained: " << 2*m*(N-m0) << " / " << nonzeroTest() << endl;
	}

	delete[] k;
}


/*** Matrix::generateLeskovec ***
 * Generates a Kronecker graph as described in Leskovec et al., 2005 (1999) type, optionally with additional self-couplings. *
 * The graph of dimension N is created from p equivalent seed graphs of dimension sqrt(N). The desired number of zeros in the seed matrices has *
 * to be specified. *
 * - int p: the Kronecker power that is used
 * - int n0: number of zeros in seed matrices (has to be <= N) *
 * - double _minS: minimum coupling strength *
 * - double _maxS: maximum coupling strength *
 * - double _minS_self: minimum coupling strength for self-coupling *
 * - double _maxS_self: standard deviation of coupling strength for self-coupling *
 * - bool is_directed [optional]: indicates if the graph shall be directed (i.e., the matrix will not be symmetric) *
 * - ofstream* f [optional]: pointer to a file for information output *
 * - return: number of non-zeros supposed to be held by the generated matrix */
int Matrix::generateLeskovec(const int p, const int n0, const double _minS, const double _maxS, const double _minS_self, const double _maxS_self, 
                             bool is_directed, ofstream* f)
{
	int N_seed = round(pow(N, 1./p)); // dimension of seed graphs
	Matrix A(N_seed); // create seed matrix for Kronecker multiplication

	int nonzeros_left = pow(N_seed,2)-n0;
	int zeros_desired;

	uniform_real_distribution<double> u_dist_strength(_minS,_maxS); // uniform distribution for coupling strengths
	uniform_real_distribution<double> u_dist_strength_self(_minS_self,_maxS_self); // uniform distribution for coupling strengths for self-coupling
	uniform_int_distribution<int> u_dist_seed_element(0, (pow(N_seed,2)-N_seed)/2-1); // uniform distribution across half of off-diagonal seed matrix elements

	if (f != NULL)
	{
		*f << "Number of zeros in seed matrix = " << n0 << endl;
		*f << "Coupling strength in [" << _minS << ":" << _maxS << "]" << endl;
		*f << "Self-coupling strength in [" << _minS_self << ":" << _maxS_self << "]" << endl << endl;
	}

	// Generate seed matrix
	A.reset(); // set all entries to zero

	while (nonzeros_left > 0)
	{
		int choice = u_dist_seed_element(rg_conn); // draw element from small distribution (to minimize number of drawings)

		// Determine index of drawn element
		int index;
		int curr_index = 0;
		for (int i=1; i<=N_seed; i++) // go through rows of seed matrix
		{
			int last_index = curr_index;
			curr_index += (N_seed-i); // (N_seed-i) is the number of elements than can be drawn from row i
			if (choice < curr_index) // chosen element is in row i
			{
				index = cNN2(i, i + (choice + 1) - last_index, N_seed);
				break;
			}
		}

		if (A.isZero(index)) // if element is not set yet
		{
			A.set(index, 1.0); // mark element by setting its value to 1.0 (coupling strength is drawn after Kronecker multiplication)
			A.set(symm2(index, N_seed), 1.0); // mark symmetric counterpart by setting its value to 1.0 
			nonzeros_left -= 2;
		}
	}

	if (f != NULL)
	{
		*f << "Non-zeros desired / obtained: " << pow(N_seed,2)-n0 << " / " << A.nonzeroTest() << endl;
		A.saveMatrix(f);
		for (int i=0; i<100; i++)
			*f << "=";
		*f << endl << endl;
	}

	// Do Kronecker multiplication
	zeros_desired = 2*pow(N_seed,2)*n0 - pow(n0,2); // zeros desired for power 2
	kronecker(A, A); // do Kronecker product for power 2 (the result is written to matrix "this")

	for (int pp=3; pp <= p; pp++) // powers higher than 2
	{
		zeros_desired = zeros_desired * pow(N_seed,2) + (pow(N_seed,4) - zeros_desired) * n0;
		kronecker(*this, A); // do Kronecker product (the result is again written to matrix "this")
	}

	for (int i=0; i<N*N; i++) // loop over all elements of final Kronecker matrix to set coupling strengths
	{
		int irow = row(i), icol = col(i);
		if (irow == icol) // elements on the main diagonal
		{
			M[i] = u_dist_strength_self(rg_strength); // draw coupling strength for self-coupling from special distribution
		}
		else if (M[i] > 0.5 && icol > irow) // marked elements in the upper triangle (except the main diagonal)
		{
			M[i] = u_dist_strength(rg_strength); // draw coupling strength for element i

			if (!is_directed)
				M[symm(i)] = M[i]; // use the same coupling strength for the symmetric counterpart of element i
			else
				M[symm(i)] = u_dist_strength(rg_strength); // draw another coupling strength for the symmetric counterpart of element i
		}
	}

	int nonzeros_desired = pow(N,2) - zeros_desired;
	if (f != NULL)
	{
		*f << "Non-zeros desired / obtained: " << nonzeros_desired << " / " << nonzeroTest() << endl;
	}

	return nonzeros_desired;
}

/*** Matrix::trace ***
 * Returns the value of the trace of the matrix *
 * - return: the trace */
double Matrix::trace() const
{
	double tr = 0.;

	for (int i=1; i<=N; i++)
	{
		tr += M[cNN(i,i)];
	}

	return tr;
}


/*** Matrix::add ***
 * Adds to the matrix another matrix *
 * - B: second matrix *
 * - return: the matrix sum */
void Matrix::add(Matrix B)
{
	for (int i=0; i<N*N; i++)
	{
		M[i] += B.get(i);
	}
}

/*** Matrix::multiply ***
 * Multiplies the matrix with another matrix B or with a scalar *
 * - B: second matrix */
void Matrix::multiply(Matrix B)
{
	Matrix ret(N);
	for (int i=1; i<=N; i++)
	{
		for (int j=1; j<=N; j++)
		{
			double v = 0.;
			for (int k=1; k<=N; k++)
			{
				v += M[cNN(i,k)]*B.get(cNN(k,j));
			}
			ret.set(cNN(i,j), v);
		}
	}

	*this = ret;
	//delete ret;
}
/* - s: scalar *
 * - return: the resulting matrix */
void Matrix::multiply(double s)
{
	for (int n=0; n<N*N; n++)
	{
		M[n] *= s;
	}
}


/*** Matrix::mpow ***
 * Takes the n-th power of the matrix *
 * - int n: the power */
void Matrix::mpow(int n)
{
	Matrix ret = *this;
	for (int i=1; i<n; i++)
	{
		ret.multiply(*this);
	}

	*this = ret;
}


/*** Matrix::transpose ***
 * Transposes the matrix */
void Matrix::transpose()
{
	Matrix ret(N);
	for (int n=0; n<N*N; n++)
	{
		ret.set(n, M[symm(n)]);
	}

	*this = ret;
}

/*** Matrix::get ***
 * Gets the specified element of the matrix *
 * - int i: the number of the element (starting from zero) *
 * - return: the value of the element */
double Matrix::get(int i) const
{
	return M[i];
}

/*** Matrix::isNonZero ***
 * Returns if the specified element of the matrix is non-zero *
 * - return: true if non-zero, false if zero */
bool Matrix::isNonZero(int i) const
{
	if (abs(M[i]) >= EPSILON)
		return true;
	else
		return false;
}

/*** Matrix::isZero ***
 * Returns if the specified element of the matrix is (approximately) zero *
 * - return: true if zero, false if non-zero */
bool Matrix::isZero(int i) const
{
	if (abs(M[i]) < EPSILON)
		return true;
	else
		return false;
}

/*** Matrix::set ***
 * Sets the elements of the matrix according to a given array *
 * - double* mat: the consecutive rows of the matrix in an array of length NxN */
void Matrix::set(double* mat)
{
	for (int i=0; i<N*N; i++)
	{
		//M[i / N][i % N] = double(mat[i]);
		M[i] = mat[i];
	}
}

/*** Matrix::set ***
 * Sets an element of the matrix to a given value *
 * - int el: the index of the element *
 * - double val: the new value for the element */
void Matrix::set(int el, double val)
{
	M[el] = val;
}

/*** Matrix::dim ***
 * Gets the dimension of the matrix *
 * - return: the dimension */
int Matrix::dim() const
{
	return N;
}

/*** Matrix::init ***
 * Initializes the object (typically run by constructor) *
 * - _N: the dimension */
void Matrix::init(int _N) //: N(_N), rg_conn(getClockSeed()), rg_strength(getClockSeed()), u_dist_element(0, _N*_N-1), u_dist_conn(0.0, 1.0)
{
	N = _N;
	M = new double [N*N];

	rg_conn = minstd_rand0(getClockSeed());
	rg_strength = minstd_rand0(getClockSeed());
	u_dist_element = uniform_int_distribution<int>(0, N*N-1);
	u_dist_conn = uniform_real_distribution<double>(0.0, 1.0);
}

/*** Matrix::(constructor) ***
 * Sets all parameters, creates matrix *
 * - int _N: the matrix dimension */
Matrix::Matrix(int _N)
{
	init(_N);
	//for (int i=0; i<N; i++)
	//	M[i] = new double [N];
}

/*** Matrix::(constructor) ***
 * Empty constructor, does not create matrix yet */
Matrix::Matrix()
{
	N = 0;
	M = NULL;
}

/*** Constructor ***
 * Copy constructor *
 * - Matrix& newMat: reference to an object of type Matrix */
Matrix::Matrix(const Matrix& newMat)
{
	N = 0;
	M = NULL;
	*this = newMat;
}


/*** Matrix::(destructor) ***
 * Cleans up the allocated memory for arrays */
Matrix::~Matrix()
{
	//for(int i=0; i<N; i++)
	//	delete[] M[i];
	if (M != NULL)
	{
		delete[] M;
		M = NULL;
	}
}

/*** Matrix::(assignment operator) ***
 * Assigns a matrix to another one (is basically a copy constructor) *
 * - Matrix& newMat: reference to an object of type Matrix */
Matrix& Matrix::operator=(const Matrix& newMat)
{
	if (this == &newMat) return *this;

	int N_new = newMat.dim();

	if (N != N_new)
	{
		if (M != NULL)
			delete[] M;
		init(N_new);
	}

	for (int i = 0; i < N*N; i++)
	{
		M[i] = newMat.get(i);
	}

	return *this;
}

/*** Matrix::(plus-assign operator) ***
 * Adds a matrix to another one *
 * - Matrix& Mat2: reference to an object of type Matrix */
Matrix& Matrix::operator+=(const Matrix& Mat2)
{
	int N_new = Mat2.dim();

	if (N != N_new)
	{
		cout << "Error: adding matrices of dimensions " << N << " and " << N_new << " is not possible!" << endl;
		return *this;
	}

	for (int i = 0; i < N*N; i++)
	{
		M[i] += Mat2.get(i);
	}

	return *this;
}

