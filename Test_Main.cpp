/********************************************
 *** Functions to test the implementation ***
 ***    (c) Jannik Luboeinski 2016-2022   ***
 ********************************************/

#include "Main_Stability.hpp"

/*** printStability ***
 * Function that prints textual information about the stability to the given file *
 * - int res: result of previous stability check *
 * -ofstream& f: file to print to */
void printStability(int res, ofstream& f)
{
	if (res == 1)
		f << "stable" << endl;
	else if (res == -1)
		f << "unstable" << endl;
	else
		f << "N/A" << endl;
}

/*** extractInformation ***
 * Function that extracts standard information (char. polynomial, zeros, stability) from a given matrix and prints it to the given file *
 * - int res: result of previous stability check *
 * - ofstream& f: file to print to *
 * - bool rh [optional]: indicates if stability shall also be checked with the Routh-Hurwitz criterion *
 * - bool get_chp_zeros [optional]: indicates if characteristic polynomial and zeros shall be determined *
 * - return: characteristic polynomial as a vector of complex numbers */
vector<complex<double>> extractInformation(Matrix& M, ofstream& f, bool rh = false, bool get_chp_zeros = false)
{
	vector<complex<double>> tmp_coeff;
	int dim = M.dim();
	double* chp_coeff = new double [dim+1]; // coefficients of the characteristic polynomial
	complex<double>* zeros = new complex<double> [dim]; // zeros of the characteristic polynomial

	M.saveMatrix(&f); // save the actual matrix and additional information

	if (get_chp_zeros)
	{
		M.chPolynomial(chp_coeff); // compute the characteristic polynomial
		for (int i=dim; i>0; i--)
			f << chp_coeff[i] << " x^" << i << " + "; // print characteristic polynomial
		f << chp_coeff[0] << endl << endl;
	
		int num_zeros = findZeros(zeros, chp_coeff, dim, &f); // find the zeros, can take very long to compute
		f << "Number of zeros found: " << num_zeros << endl;
		if (num_zeros > 0)
		{
			for (int i=0; i<num_zeros-1; i++)
				f << zeros[i] << "; ";
			f << zeros[num_zeros-1] << endl;
		}
	}

	if (rh)
	{
		printStability(M.isStableRH(&f), f); // Routh-Hurwitz criterion of stability, can take very long to compute
	}

	printStability(M.isStable(&f), f); // numerical assessment of stability

	f << endl;

	for (int i=0; i<=dim; i++) // transfers the coefficients of the characteristic polynomial to a vector of complex numbers
		tmp_coeff.push_back(chp_coeff[i]);

	delete[] chp_coeff;
	delete[] zeros;

	return tmp_coeff;
}

/*** Test ***
 * General purpose test function */
void Test()
{
	const double conn = 0.35; // connectance
	const int dim = 10; // dimension
	
	Matrix A(dim);
	Matrix B(dim);
	Matrix L(100); // for Kronecker-Leskovec matrix
	
	ofstream f("Test_Main_out.txt");
	
	if (!f.is_open())
	{
		cout << "Failed to open file!" << endl;
		return;
	}
	
	////////////////////////////
	// Pre-determined matrix
	f << "----------------------------------------------------" << endl;
	f << "Pre-determined matrix:" << endl << endl;
	//double mat [] = {5., 45., -66., -10.}; // 2 x 2, has complex eigenvalues
	//double mat [] = {2., 1., 0., 0., 6., 0., 0., 0., 12.}; // 3 x 3
	//double mat [] = {-0.177,0.978,-0.920, 0.000,0.000,-0.282,-0.999,0.000,0.000,0.000,-0.765,0.790,0.000,-0.639,0.467,-0.413}; // 4 x 4
	//double mat [] = {-0.119,0.000, 0.000, 0.000,-0.093,0.000,-0.280, 0.000, 0.310,-0.004,0.570, 0.000,-0.900, 0.000, 0.000,0.000, 0.000,-0.412,-0.867,0.000,0.000,0.361,0.452, 0.000, -0.699}; // 5 x 5
	//double mat [] = {-0.455, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,-0.212, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,-0.194, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,-0.139, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,-0.603, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,-0.342, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,-0.575, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,-0.257, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,-0.873, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,-0.133}; // 10 x 10
	double mat[] = {-0.914,-0.407, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,-0.158,-0.407,-0.938,-0.071, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,-0.071,-0.651,-0.154, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,-0.154,-0.635, 0.032, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.032,-0.724,-0.284, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,-0.284,-0.985,-0.756, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,-0.756,-0.979,-0.115, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,-0.115,-0.135,-0.042, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,-0.042,-0.600, 0.001,-0.158, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.001,-0.837};
	A.set(mat);
	vector<complex<double>> tmp_coeff = extractInformation(A, f, true, true);
	vector<complex<double>> tmp_coeff2 = tmp_coeff;

	derivePolynomial(&tmp_coeff[0], dim);
	for (int i=dim; i>0; i--) // print derived polynomial
		f << tmp_coeff[i] << " x^" << i << " + ";
	f << tmp_coeff[0] << endl << endl;
	
	if (polynomialLongDivision(&tmp_coeff2[0], -0.129018, dim))
	{
		for (int i=dim-1; i>0; i--)
			f << tmp_coeff2[i] << " x^" << i << " + ";
		f << tmp_coeff2[0] << endl << endl;
	}

	////////////////////////////
	// Small-world matrices, undirected
	f << "----------------------------------------------------" << endl;
	f << "Small-world matrices, undirected:" << endl << endl;
	B.generateSmallWorld(4, 0.05, -1., +1, -2.8, -1.9, true, &f);
	extractInformation(B, f, true, true);

	B.generateSmallWorld(4, 0.05, -1., +1, -2.8, -1.9, true, &f);
	extractInformation(B, f, true, true);

	////////////////////////////
	// Small-world matrices, bidirectional
	f << "----------------------------------------------------" << endl;
	f << "Small-world matrices, bidirectional:" << endl << endl;
	B.generateSmallWorld(4, 0.05, -1., +1, -2.8, -1.9, false, &f);
	extractInformation(B, f, true, true);

	B.generateSmallWorld(4, 0.05, -1., +1, -2.8, -1.9, false, &f);
	extractInformation(B, f, true, true);

	////////////////////////////
	// Scale-free matrices, undirected
	f << "----------------------------------------------------" << endl;
	f << "Scale-free matrices, undirected:" << endl << endl;
	B.generateScaleFree(0, 2, -1., +1, -2.8, -1.9, 0, true, &f);
	extractInformation(B, f, true, true);

	B.generateScaleFree(0, 2, -1., +1, -2.8, -1.9, 0, true, &f);
	extractInformation(B, f, true, true);

	////////////////////////////
	// Scale-free matrices, bidirectional
	f << "----------------------------------------------------" << endl;
	f << "Scale-free matrices, bidirectional:" << endl << endl;
	B.generateScaleFree(0, 2, -1., +1, -2.8, -1.9, 0, false, &f);
	extractInformation(B, f, true, true);

	B.generateScaleFree(0, 2, -1., +1, -2.8, -1.9, 0, false, &f);
	extractInformation(B, f, true, true);

	//////////////////////////////////////////////////////
	// Uniform random matrix (exact number of non-zeros)
	f << "----------------------------------------------------" << endl;
	f << "Uniform random matrix (exact number of non-zeros):" << endl << endl;
	B.generateUniformRandom(conn, -1., 1., -1., -0.1, false, &f);
	extractInformation(B, f, true);

	////////////////////////////
	// Uniform random matrix
	f << "----------------------------------------------------" << endl;
	f << "Uniform random matrix:" << endl << endl;
	B.generateUniformRandomNonExact(conn, -1., 1., -1., -0.1, false, &f);
	extractInformation(B, f, true);

	////////////////////////////
	// Kronecker-Leskovec matrix, bidirectional
	f << "----------------------------------------------------" << endl;
	f << "Kronecker-Leskovec matrix, bidirectional:" << endl << endl;
	L.generateLeskovec(2, 82, -1., +1, -2.8, -1.9, false, &f);
	extractInformation(L, f);

	////////////////////////////
	// Gaussian random matrix
	f << "----------------------------------------------------" << endl;
	f << "Gaussian random matrix:" << endl << endl;
	B.generateGaussRandom(conn, 3, 2, -5, 3, &f);
	extractInformation(B, f);

	f.close();	
}

int main()
{

	Test();

}
