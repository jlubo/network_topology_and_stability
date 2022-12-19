/********************************************
 *** Useful generic functions and classes ***
 ***    (c) Jannik Luboeinski 2015-2022   ***
 ********************************************/

#include "Tools.hpp"

/*** node_and_degree::(constructor) ***
 * _node: number of the node *
 * _degree: the degree of the node *
 * return: an instance of the class */
node_and_degree::node_and_degree(int _node, int _degree)
{
	node = _node;
	degree = _degree;
}

/*** concat ***
 * Concatenates two character arrays/strings and returns a string *
 * - str1: the first character array
 * - str2: the second character array
 * - return: the concatenated character array in a new string variable */
string concat(const char *str1, const char *str2)
{
	string s = string(str1) + string(str2);
	return s;
}
string concat(string str1, const char *str2)
{
	string s = str1 + string(str2);
	return s;
}
string concat(const char *str1, string str2)
{
	string s = string(str1) + str2;
	return s;
}

/*** sortvec ***
 * Function to sort a vector of node_and_degree objects in descending degree order *
 * (to be used by std::sort()). *
 * - a: the first element to be compared
 * - b: the second element to be compared
 * - return: true if the degree of element a is greater than that of b */
bool sortvec(const node_and_degree &a,
              const node_and_degree &b)
{
    return (a.degree > b.degree);
}

/*** dtos ***
 * Converts a double variable to a string, rounded to the given number of internal decimal places, *
 * using the C method sprintf *
 * - num: the double floating point number to be converted *
 * - n: the number of internal decimal places (has to be greater than or equal to 0) *
 * - dont_allow_zero [optional]: specifies if zeros shall be replaced by "nan" values *
 * - return: the string containing the rounded floating point number */
string dtos(double num, int n, bool dont_allow_zero)
{
	string s; // the string to be returned
	char *buf; // the buffer the converted number is written to
	char *format; // the format string for sprintf
	int prepoint; // number of places before point, has to be at least 1
	int dn = (n > 0 ? dn = (int) floor(log10(n)) + 1 : dn = 1); // decimal places of the variable n (needed to generate the format string)
	char nanstring [] = {"nan"};

	if (!std::isnan(num))
	{
		if (num <= 0.0 || floor(log10(num)) < 0.0)
			prepoint = 1;
		else
			prepoint = (int) floor(log10(num)) + 1;

		buf = new char[prepoint + n + 2]; // decimal places plus internal decimal places of num plus 2 for point and zero-terminator
		format = new char[dn + 4]; // decimal places of n plus 4 for the characters "%.f\0"

		sprintf(format, "%%.%df", n); // create format string
		sprintf(buf, (const char*) format, num); // create character array from double
		s.assign(buf); // create buffer

		if (num == 0. && dont_allow_zero)
			s.assign(nanstring);

		delete[] buf;
		delete[] format;
	}
	else
		s.assign(nanstring);

	return s;
}

/*** dateStr ***
 * Adds a string containing a time stamp of 18 char's, *
 * optionally added to the start of another string *
 * - str [optional]: a string to which the time stamp will be added  *
 * - fixdate [optional]: true if this time shall be fixed, false by default *
 * - return: the date string plus the specified string, if stated */
static string datestr = "";
string dateStr(string str, bool fixdate)
{
	const int len = 19; // time stamp length: 18 characters (plus terminating \0)

	if (fixdate == true) // if this date is to be fixed, update date
	{
		char buf[len];
		string new_datestr;
		do
		{
			time_t t = time(NULL);
			struct tm *now  = localtime(&t);
			strftime(buf, len, "%y-%m-%d_%H-%M-%S", now);
			new_datestr.assign(buf);
		} while (new_datestr == datestr);
		datestr = new_datestr;
	}

	return datestr + str;
}

/*** writeLineToFile ***
 * Writes a line of text to a file, if provided *
 * - string line: the text *
 * - ofstream* f: pointer to file for output */
void writeLineToFile(string line, ofstream *f)
{
	if (f == NULL)
		return;

	*f << line << endl;
}

/*** sgn ***
 * Returns the sign of a double number *
 * - num: the number of which the sign shall be determined
 * - return: -1 for negative sign, +1 for positive sign and zero */
int sgn (double num)
{
	if (num < -0.0)
		return -1;
	else
		return +1;
}

/*** timeMeasure ***
 * Starts or stops a time measurement (accuracy is one second) *
 * - start: boolean which is true when intending to start and false when intending to stop measurement
 * - return: the seconds elapsed */
int timeMeasure (bool start)
{
	static int start_time;
	if (start == true)
	{
		start_time = time(NULL);
		return 0;
	}
	else
	{
		return (time(NULL) - start_time);
	}
}

/*** showChDirErrMessage ***
 * Prints the last errno message out of the possible chdir() errors via cout */
void showChDirErrMessage()
{
	cout << "Error while changing directory: ";
	switch(errno)
	{
		case EACCES:
			cout << "EACCES" << endl;
			break;
		case ENAMETOOLONG:
			cout << "ENAMETOOLONG" << endl;
			break;
		case ENOENT:
			cout << "ENOENT" << endl;
			break;
		case ENOTDIR:
			cout << "ENOTDIR" << endl;
			break;
		case ELOOP:
			cout << "ELOOP" << endl;
			break;
	}
}

/*** getClockSeed ***
 * Gets a random generator seed from the computer clock, guarantees not to return *
 * the same seed in two subsequent calls (very important!) *
 * - return: clock counts since epoch of the clock */
static unsigned int last_seed = 0;
unsigned int getClockSeed()
{
	int seed;

	while ( (seed = system_clock::now().time_since_epoch().count()) == last_seed ) {}
	last_seed = seed;

	return seed;
}

/*** copyFile ***
 * Copies the content of a file from a source file to a target file *
 * - src: file name of the source file
 * - dest: file name of the target file
 * - return: false if one of the files could not be opened, true if succesful */
bool copyFile(string src, string dest)
{
	ifstream ssrc(src, ios::binary);
	ofstream sdest(dest, ios::binary);

	if (!ssrc.is_open() || !sdest.is_open())
		return false;

	sdest << ssrc.rdbuf();

	ssrc.close();
	sdest.close();

	return true;
}


/*** circle_dist ***
 * Calculates the distance in units of integer numbers between two elements in a circle *
 * of numbered elements, starting from 1 *
 * - i: the number of the first element
 * - j: the number of the second element
 * - len: the length of the circle (the highest element number in the circle)
 * - return: the distance between the elements */
int circle_dist(int i, int j, int len)
{
	int first_dist = abs(i - j); // the simple distance on the linear line
	int second_dist = len - first_dist; // the complemetary distance

	return (first_dist < second_dist ? first_dist : second_dist);
}


/*** derivePolynomial ***
 * Takes the derivative of a polynomial of given order *
 * - pv: consumes the coefficients of the original polynomial and returns the coefficients of the derivative (coefficient for x^i with array index i) *
 * - n: the order of the polynomial to be derived */
void derivePolynomial(complex<double>* pv, int n)
{
	for (int i=0; i<n; i++)
	{
		pv[i] = (complex<double>)(i+1) * pv[i+1];
	}
	pv[n] = 0.;
}
void derivePolynomial(double* pv, int n)
{
	for (int i=0; i<n; i++)
	{
		pv[i] = (double)(i+1) * pv[i+1];
	}
	pv[n] = 0.;
}

/*** polynomialLongDivision ***
 * Divides a complex polynomial by (x - div), eliminating one zero *
 * - pv: passes the polynomial coefficients, returns the coefficients of the resulting polynomial (coefficient for x^i with array index i) *
 * - div: the zero of the original polynomial that should be eliminated *
 * - n: the order of the polynomial to be divided *
 * - err_output [optional]: pointer to a file for error output *
 * - return: true if successful (no remarkable residue remains), false if not */
bool polynomialLongDivision(complex<double>* pv, complex<double> div, int n, ofstream* err_output)
{
	complex<double>* ret;

	if (n <= 0)
		return false;
	ret = new complex<double>[n];

	ret[n-1] = pv[n];

	for (int i=n-1; i>0; i--) // perform division
	{
		ret[i-1] = pv[i] - ((-1.) * div * ret[i]);
	}

	if (abs((pv[0] - ((-1.) * div * ret[0]))) > EPSILON2) // check - no residue remains?
	{
		if (err_output)
		{
			int prec = err_output->precision();
			*err_output << "Residue: " << (pv[0] - ((-1.) * div * ret[0]));
			err_output->precision(log10(10./EPSILON2)); // set precision of the order of magnitude of EPSILON2
			*err_output << ", abs = " << abs((pv[0] - ((-1.) * div * ret[0]))) << endl;
			err_output->precision(prec); // set back to previous precision
		}
		delete[] ret;

		return false;
	}

      pv[n] = 0.;
	memcpy(pv, ret, n*sizeof(complex<double>)); // copy coefficients of resulting polynomial

	delete[] ret;

	return true;
}

/*** findZeros ***
 * Computes the zeros of a given real polynomial employing Horner's/Newton's method *
 * - zeros: returns the zeros in an array of complex numbers (has to have length n) *
 * - pv: passes the polynomial coefficients (coefficient for x^i with array index i) *
 * - n: the order of the polynomial *
 * - err_output [optional]: pointer to a file for error output *
 * - return: the number of found zeros, or -1 if not successful */
int findZeros(complex<double>* zeros, double* pv, int n, ofstream* err_output)
{
	int ret;
	complex<double> x = complex<double>(0.,1.0); // current value for the zero, use 0+i as first guess for the starting point
	complex<double> *pv2, *pv3; // the polynomial (pv2) and its derivative (pv3)

	if (n <= 0)
		return -1;

	pv2 = new complex<double>[n+1];
	pv3 = new complex<double>[n+1];

	for (int i=0; i<n+1; i++) // transfer real coefficients into complex numbers
		pv2[i] = complex<double>(pv[i], 0.);

	// Horner's method with complex coefficients to find all zeros (possible alternatives: Bairstow's method, Jenkins-Traub RPOLY)
	for (int j=n; j>2; j--) // decrease the order to find all zeros
	{
		// Newton's method
		int it; // number of iterations
		complex<double> x_initial = x; // use previously obtained zero as starting point for the next one
		complex<double> sol1, sol2; // polynomial at x, derivative of the polynomial at x


		memcpy(pv3, pv2, (n+1)*sizeof(complex<double>));
		derivePolynomial(pv3, j);
#ifdef ENFORCE_CONVERGENCE
		int it2 = 0;
		complex<double> enforce_step = complex<double>(1.,0.); // step size for enforcing convergence
		bool negative_enforcement = false; // specifies if enforcing convergence by decreasing the number has already been done

		do // loop over directions of variation of x_initial
		{
			if (it2 > 0) // positive enforcement has already been done and was not successful
			{
				negative_enforcement = true;
				x_initial = x_initial - (double) MAX_ITERATIONS * complex<double>(1.,1.); // return to initial value
				enforce_step = complex<double>(-1.,0.); // negative direction
				it2 = 0;
			}

			do // loop over x_initial variations
			{
#endif
				x = x_initial;
				it = 0;

				do
				{
					sol1 = 0.;
					sol2 = 0.;
					for (int i=0; i<=j; i++) // compute values of the polynomial and of its derivative at x
					{
						sol1 += pv2[i] * pow(x,i);

						if (i < j)
							sol2 += pv3[i] * pow(x,i);
					}

#ifdef ENFORCE_CONVERGENCE
					if (sol2 != 0.)
#endif
						x = x - (sol1 / sol2);
#ifdef ENFORCE_CONVERGENCE
					else  // derivative vanishes - modify x
					{
						x += enforce_step;
						if (err_output) *err_output << "Zero #" << n-j+1 << ": conv. enforced, x = " << x << endl;
					}
#endif
					//if (err_output) *err_output << "sol1 = " << sol1 << "; sol2 = " << sol2 << endl;
				} while (abs(sol1) > EPSILON && it++ < MAX_ITERATIONS);

#ifdef ENFORCE_CONVERGENCE
				x_initial += enforce_step; // increase x to escape infinite cycle

			} while (it > MAX_ITERATIONS && it2++ < MAX_ITERATIONS); // do while process keeps being caught in an infinite cycle

		} while (it2 > MAX_ITERATIONS && !negative_enforcement); // do while process is assumed to be caught in an infinite cycle

		if (it2 > 0)
		{
			if (err_output)
				*err_output << "Zero #" << n-j+1 << ": conv. enforced, x_initial = " << x_initial-enforce_step << endl;
		}
#endif

		if (!polynomialLongDivision(pv2, x, j, err_output)) // perform polynomial long division using the found zero to decrease the order
		{
			if (err_output) 
				*err_output << "Iterations: it = " << it 
#ifdef ENFORCE_CONVERGENCE
				            << ", it2 = " << it2
#endif
				            << endl;
			return n-j; // if it fails, this zero is not saved
		}

		zeros[n-j] = x;
		//if (err_output) *err_output << "Zero #" << n-j+1 << ": " << x << endl;
	}

	// Apply "abc" formula for the last two zeros
	zeros[n-2] = (-pv2[1]-sqrt(pow(pv2[1], 2) - 4.*pv2[2]*pv2[0])) / (2.*pv2[2]);
	zeros[n-1] = (-pv2[1]+sqrt(pow(pv2[1], 2) - 4.*pv2[2]*pv2[0])) / (2.*pv2[2]);

	return n;
}

