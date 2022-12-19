/**************************************************
 *** Simulations for assessing the stability of ***
 ***      networks of different topology        ***
 ***     (c) Jannik Luboeinski 2016-2022        ***
 **************************************************/

#include "Main_Stability.hpp"

/*** SmallWorld ***
 * Evaluates the stability of small world graphs similar to Watts & Strogatz, 1998 *
 * - N: number of system constituents *
 * - minS_self: minimal self-coupling strength *
 * - maxS_self: maximal self-coupling strength *
 * - sample_number: number of samples per parameter setting *
 * - _k: number of nearest-neighbor connections for initial regular graph *
 * - p_rew: fraction of connections to be rewired (if <0, a parameter scan will be performed) *
 * - int graph_struct: indicates if the graph shall be undirected (0) or asymmetrically bidirectional (1) */
void SmallWorld(int N, double minS_self, double maxS_self, int sample_number, int _k, double p_rew, int graph_struct)
{
	const double minS = -1., maxS = 1.; // minimal and maximal coupling strength between different parts of the system
	
	int k_min; // minimum number of nearest-neighbor connections
	int k_max; // maximum number of nearest-neighbor connections
	const int k_step = 2; // step size for variation of nearest-neighbor couplings (has to be even!!!)
	
	int rew_num_min; // minimum number of edges to be rewired
	int rew_num_max; // maximum number of edges to be rewired
	const int rew_num_step = 2; // step size for variation of rewiring number
	
	Matrix M (N);
	int* degree_dist = new int[N+1]; // counters for all degrees of vertices in the graph

#ifdef DETAILED_OUTPUT
	ofstream f(dateStr("_SmallWorld_N" + dtos((double)N,0) + ".txt")); // file containing parameter information
#endif
	ofstream data(dateStr("_SmallWorld.txt")); // file containing final stability data
	ofstream metadata("../stability_over_self_coupling.txt", ofstream::app); // file containing final stability data of several settings (but in this run only one setting is written to it)
	
	if (
#ifdef DETAILED_OUTPUT
		!f.is_open() ||
#endif
		!data.is_open() || !metadata.is_open())
	{
		cout << "Failed to open file!" << endl;
		return;
	}

	if (_k > 0) // if number of nearest neighbors was specified, use it
	{
		k_min = _k;
		k_max = k_min;
	}
	else // else, loop over different nearest neighbor numbers
	{
		k_min = k_step; // minimum number of nearest-neighbor connections
		k_max = N-k_step; // maximum number of nearest-neighbor connections
	}
	
	if (p_rew >= 0.) // if rewiring fraction was specified, use it
	{
		rew_num_max = (int) round(p_rew*k_min*N/2);
		rew_num_min = rew_num_max;
	}
	else // else, loop over different rewiring numbers
	{
		rew_num_min = 0; // minimum number of edges to be rewired
		rew_num_max = k_min*N/2; // maximum number of edges to be rewired
	}
	
	for (int rew_num = rew_num_min; rew_num<=rew_num_max; rew_num+=rew_num_step) // loop over rewiring numbers
	{
		if (rew_num != rew_num_min)
			printf("\n");
		
		for (int k = k_min; k<=k_max; k+=k_step) // loop over connectances (here introduced through the even number k of neighbor connections)
		{
			int stable = 0;
			int samples = sample_number;
		
			printf("\rrew_num = %d, k = %d", rew_num, k);
			fflush(stdout);
			
			memset(degree_dist, 0, (N+1)*sizeof(int)); // set counters for all degrees to zero

			for (int j=0; j<sample_number; j++)
			{
				int ret;
#ifdef DETAILED_OUTPUT	
				M.generateSmallWorld(k, rew_num, minS, maxS, minS_self, maxS_self, (bool) graph_struct, &f);

				if ((ret = M.isStable(&f)) == 1) 
				{
					stable++;
					f << "stable" << endl;
				}
				else if (ret == -1)
					f << "unstable" << endl;
				else
				{
					samples--; // this sample cannot be used
					f << "N/A" << endl;
				}
#ifdef EV_OUTPUT
				complex<double>* ev = new complex<double>[N]; // array of eigenvalues
				if (M.eigenvalues(ev, EV_METHOD_GSL, &f))
				{
					f << "EV = {";
					for (int n_ev=0;n_ev<N;n_ev++)
					{
						f << ev[n_ev];
						if (n_ev < N-1)
							f << ", ";
					}
					f << "}" << endl;
				}
				Matrix MT = M;
				Matrix B = M;
				Matrix BSQ;

				MT.transpose();
				B.add(MT);
				B.multiply(0.5);
				BSQ = B;
				BSQ.mpow(2);

				double mb = B.trace()/N;
				double sb = BSQ.trace()/N - pow(mb/N, 2);
				
				f << "EV mean = " << mb << endl // mean of eigenvalues
				  << "EV std. dev. = " << sb << endl // standard deviation of eigenvalues
				  << "Lower bound = " << mb - sb * sqrt(N-1) << endl // lower bound
				  << "Upper bound = " << mb + sb * sqrt(N-1) << endl << endl; // upper bound
#endif
#ifdef MAT_OUTPUT
				M.saveMatrix(&f);
#endif
				for (int i=0; i<100; i++)
					f << "=";
				f << endl << endl;
#else
				M.generateSmallWorld(k, rew_num, minS, maxS, minS_self, maxS_self, (bool) graph_struct);
				if ((ret = M.isStable()) == 1) 
					stable++;
				else if (ret == 0)
					samples--; // this sample cannot be used
#endif
				M.sumDegreeDist(degree_dist); // adds data from current matrix to the degree distribution
			}
			data.precision(1);
			metadata.precision(1);
			if (rew_num_min != rew_num_max)
				data << fixed << rew_num << "\t\t" << 100. * k/double((N-1)) << "\t\t" << stable << "\t\t" << samples << endl; // k / (N-1) = (k*N) / (N*N - N) is the connectance
			else
			{
				data << fixed << 100. * k/double((N-1)) << "\t\t" << stable << "\t\t" << samples << endl;
				metadata << fixed << minS_self << "\t\t" << maxS_self << "\t\t" << stable << "\t\t" << samples << endl;
				
				string fstr = concat("SmallWorld_dist_k",dtos((double)k,0) + string("_")+dtos((double)rew_num,0)); // string describing the data
				plot_degree_dist(fstr.c_str(), "Degree distribution", degree_dist, N, N*sample_number);
			}
		}
	}
	cout << endl;
#ifdef DETAILED_OUTPUT
	f.close();
#endif
	data.close();
	metadata.close();
	delete[] degree_dist;
	
	if (k_min != k_max)
		plot_over_connectance("SmallWorld", "Results for small world network", N);
	if (rew_num_min != rew_num_max)
		plot_over_rewirings("SmallWorld", "Results for small world network", N, rew_num_min, rew_num_max);
}

/*** ScaleFree ***
 * Evaluates the stability of scale-free graphs similar to Barabási & Albert, 1999 *
 * - N: number of system constituents *
 * - minS_self: minimal self-coupling strength *
 * - maxS_self: maximal self-coupling strength *
 * - sample_number: number of samples per parameter setting *
 * - m0: size of initial vertex set *
 * - m: number of edges for each new vertex *
 * - stabilize_hubs: indicates if only the largest hubs, and how many of them, shall receive self-stabilization (0 to disable) *
 * - int graph_struct: indicates if the graph shall be undirected (0) or asymmetrically bidirectional (1) */
void ScaleFree(int N, double minS_self, double maxS_self, int sample_number, int m0, int m, int stabilize_hubs, int graph_struct)
{
	const double minS = -1., maxS = 1.; // minimal and maximal coupling strength between different parts of the system
	
	Matrix M (N);
#ifdef DETAILED_OUTPUT
	ofstream f(dateStr("_ScaleFree_N" + dtos((double)N,0) + ".txt")); // file containing parameter information
#endif
	ofstream data(dateStr("_ScaleFree.txt")); // file containing final stability data
	ofstream metadata("../stability_over_self_coupling.txt", ofstream::app); // file containing final stability data of several settings (but in this run only one setting is written to it)
	
	int* degree_dist = new int[N+1]; // counters for all degrees of vertices in the graph

	if (
#ifdef DETAILED_OUTPUT
		!f.is_open() ||
#endif
		!data.is_open() || !metadata.is_open())
	{
		cout << "Failed to open file!" << endl;
		return;
	}
	
	int stable = 0; // total number of stable systems
	int samples = sample_number; // total number of valid samples
	
	memset(degree_dist, 0, (N+1)*sizeof(int)); // set counters for all degrees to zero
	
	printf("\rm = %d, m0 = %d", m, m0);
	fflush(stdout);

	for (int j=0; j<sample_number; j++)
	{
		int ret;
#ifdef DETAILED_OUTPUT	
		M.generateScaleFree(m0, m, minS, maxS, minS_self, maxS_self, stabilize_hubs, (bool) graph_struct, &f);

		if ((ret = M.isStable(&f)) == 1) 
		{
			stable++;
			f << "stable" << endl;
		}
		else if (ret == -1)
			f << "unstable" << endl;
		else
		{
			samples--; // this sample cannot be used
			f << "N/A" << endl;
		}
#ifdef EV_OUTPUT
		complex<double>* ev = new complex<double>[N]; // array of eigenvalues
		if (M.eigenvalues(ev, EV_METHOD_GSL, &f))
		{
			f << "EV = {";
			for (int n_ev=0;n_ev<N;n_ev++)
			{
				f << ev[n_ev];
				if (n_ev < N-1)
					f << ", ";
			}
			f << "}" << endl;
		}
		Matrix MT = M;
		Matrix B = M;
		Matrix BSQ;

		MT.transpose();
		B.add(MT);
		B.multiply(0.5);
		BSQ = B;
		BSQ.mpow(2);

		double mb = B.trace()/N;
		double sb = BSQ.trace()/N - pow(mb/N, 2);
		
		f << "EV mean = " << mb << endl // mean of eigenvalues
		  << "EV std. dev. = " << sb << endl // standard deviation of eigenvalues
		  << "Lower bound = " << mb - sb * sqrt(N-1) << endl // lower bound
		  << "Upper bound = " << mb + sb * sqrt(N-1) << endl << endl; // upper bound	
#endif
#ifdef MAT_OUTPUT
		M.saveMatrix(&f);
#endif
		for (int i=0; i<100; i++)
			f << "=";
		f << endl << endl;
#else
		M.generateScaleFree(m0, m, minS, maxS, minS_self, maxS_self, stabilize_hubs, (bool) graph_struct);
		if ((ret = M.isStable()) == 1) 
			stable++;
		else if (ret == 0)
			samples--; // this sample cannot be used
#endif
		M.sumDegreeDist(degree_dist); // adds data from current matrix to the degree distribution
	}
	data.precision(1);
	metadata.precision(1);

	data << fixed << (2*m*(N-m0)) / double(N*N - N) << "\t\t" << stable << "\t\t" << samples << endl;
	metadata << fixed << minS_self << "\t\t" << maxS_self << "\t\t" << stable << "\t\t" << samples << endl;
	
	string fstr = concat("ScaleFree_dist_m",dtos((double)m,0) + string("_")+dtos((double)m0,0)); // string describing the data
	plot_degree_dist(fstr.c_str(), "Degree distribution", degree_dist, N, N*sample_number);

	cout << endl;
	
#ifdef DETAILED_OUTPUT
	f.close();
#endif
	data.close();
	metadata.close();
	
	delete[] degree_dist;
}


/*** Leskovec ***
 * Evaluates the stability of Kronecker graphs as described in Leskovec et al., 2005 *
 * - N_seed: dimension of the seed matrix *
 * - p: Kronecker power that is used (leading to a matrix of dimension N^p) *
 * - n0_set: fixed value of number of off-diagonal zeros in seed matrix *
 * - minS_self: minimal self-coupling strength *
 * - maxS_self: maximal self-coupling strength *
 * - sample_number: number of samples per parameter setting *
 * - int graph_struct: indicates if the graph shall be undirected (0) or asymmetrically bidirectional (1) */
void Leskovec(int N_seed, int p, int n0_set, double minS_self, double maxS_self, int sample_number, int graph_struct)
{
	const double minS = -1., maxS = 1.; // minimal and maximal coupling strength between different parts of the system

	const int n0_min = n0_set; // minimum number of zeros in seed matrix
	const int n0_max = n0_set; // maximum number of zeros in seed matrix
	const int n0_step = 1; // step size for variation of number of zeros in seed matrix
	int nonzeros_final; // number of zeros in final matrix
	int N = pow(N_seed,p); // dimension of the Kronecker matrix

	Matrix M(N); // the Kronecker matrix
	int* degree_dist = new int[N+1]; // counters for all degrees of vertices in the graph

#ifdef DETAILED_OUTPUT
	ofstream f(dateStr("_Leskovec_N" + dtos((double)N,0) + "_p" + dtos((double)p,0) + ".txt")); // file containing parameter information
#endif
	ofstream data(dateStr("_Leskovec.txt")); // file containing final stability data
	ofstream metadata("../stability_over_self_coupling.txt", ofstream::app); // file containing final stability data of several settings (but in this run only one setting is written to it)

	if (
#ifdef DETAILED_OUTPUT
		!f.is_open() ||
#endif
		!data.is_open() || !metadata.is_open())
	{
		cout << "Failed to open file!" << endl;
		return;
	}	

	for (int n0 = n0_min; n0<=n0_max; n0+=n0_step) // loop over number of zeros in seed matrix
	{
		int stable = 0;
		int samples = sample_number;

		if (n0 != n0_min)
			printf("\n");
			
		printf("\rn0 = %d", n0);
		fflush(stdout);

		memset(degree_dist, 0, (N+1)*sizeof(int)); // set counters for all degrees to zero

		for (int j=0; j<sample_number; j++)
		{
			int ret;
#ifdef DETAILED_OUTPUT	
			nonzeros_final = M.generateLeskovec(p, n0, minS, maxS, minS_self, maxS_self, (bool) graph_struct, &f);

			if ((ret = M.isStable(&f)) == 1) 
			{
				stable++;
				f << "stable" << endl;
			}
			else if (ret == -1)
				f << "unstable" << endl;
			else
			{
				samples--; // this sample cannot be used
				f << "N/A" << endl;
			}
#ifdef EV_OUTPUT
			complex<double>* ev = new complex<double>[N]; // array of eigenvalues
			if (M.eigenvalues(ev, EV_METHOD_GSL, &f))
			{
				f << "EV = {";
				for (int n_ev=0;n_ev<N;n_ev++)
				{
					f << ev[n_ev];
					if (n_ev < N-1)
						f << ", ";
				}
				f << "}" << endl;
			}
			Matrix MT = M;
			Matrix B = M;
			Matrix BSQ;

			MT.transpose();
			B.add(MT);
			B.multiply(0.5);
			BSQ = B;
			BSQ.mpow(2);

			double mb = B.trace()/N;
			double sb = BSQ.trace()/N - pow(mb/N, 2);
			
			f << "EV mean = " << mb << endl // mean of eigenvalues
			  << "EV std. dev. = " << sb << endl // standard deviation of eigenvalues
			  << "Lower bound = " << mb - sb * sqrt(N-1) << endl // lower bound
			  << "Upper bound = " << mb + sb * sqrt(N-1) << endl << endl; // upper bound	
#endif
#ifdef MAT_OUTPUT
			M.saveMatrix(&f);
#endif
			for (int i=0; i<100; i++)
				f << "=";
			f << endl << endl;
#else
			nonzeros_final = M.generateLeskovec(p, n0, minS, maxS, minS_self, maxS_self, (bool) graph_struct);
			if ((ret = M.isStable()) == 1) 
				stable++;
			else if (ret == 0)
				samples--; // this sample cannot be used
#endif
			M.sumDegreeDist(degree_dist); // adds data from current matrix to the degree distribution
		}
		data.precision(1);
		metadata.precision(1);
		if (n0_min != n0_max)
			data << fixed << n0 << "\t\t" << nonzeros_final / double(N*N - N) << "\t\t" << stable << "\t\t" << samples << endl; // nonzeros_final / (N*N - N) is the connectance
		else
		{
			data << fixed << nonzeros_final / double(N*N - N) << "\t\t" << stable << "\t\t" << samples << endl;
			metadata << fixed << minS_self << "\t\t" << maxS_self << "\t\t" << stable << "\t\t" << samples << endl;

			string fstr = string("Leskovec_dist_") + dtos((double)n0_set,0); // string describing the data
			plot_degree_dist(fstr.c_str(), "Degree distribution", degree_dist, N, N*sample_number);
		}
	}
	cout << endl;
	
#ifdef DETAILED_OUTPUT
	f.close();
#endif
	data.close();
	metadata.close();
	delete[] degree_dist;
	
	if (n0_min != n0_max)
		plot_over_connectance("Leskovec", "Results for Leskovec network of power " + dtos(p,0), N);
}


/*** RandomGraph ***
 * Evaluates the stability of Erdős-Rényi(-Gilbert) random graphs of arbitrary size *
 * - N: number of system constituents *
 * - minS_self: minimal self-coupling strength *
 * - maxS_self: maximal self-coupling strength *
 * - sample_number: number of samples per parameter setting *
 * - int graph_struct: indicates if the graph shall be undirected (0), asymmetrically bidirectional (1), or freely directed (2) */
void RandomGraph(int N, double minS_self, double maxS_self, int sample_number, int graph_struct)
{
	const int conn_step = 5; // step size for variation of connectances
	const double minS = -1., maxS = 1.; // minimal and maximal coupling strength between different parts of the system
	
	Matrix M(N);
	int* degree_dist = new int[N+1]; // counters for all degrees of vertices in the graph

#ifdef DETAILED_OUTPUT
	ofstream f(dateStr("_RandomGraph_N" + dtos((double)N,0) + ".txt")); // file containing parameter information
#endif
	ofstream data(dateStr("_RandomGraph.txt")); // file containing final stability data
	
	if (
#ifdef DETAILED_OUTPUT
		!f.is_open() ||
#endif
		!data.is_open())
	{
		cout << "Failed to open file!" << endl;
		return;
	}
	
	for (int i=0; i<=40; i+=conn_step) // loop over connectances
	{
		int stable = 0;
		int samples = sample_number;
		
		printf("\rC = %d %%", i);
		fflush(stdout);

		memset(degree_dist, 0, (N+1)*sizeof(int)); // set counters for all degrees to zero

		for (int j=0; j<sample_number; j++)
		{
			int ret;
#ifdef DETAILED_OUTPUT	
#ifdef NONZEROS_EXACT
			M.generateUniformRandom(double(i)/100., minS, maxS, minS_self, maxS_self, graph_struct, &f);
#else
			M.generateUniformRandomNonExact(double(i)/100., minS, maxS, minS_self, maxS_self, graph_struct, &f);
#endif // NONZEROS_EXACT

			if ((ret = M.isStable(&f)) == 1) 
			{
				stable++;
				f << "stable" << endl;
			}
			else if (ret == -1)
				f << "unstable" << endl;
			else
			{
				samples--; // this sample cannot be used
				f << "N/A" << endl;
			}
#ifdef MAT_OUTPUT
		M.saveMatrix(&f);
#endif
			for (int i=0; i<100; i++)
				f << "=";
			f << endl << endl;
#else
#ifdef NONZEROS_EXACT
			M.generateUniformRandom(double(i)/100., minS, maxS, minS_self, maxS_self, graph_struct);
#else
			M.generateUniformRandomNonExact(double(i)/100., minS, maxS, minS_self, maxS_self, graph_struct);
#endif
			if ((ret = M.isStable()) == 1) 
				stable++;
			else if (ret == 0)
				samples--; // this sample cannot be used
#endif // DETAILED_OUTPUT
			M.sumDegreeDist(degree_dist); // adds data from current matrix to the degree distribution
		}
		data.precision(1);
		data << fixed << i << "\t\t" << stable << "\t\t" << samples << endl;

		string fstr = "RandomGraph_dist_C" + dtos((double)i,0); // string describing the data
		plot_degree_dist(fstr.c_str(), "Degree distribution", degree_dist, N, N*sample_number);
	}
	cout << endl;
#ifdef DETAILED_OUTPUT
	f.close();
#endif
	data.close();
	delete[] degree_dist;
	
	plot_over_connectance("RandomGraph", "Results for random graph network", N);
}


/*** GardnerAshby ***
 * Reconstructs the results of Gardner & Ashby, 1970 with Erdős-Rényi(-Gilbert) random graphs *
 * - sample_number: number of samples per parameter setting *
 * - int graph_struct: indicates if the graph shall be freely directed as in the original study (2), or undirected (0), or asymmetrically bidirectional (1)  */
void GardnerAshby(int sample_number, int graph_struct)
{
	const int conn_step = 2; // step size for variation of connectances
	const double minS = -1., maxS = 1.; // minimal and maximal coupling strength between different parts of the system
	const double minS_self = -1., maxS_self = -0.1; // minimal and maximal self-coupling strength
	
	Matrix N4 (4);
	Matrix N7 (7);
	Matrix N10 (10);
#ifdef DETAILED_OUTPUT
	ofstream f1(dateStr("_GardnerAshby_N4.txt")); // file containing parameter information
	ofstream f2(dateStr("_GardnerAshby_N7.txt")); // file containing parameter information
	ofstream f3(dateStr("_GardnerAshby_N10.txt")); // file containing parameter information
#endif
	ofstream data(dateStr("_GardnerAshby.txt")); // file containing final stability data
	
	if (
#ifdef DETAILED_OUTPUT
		!f1.is_open() || !f2.is_open() || !f3.is_open() || 
#endif
		!data.is_open())
	{
		cout << "Failed to open file!" << endl;
		return;
	}
	
	for (int i=0; i<=100; i+=conn_step) // loop over connectances
	{
		int stable1 = 0, stable2 = 0, stable3 = 0;
		int samples1 = sample_number, samples2 = sample_number, samples3 = sample_number;
		
		printf("\rC = %d %%", i);
		fflush(stdout);

		for (int j=0; j<sample_number; j++)
		{
			int ret;
#ifdef DETAILED_OUTPUT
			ofstream* fout_handle_N4 = &f1;
			ofstream* fout_handle_N7 = &f2;
			ofstream* fout_handle_N10 = &f3;
#else
			ofstream* fout_handle_N4 = NULL;
			ofstream* fout_handle_N7 = NULL;
			ofstream* fout_handle_N10 = NULL;
#endif
			
#ifdef NONZEROS_EXACT // Erdős-Rényi-Gilbert random graph
			N4.generateUniformRandom(double(i)/100., minS, maxS, minS_self, maxS_self, graph_struct, fout_handle_N4);
#else // Erdős-Rényi random graph
			N4.generateUniformRandomNonExact(double(i)/100., minS, maxS, minS_self, maxS_self, graph_struct, fout_handle_N4);
#endif
			if ((ret = N4.isStable(fout_handle_N4)) == 1) 
			{
				stable1++;
				writeLineToFile("stable", fout_handle_N4);
			}
			else if (ret == -1)
				writeLineToFile("unstable", fout_handle_N4);
			else
			{
				samples1--; // this sample cannot be used
				writeLineToFile("N/A", fout_handle_N4);
			}
			N4.saveMatrix(fout_handle_N4);
			
#ifdef NONZEROS_EXACT // Erdős-Rényi-Gilbert random graph
			N7.generateUniformRandom(double(i)/100., minS, maxS, minS_self, maxS_self, graph_struct, fout_handle_N7);
#else // Erdős-Rényi random graph
			N7.generateUniformRandomNonExact(double(i)/100., minS, maxS, minS_self, maxS_self, graph_struct, fout_handle_N7);
#endif
			if ((ret = N7.isStable(fout_handle_N7)) == 1) 
			{
				stable2++;
				writeLineToFile("stable", fout_handle_N7);
			}
			else if (ret == -1)
				writeLineToFile("unstable", fout_handle_N7);
			else
			{
				samples2--; // this sample cannot be used
				writeLineToFile("N/A", fout_handle_N7);
			}
			N7.saveMatrix(fout_handle_N7);
		
#ifdef NONZEROS_EXACT // Erdős-Rényi-Gilbert random graph
			N10.generateUniformRandom(double(i)/100., minS, maxS, minS_self, maxS_self, graph_struct, fout_handle_N10);
#else // Erdős-Rényi random graph
			N10.generateUniformRandomNonExact(double(i)/100., minS, maxS, minS_self, maxS_self, graph_struct, fout_handle_N10);
#endif
			if ((ret = N10.isStable(fout_handle_N10)) == 1) 
			{
				stable3++;
				writeLineToFile("stable", fout_handle_N10);
			}
			else if (ret == -1)
				writeLineToFile("unstable", fout_handle_N10);
			else
			{
				samples3--; // this sample cannot be used
				writeLineToFile("N/A", fout_handle_N10);
			}
			N10.saveMatrix(fout_handle_N10);
		}
		data.precision(1);
		data << fixed << i << "\t\t" << stable1 << "\t\t" << samples1 << "\t\t" << stable2 << "\t\t" << samples2 << "\t\t" << stable3 << "\t\t" << samples3 << endl;
	}
	cout << endl;
#ifdef DETAILED_OUTPUT
	f1.close();
	f2.close();
	f3.close();
#endif
	data.close();
	
	plot_over_connectance("GardnerAshby", "Reproduction of the results of Gardner \\& Ashby, 1970", 4, 7, 10);
}


int main(int argc, char** argv) 
{ 
	dateStr("", true);
	string path = string("./SC_") + dateStr("", true); // path to working directory, get a new timestamp just to be sure that it is unique

	ofstream infof;
	string goal = ""; // goal of investigation (is not asked for if set to anyting but "")
	double minS_self = -1.3, maxS_self = -0.4; // minimal and maximal self-coupling strength
	int N = 10; // number of system constituents (seed graph constituents in case of Leskovec simulation)
	int p = 2; // Kronecker power (in case of Leskovec simulation)
	int n0 = 82; // number of off-diagonal zeros in seed matrix (in case of Leskovec simulation; should be even for undirected graphs)
	int trials = 1; // number of trials
	int simtype = SIMTYPE_SF; // type of the simulation (small world, scale-free, ...)
	double k = 4; // number of nearest-neighbor connections for initial regular graph (only for small-world networks)
	double p_rew = -1.; // if set, this fraction of connections is rewired (only for small-world networks)
	int m0 = 0; // size of initial vertex set (only for scale-free networks)
	int m = 2; // number of edges for each new vertex (only for scale-free networks)
	int stabilize_hubs = 0; // indicates if only the largest hubs, and how many of them, shall receive self-stabilization (0 to disable)
	int graph_struct = 1; // indicates if the considered graph shall be undirected (0), asymmetrically bidirectional (1), or fully directed (2)
	
	// process call arguments
	for(int i=1; i<argc; i++)
	{
		char* pt;
		if( (pt = strstr(argv[i], "=")) != NULL )
		{
			pt++;
		}
		
		if ( strstr(argv[i], "-minS_self=") == argv[i] || strstr(argv[i], "-f=") == argv[i] )
			minS_self = atof(pt);
		else if ( strstr(argv[i], "-maxS_self=") == argv[i] || strstr(argv[i], "-f=") == argv[i] )
			maxS_self = atof(pt);
		else if (strstr(argv[i], "-N=") == argv[i])
			N = atoi(pt);
		else if (strstr(argv[i], "-p=") == argv[i])
			p = atoi(pt);
		else if (strstr(argv[i], "-n0=") == argv[i])
			n0 = atoi(pt);
		else if (strstr(argv[i], "-trials=") == argv[i])
			trials = atoi(pt);
		else if (strstr(argv[i], "-smallworld") == argv[i])
			simtype = SIMTYPE_SW;
		else if (strstr(argv[i], "-scalefree") == argv[i])
			simtype = SIMTYPE_SF;
		else if (strstr(argv[i], "-leskovec") == argv[i])
			simtype = SIMTYPE_LE;
		else if (strstr(argv[i], "-randomgraph") == argv[i])
			simtype = SIMTYPE_RG;
		else if (strstr(argv[i], "-gardnerashby") == argv[i])
			simtype = SIMTYPE_GA;
		else if (strstr(argv[i], "-nogoal") == argv[i])
			goal = string("-");
		else if (strstr(argv[i], "-goal") == argv[i])
			goal = string(pt);
		else if (strstr(argv[i], "-k=") == argv[i])
			k = atoi(pt);
		else if (strstr(argv[i], "-p_rew=") == argv[i])
			p_rew = atof(pt);
		else if (strstr(argv[i], "-undirected=") == argv[i]) // (deprecated, 'graph_struct' should be used instead)
		{
			if ((bool) atoi(pt))
				graph_struct = 0;
		}
		else if (strstr(argv[i], "-graph_struct=") == argv[i])
			graph_struct = atoi(pt);
		else if (strstr(argv[i], "-m0=") == argv[i])
			m0 = atoi(pt);
		else if (strstr(argv[i], "-m=") == argv[i])
			m = atoi(pt);
		else if (strstr(argv[i], "-onlyhubs=") == argv[i])
			stabilize_hubs = atoi(pt);
	}
	
	mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); // create working directory

	if (chdir(path.c_str()) == -1) { // try to change directory
		showChDirErrMessage();
		return -1;
	}
	system("cp -r ../*.cpp ."); 
	system("cp -r ../*.hpp .");
	system("cp -r ../compile* ."); 
	
	infof.open("info.txt");
	if (!infof.is_open())
	{
		cout << "Failed to open info file!" << endl;
		return 0;
	}
	
	cout << "Simulation " << dateStr("") << endl;
	if (goal == "")
	{
		cout << "Goal of this investigation: "; 
		getline(cin, goal);
	}
	infof << "Goal of this investigation:" << endl << " " << goal << endl << endl 
	      << "Parameters:" << endl;
	      
	if (simtype != SIMTYPE_LE)
		infof << " N = " << N << endl;
	else
		infof << " N_seed = " << N << endl
		      << " p = " << p << endl
		      << " n0 = " << n0 << endl;
		
	infof << " minS_self = " << minS_self << endl
	      << " maxS_self = " << maxS_self << endl
	      << " graph_struct = ";
	switch (graph_struct)
	{
	case 0:
		infof << "0 (undirected)" << endl;
		break;
	case 1:
		infof << "1 (asymmetrically bidirectional)" << endl;
		break;
	case 2:
		infof << "2 (freely directed)" << endl;
		break;
	}
	infof << " Number of trials = " << trials << endl << endl;
	
	timeMeasure(true);

	if (graph_struct > 1
	   && simtype != SIMTYPE_RG && simtype != SIMTYPE_GA)
	{
		cout << "Error: freely directed graphs are not available for the selected type of simulation." << endl;
		return 0;
	}
	
	if (simtype == SIMTYPE_SW) // small world simulation
	{
		SmallWorld(N, minS_self, maxS_self, trials, k, p_rew, graph_struct);
	}
	else if (simtype == SIMTYPE_SF) // scale-free simulation
	{
		ScaleFree(N, minS_self, maxS_self, trials, m0, m, stabilize_hubs, graph_struct);
	}
	else if (simtype == SIMTYPE_LE) // Leskovec graph simulation
	{
		Leskovec(N, p, n0, minS_self, maxS_self, trials, graph_struct);
	}
	else if (simtype == SIMTYPE_RG) // random graph simulation
	{
		RandomGraph(N, minS_self, maxS_self, trials, graph_struct);
	}
	else if (simtype == SIMTYPE_GA) // reproduction of the simulations of Gardner and Ashby, 1970
	{
		GardnerAshby(trials, graph_struct);
	}
	
	
	int tsec = timeMeasure(false);
	
	char el_time [32];
	sprintf(el_time, "Elapsed time:\n %d min %d sec", int(floor(tsec / 60)), int(tsec % 60));
	cout << el_time << endl;
	
	infof << el_time << endl;
	infof.close();
}
