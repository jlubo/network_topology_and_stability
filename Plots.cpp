/********************************************
 ***  Functions for plotting via gnuplot  ***
 ***   (c) Jannik Luboeinski 2016-2022    ***
 ********************************************/

#include "Plots.hpp"

/*** plot_over_connectance ***
 * Creates a gnuplot plot of data that was created for several values of the connectance *
 * - name: filename for the plot file *
 * - title: title that appears in the plot *
 * - N: the system size *
 * - N2 [optional]: the system size *
 * - N3 [optional]: the system size */
void plot_over_connectance(string name, string title, int N, int N2, int N3)
{
	ofstream f(name + ".gpl");
	if (!f.is_open())
	{
		cout << "Failed to open plot file!" << endl;
		return;
	}
	f << "set terminal pdfcairo" << endl;
	f << "set output '" << name << ".pdf'" << endl;

	f << "set key top right" << endl;

	f << "set title \"" << title << "\"" << endl;
	f << "set size square" << endl;

	f << "set xlabel \"Connectance (%)\"" << endl;
	f << "set ylabel \"Probability of stability\"" << endl;
	
	f << "set xrange [0:100]" << endl;
	f << "set yrange [0:1]" << endl;

	f << "plot  \"" << dateStr("_") << name << ".txt\" using 1:($2/$3):(1/sqrt($3)) lc 1 title \"N=" << N << "\" with yerrorbars";
	if (N2 > 0)
		f << ", \\" << endl << "\t\"" << dateStr("_") << name << ".txt\" using 1:($4/$5):(1/sqrt($5)) lc 1 title \"N=" << N2 << "\" with yerrorbars";
	if (N3 > 0)
		f << ", \\" << endl << "\t\"" << dateStr("_") << name << ".txt\" using 1:($6/$7):(1/sqrt($7)) lc 1 title \"N=" << N3 << "\" with yerrorbars";
	
	f << endl;		
	
	f.close();

	system(concat("gnuplot " + name, ".gpl").c_str());
}

/*** plot_over_rewirings ***
 * Creates a gnuplot plot of data that was created for several values of the rewiring number *
 * - name: filename for the plot file *
 * - title: title that appears in the plot *
 * - N: the system size *
 * - min_rew_num: the smallest rewiring number used *
 * - max_rew_num: the largest rewiring number used */
void plot_over_rewirings(string name, string title, int N, int min_rew_num, int max_rew_num)
{
	ofstream f(name + ".gpl");
	if (!f.is_open())
	{
		cout << "Failed to open plot file!" << endl;
		return;
	}
	f << "set terminal pdfcairo" << endl;
	f << "set output '" << name << ".pdf'" << endl;

	f << "unset key" << endl;

	f << "set title \"" << title << "\"" << endl;
	f << "set size square" << endl;

	f << "set xlabel \"# of rewirings\"" << endl;
	f << "set ylabel \"Probability of stability\"" << endl;
	
	f << "set xrange [" << min_rew_num-1 << ":" << max_rew_num+1 << "]" << endl;

	f << "plot  \"" << dateStr("_") << name << ".txt\" using 1:($3/$4):(1/sqrt($4)) lc 1 title \"N=" << N << "\" with yerrorbars" << endl;	
	
	f.close();
	
	system(concat("gnuplot " + name, ".gpl").c_str());
}

/*** plot_degree_dist ***
 * Evaluates a given degree distribution, writes it to a data file, *
 * and creates a gnuplot plot of the degree distribution *
 * - fstr: string describing the data *
 * - title: title that appears in the plot *
 * - degree_dist: the degree distribtuion as an array of counters *
 * - N: the dimension of the graph (the maximum possible degree) *
 * - norm: normalization constant *
 * - return: true if successful, false if not */
bool plot_degree_dist(string fstr, string title, int* degree_dist, int N, int norm)
{
	int k_min = 0; // lowest degree value that has occurred (minimally possible: 0)
	int p_at_k_min = 0; // frequency of nodes with degree k_min
	int k_fit_min; // minimum degree value for power-law fit (for large N, equals the lowest degree value that has occurred)
	int k_max = N; // highest degree value that has occurred (maximally possible: N)
	int p_at_k_max = 0; // frequency of nodes with degree k_max
	
	ofstream distf(dateStr("_") + fstr + string(".txt"));
	if (!distf.is_open())
	{
		cout << "Failed to open file!" << endl;
		return false;
	}
	for (int i=0; i<=N; i++) 
	{
		if (degree_dist[i] > p_at_k_min) // encounter k_min
		{
			k_min = i;
			p_at_k_min = degree_dist[i];
		}
		
		if (degree_dist[i] > 0) // encounter k_max
		{
			k_max = i;
			p_at_k_max = degree_dist[i];
		}
		distf << i << "\t\t" << degree_dist[i] << endl;
	}
	distf.close();


	ofstream gplf(fstr + ".gpl");
	if (!gplf.is_open())
	{
		cout << "Failed to open plot file!" << endl;
		return false;
	}
	gplf << "set terminal pdfcairo" << endl;
	gplf << "set output '" << fstr << ".pdf'" << endl;

	gplf << "set key outside right reverse width -4 box" << endl;
	gplf << "set key samplen 1.5" << endl;
	
	gplf << "#set title \"" << title << "\"" << endl;
	gplf << "set size square" << endl;

	gplf << "set xlabel \"Node degree k\"" << endl;
	gplf << "set ylabel \"Relative frequency P(k)\"" << endl << endl;
	
	gplf << "norm = " << norm << endl;
	gplf << "k_min = " << k_min << endl;
	gplf << "k_max = " << k_max << endl;
	gplf << "gamma = 3" << endl;
	gplf << "#f(x) = C * x**(-gamma)" << endl;
	gplf << "f(x) = (gamma - 1) * k_min**(gamma-1) * x**(-gamma)" << endl << endl;
	
	gplf << "set xrange [" << -1 << ":k_max+1]" << endl;
	gplf << "set log x" << endl;
	gplf << "set format y \"\%.0e\"" << endl;
	gplf << "set log y" << endl << endl;

	gplf << "set fit quiet" << endl;
	gplf << "fit [k_min:" << N << "] f(x) \"" << dateStr("_") << fstr << ".txt\" using 1:($2/norm) via gamma" << endl << endl;

	gplf << "plot \"" << dateStr("_") << fstr << ".txt\" using 1:($2/norm) lc 1 title \"data\", \\" << endl
	     << "     f(x) title sprintf(\"fit with {/Symbol g}=\%.2f\\nand k_{min}=\%.2f\", gamma, k_min)" << endl;	
	
	gplf.close();
	
	system(concat("gnuplot " + fstr, ".gpl").c_str());

	return true;
}
