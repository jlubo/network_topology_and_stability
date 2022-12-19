#ifndef PLOTS_HPP
#define PLOTS_HPP

// Includes
#include <iostream>
#include "Tools.hpp"
using ::std::string;

// Function declarations
void plot_over_connectance(string, string, int, int = 0, int = 0);
void plot_over_rewirings(string, string, int, int, int);
bool plot_degree_dist(string, string, int*, int, int);

#endif
