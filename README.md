# Networks of different topology and their stability 


## Outline

The code provided here serves to randomly generate __networks of different topology__ (in particular, __small-world__ and __scale-free__ topologies) and to analyze their asymptotic stability. 
In particular, it enables to investigate the __stabilization of a network by adding self-coupling to its nodes__.

The code has been developed and used for the following publication:

   Luboeinski J, Claro L, Pomi A, Mizraji E. Stabilization through self-coupling in networks of small-world and scale-free topology. Sci. Rep. (2023). https://doi.org/10.1038/s41598-023-27809-8

If you use parts of the software package or the model for your research, please cite appropriately (see [here](BIBTEX.md) for BibTeX references). 
Please feel free to contact us for any questions.


## Build and run

This program uses [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/) v2.5 (it is therefore provided under the GNU General Public License v3). Please make sure that you have installed a suitable version of GSL before building the code. If you would like to have plots created automatically, also make sure that you have installed [gnuplot](http://www.gnuplot.info/).

The code comes with the following shell scripts for building:

* __build_main__ - builds a program to generate and stabilize networks (to reproduce the results of the paper mentioned above)
* __build_test__ - builds a test program to test specific aspects of the code

The code comes with the following shell scripts for running:

* __run_main__ - runs different program calls to generate and stabilize networks (to reproduce the results of the paper mentioned above; this will create "SC\_\*" data directories and may take some time)
* __run_test1__ - runs different program calls to quickly test the functionality of the code (this will create "SC\_\*" data directories)
* __run_test2__ - runs a test program to test specific aspects of the code (this will create a file "Test\_Main\_out.txt")


## License

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along with this program (in "LICENSE.txt").
If not, see <https://www.gnu.org/licenses/>.
