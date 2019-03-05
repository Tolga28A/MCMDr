//
//  util_functions.h
//
//
//  Created by Srikanth Patala on 11/15/17.
//
//
/*! \file util_functions.h
   \brief Other functions

   This header consists of functions with different purposes as will be mentioned below.
*/


#ifndef util_functions_h
#define util_functions_h

/** @cond */
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <iomanip>
#include <dirent.h>
#include <glob.h>
#include <bits/stdc++.h>

// For MATLAB-like functionalities
#include <armadillo>
/** @endcond */

using namespace arma;
using std::vector;
using namespace std;

/**
/ Read the dump file for all the atom attributes and the coordinates.
*/
void read_dump_file(string fname, int &num_atoms, int &num_attr, vec &box_period,
		    mat &box_vecs, vector<string>& attribute_labels, mat &atm_coor,
		    mat &Atm_attr, uvec &atm_ids);

/**
/ Calculate the atomic environment for a given atom.
*/
mat Find_Box_Neighbors(mat nn_coors, rowvec pt, double r_cut);
vector<string> split(string str, char delimiter);
vector<string> globVector(const string& pattern);

/**
/ Determines the inside (rigid) and outside atoms for the relaxation.
/ Eliminates the atoms which are very close to box boundary to avoid overlapping.
/ Prepares the atom IDs, types and coordinates for the .dat file.
*/

mat Add_ID_Type(mat Box_Coors,  double r_cut);

/**
/ Write the data file of atomic environment for LAMMPS.
*/
void Write_Data_File(mat Box_Coors, double r_cut, int atom, int file, string df_name);

/**
/ Write the resulting (relaxed) configurations in POSCAR files for VASP. \n\n

*/
void Write_POSCAR(int atom, int file, int minmc, string fdel);

/**
/ Write the resulting (relaxed) configurations in a .cfg file for MLIP. \n\n

*/
void Write_ML(int atom, int file, int minmc, string fdel);

#endif /* util_functions_h */
