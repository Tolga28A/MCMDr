//
//  mlpack_util.h
//
//
//  Created by Srikanth Patala on 11/15/17.
//
//
/*! \file mlpack_util.h
   \brief Neighbor atom calculations for the atomic environment.

   These functions find the neighbor atoms and the corresponding coordinates for an atomic environment using RangeSearch function of MLpack library.
*/

#ifndef mlpack_util_h
#define mlpack_util_h

/** @cond */
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>

// For MATLAB-like functionalities
#include <armadillo>

// For near-neighbor calculations
#include <mlpack/core.hpp>
#include <mlpack/methods/neighbor_search/neighbor_search.hpp>
#include <mlpack/methods/range_search/range_search.hpp>
/** @endcond */

using namespace arma;
using namespace std;

/**
/ Find all the neighbors for the atomic environment.
*/
void Find_Neighbors_mlpack(mat Img_Coors, uvec atm_inds, double r_cut,
			   vector<vector<size_t> >& resultingNeighbors,
			   vector<vector<double> >& resultingDistances);

/**
/ Find the atomic environment coordinates.
*/
mat Find_Neighbors_mlpack(mat Img_Coors, rowvec pt, double r_cut);

#endif /* mlpack_util_h */
