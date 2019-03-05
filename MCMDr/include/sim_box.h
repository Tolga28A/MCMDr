//
//  sim_box.h
//  
//
//  Created by Srikanth Patala on 7/11/17.
//
//
/*! \file sim_box.h
   \brief Class for calculating the box geometry and image coordinates.

   This class reads the box properties from given dump files and and calculates the resulting image coordinates for the atomic environment calculations.
*/

#ifndef sim_box_h
#define sim_box_h

/** @cond */
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <armadillo>
/** @endcond */

using namespace std;
using namespace arma;

class Sim_Box {
public:
  /// Atom coordinates
  mat atm_coors;
  /// Number of atoms
  int num_atoms;
  /// The periodicity flags for each direction
  vec box_period;
  /// The supercell vectors
  mat box_vecs;

  /* Member Functions */
  /**
  / Get coordinates of the atoms in the simulation box.
  */
  void set_atm_coors(mat Coors);
  /**
  / Set box periodicity and box vectors (origin is (0,0,0))
  */
  void set_box_props(vec bp, mat bv);
  /**
  / Calculate the number of minimum number of images along box dimensions.
  */
  vec calc_num_imgs(double r_cut);
  /**
  / Create Images (3D, 2D or 1D depends on the values of nx, ny, nz).
  */
  mat create_img_coors(double r_cut);
  /**
  / Calculate the number of box-images required.
  */
  int calc_nproj (vec vec1, vec vec2, vec vec3, double r_cut);

};
	
#endif /*sim_box_h*/

