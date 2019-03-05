//
//  monte_carlo.h
//  
//
//  Created by Tolga Akiner on 09/11/2018
//
//
/*! \file monte_carlo.h
   \brief Class for executing hybrid Monte Carlo-Molecular Dynamics relaxation.

   This class executes the Monte-Carlo algorithm for relaxing the given atomic structure by removing the highest energy atoms. It writes input scripts for LAMMPS and runs it, and decides on atom removal in a stochastic-iterative manner until the system is relaxed.
*/

#ifndef monte_carlo_h
#define monte_carlo_h

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

class Monte_Carlo {
public:
   /// Counter for the MC perturbation
   int count; 
   /// Atom ID to be removed in a MC perturbation
   int atomrem; 
   /// Boolen value to accept/reject the perturbation
   int check; 
   /// Flag to stop the MC perturbations
   int flag; 
   /// Per-atom energy
   double mce; 
   /// The MC perturbation step with the minimum energy value of maximum-energy-atom 
   double minmc; 
   /// LJ parameters
   double eps, sigma; 

  /**
  / Write LAMMPS input script depending on the MC flag that is to accept/reject the atom removal.
  */    
  void Write_Lammps_In(int atom, int file, int atomrem, int flag, int count, double eps, double sigma);

  /**
  / Find the atom with maximum energy to remove.
  */
  int Atom_Remove(vector<double>& maxen, vector<int>& numofa);

  /**
  / Read the file including average atomic potential energy and the length of the simulation box and determine the acceptance/removal with a MC formulation.
  */
  bool MC_Energy(vector<double>& lenx, double& mce);

  /**
  / Single Monte-Carlo move. \n\n
  */
  void MC_Perturbation(int atom, int file, int atomrem, int flag, int count, vector<double>& maxen, vector<int>& numofa, double eps, double sigma);

  /**
  / Delete the .dat files that are not needed.
  */
  void Delete_Files(int atom, int file, int count, int minmc);

  /**
  / Iterative Monte-Carlo moves with atom removal/acceptance flag.
  */
  void MC_Iteration(int atom, int file, int atomrem, vector<double>& lx, int& count, vector<double>& maxen, vector<int>& numofa, double r_cut, int check, double eps, double sigma);

};

#endif /*monte_carlo_h*/

