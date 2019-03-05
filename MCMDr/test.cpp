
/** @cond */
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <dirent.h>
#include <glob.h>

#include "util_functions.h"
#include "sim_box.h"
#include "monte_carlo.h"
#include "mlpack_util.h"
/** @endcond */

/*! \file test.cpp
   \brief Main function for executing the code.

   Hybrid Monte Carlo-Molecular Dynamics relaxation scheme is executed in this function by calling all the other member (in classes) or non-member functions with the steps given below.
*/

int main()  {
	

  /// 1) Read the dump files in the directory and store them in a string-vector.
  string fdump = "dump.*";
  vector<string> files = globVector(fdump); 
  int ndump = files.size(); 

  /// 2) Input section
  // Cutoff distance for the atomic environment is two times the lattice_constant.
  int lat = 4.05; // input-1 --> Lattice constant
  int r_cut = 2*lat; // Reference might be Botu et al. (2017)
  double amass = 26.98; // input-2 --> Atomic mass
  
  double eps = 0.5, sigma = 2.755; // LJ parameters, I'm gonna link this input with new LJ-class!

  // Iterate over the dump files.
  for (int file=0; file<1; file++) {  

      int num_atoms; // Total # of atoms in a dump file.
      int num_attr; // Total # of attributes.
      mat box_vecs = zeros<mat>(3,4); // Box vectors.
      vec box_period = zeros<vec>(3); // Flags for box periodicity. 
      mat atm_coor, Atm_attr; // Array storing the coordinates and other attributes.
      uvec atm_ids; // Vector for the atom IDs.
      vector<string> attribute_labels; // String for the Label of attributes.

      /// 3) Read the dump file and store the data.
      read_dump_file(files[file], num_atoms, num_attr, box_period, box_vecs, attribute_labels, atm_coor, Atm_attr, atm_ids);

      /// 4) Create the image coordinates.
      Sim_Box s1;
      s1.set_atm_coors(atm_coor);
      s1.set_box_props(box_period, box_vecs);
      mat Img_Coors = s1.create_img_coors(r_cut);

      /// 5) Calculate the neighbor IDs and distances for the atomic environment of each atom.
      uvec t_ids = atm_ids; // Atom IDs in the dump file.
      vector<vector<size_t> > resultingNeighbors; // Neighbor IDs for each atom for the atomic environment.
      vector<vector<double> > resultingDistances; // Neighbor distances for each atom for the atomic environment.
      Find_Neighbors_mlpack(Img_Coors, t_ids, sqrt(3)*1.15*r_cut, resultingNeighbors, resultingDistances); // Find the neighbors for each atom.


      // The output file for the box length, maximum-energy-atom and the number of atoms (outside region) for the selected configuration.
      ofstream lenfile;
      lenfile.open ("le_f" +  to_string(file) + ".dat");
      lenfile << "atom count length maxen numofao \n";

      // Iterate over the atoms.
      for (int atom=0; atom<resultingNeighbors.size(); atom++)  {
//      for (int atom=0; atom<1; atom++)  {

        /// 6) Find the atomic environment coordinates for each atom.
        uvec t1_ids = conv_to< uvec >::from(resultingNeighbors[atom]);
        mat nn_coors = Img_Coors.rows(t1_ids);
        rowvec pt = Img_Coors.row(t_ids(atom));  
        mat Box_Coors = Find_Box_Neighbors(nn_coors, pt, r_cut);

        /// 7) Write the initial data file.
        string df_name = "al_a" + to_string(atom) + "f" + to_string(file) + "c0.dat";
        Write_Data_File(Box_Coors,r_cut,atom,file,df_name); 

        int count=0; // Counter for the MC perturbation
        int atomrem; // Atom ID to be removed in a MC perturbation
        int check; // Boolen value to accept/reject the perturbation
        int flag; // Flag to stop the MC perturbations
        double mce; // Per-atom energy
        int minmc; // The # of MC perturbation with min energy of max-energy-atom   

        // Vectors for max-energy-atom, # of atom, and sim box length for MC perturbations, will be used as output
        vector<double> maxen; 
        vector<int> numofa; 
        vector<double> lx; 

        //Monte Carlo (MC) class
        Monte_Carlo mc; 

        /// 8) Start with the initial MC move.
        mc.MC_Perturbation(atom, file, atomrem, 0, count, maxen, numofa, eps, sigma);

        /// 9) Boolen check for accept/reject the move and increase the MC counter.
        check =  mc.MC_Energy(lx, mce); // MC check
        count++;

        /// 10) Execute the MC iterations until it converges. 
        mc.MC_Iteration(atom,file,atomrem,lx,count,maxen,numofa,r_cut,check,eps,sigma);

        /// 11) Find the MC move number for the minimum energy configuration among all the moves.
        auto minmca = min_element(maxen.begin(), maxen.end());
        minmc = distance(begin(maxen), minmca);


        // Write the output file including box length, energy of max-energy-atom and #of atom at the outside group
        lenfile << atom << " " << minmc << " " << lx[minmc] << " " << maxen[minmc] << " " << numofa[minmc] << " \n";   

        /// 12) Delete the unused .dat files that are used during the relaxation.
        mc.Delete_Files(atom,file,count,minmc);

        /// 13) Write the .cfg file for MLIP in the desired file structure.
        string fdel = "al_a" + to_string(atom) + "f" + to_string(file) + "c" + to_string(minmc) + ".dat";
        Write_POSCAR(atom,file,minmc,fdel); 

     }     
                             
     lenfile.close();     
                   
  }

  return 0;

}

