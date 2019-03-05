#include "monte_carlo.h"

void Monte_Carlo::Write_Lammps_In(int atom, int file, int atomrem, int flag, int count, double eps, double sigma)  {

    ofstream infile;
    infile.open ("mcmd.in");

    if (flag == 0)  {

        infile << "units	     metal" << endl;
        infile << "atom_style    atomic" << endl;
        infile << "dimension     3" << endl;
        infile << "boundary      p p p" << endl;
        infile << "timestep      0.001" << endl;
        infile << "variable      a equal " << atom << endl;
        infile << "variable      f equal " << file << endl;
        infile << "variable      i equal " << count << endl;
        infile << "variable      l equal " << count+1 << endl;
        infile << "log	         al_a$af$f.lammps" << endl;
        infile << "read_data     al_a$af$fc$i.dat" << endl;
        infile << "group         outside type 1" << endl;
        infile << "group         inside type 2" << endl;
        infile << "mass          * 26.98" << endl;
        infile << "pair_style    lj/cut 10" << endl;
        infile << "pair_coeff    * *  " << eps << " " << sigma << endl;        
        infile << "neighbor      1.5 bin" << endl;

        infile << "fix           NPTr inside rigid/npt single temp 1.0 1.0 50.0 iso 0.0 0.0 50.0" << endl;
        infile << "fix           NPT outside npt temp 1.0 1.0 10.0 iso 0.0 0.0 10.0" << endl;

        infile << "variable      ape equal pe/atoms" << endl;

        infile << "compute       atomenergyc all pe/atom" << endl;
        infile << "dump          atomenergyd outside custom 10000 atomenergy.dat id c_atomenergyc[*]" << endl;

        infile << "thermo        10000" << endl;
        infile << "thermo_style  custom step temp press etotal density v_ape" << endl;
        infile << "run           10000" << endl;

        infile << "print         '${ape}' file energy.dat" << endl;
        infile << "write_data    al_a$af$fc$i.dat nocoeff" << endl;
        infile << "write_data    al_a$af$fc$l.dat nocoeff" << endl;

    }

    else if (flag == 1)  {

        infile << "units	     metal" << endl;
        infile << "atom_style    atomic" << endl;
        infile << "dimension     3" << endl;
        infile << "boundary      p p p" << endl;
        infile << "timestep      0.001" << endl;
        infile << "variable      a equal " << atom << endl;
        infile << "variable      f equal " << file << endl;
        infile << "variable      i equal " << count << endl;
        infile << "variable      l equal " << count+1 << endl;
        infile << "log	         al_a$af$f.lammps" << endl;
        infile << "read_data     al_a$af$fc$i.dat" << endl;
	    infile << "group         remove id " << atomrem << endl;
        infile << "delete_atoms  group remove" << endl;
        infile << "group         outside type 1" << endl;
        infile << "group         inside type 2" << endl;
        infile << "pair_style    lj/cut 10" << endl;
        infile << "pair_coeff    * *  " << eps << " " << sigma << endl;
        infile << "mass          * 26.98" << endl;
        infile << "neighbor      1.5 bin" << endl;

        infile << "fix           NPTr inside rigid/npt single temp 1.0 1.0 50.0 iso 0.0 0.0 50.0" << endl;
        infile << "fix           NPT outside npt temp 1.0 1.0 10.0 iso 0.0 0.0 10.0" << endl;

        infile << "variable      ape equal pe/atoms" << endl;

        infile << "compute       atomenergyc all pe/atom" << endl;
        infile << "dump          atomenergyd outside custom 10000 atomenergy.dat id c_atomenergyc[*]" << endl;

        infile << "thermo        10000" << endl;
        infile << "thermo_style  custom step temp press etotal density v_ape" << endl;
        infile << "run           10000" << endl;

        infile << "print         '${ape}' append energy.dat" << endl;
        infile << "variable	     lenx equal lx" << endl;
        infile << "print         '${lenx}' append energy.dat" << endl;
        infile << "write_data    al_a$af$fc$l_rem.dat nocoeff" << endl;

    }           

    else  if (flag == 2) {

        infile << "units	     metal" << endl;
        infile << "atom_style    atomic" << endl;
        infile << "dimension     3" << endl;
        infile << "boundary      p p p" << endl;
        infile << "timestep      0.001" << endl;
        infile << "variable      a equal " << atom << endl;
        infile << "variable      f equal " << file << endl;
        infile << "variable      i equal " << count << endl;
        infile << "variable      l equal " << count+1 << endl;
        infile << "log	         al_a$af$f.lammps" << endl;
        infile << "read_data     al_a$af$fc$i_rem.dat" << endl;
        infile << "group         outside type 1" << endl;
        infile << "group         inside type 2" << endl;
        infile << "pair_style    lj/cut 10" << endl;
        infile << "pair_coeff    * *  " << eps << " " << sigma << endl;        
        infile << "mass          * 26.98" << endl;
        infile << "neighbor      1.5 bin" << endl;
        infile << "neigh_modify  every 1 delay 0 check yes" << endl;

        infile << "fix           NPTr inside rigid/npt single temp 1.0 1.0 50.0 iso 0.0 0.0 50.0" << endl;
        infile << "fix           NPT outside npt temp 1.0 1.0 10.0 iso 0.0 0.0 10.0" << endl;

        infile << "variable      ape equal pe/atoms" << endl;

        infile << "compute       atomenergyc all pe/atom" << endl;
        infile << "dump          atomenergyd outside custom 10000 atomenergy.dat id c_atomenergyc[*]" << endl;

        infile << "thermo        10000" << endl;
        infile << "thermo_style  custom step temp press etotal density v_ape" << endl;
        infile << "run           10000" << endl;

        infile << "print         '${ape}' file energy.dat" << endl;
        infile << "write_data    al_a$af$fc$i.dat nocoeff" << endl;
        infile << "write_data    al_a$af$fc$l.dat nocoeff" << endl;

    }                                                                                            

    infile.close();

}


int Monte_Carlo::Atom_Remove(vector<double>& maxen, vector<int>& numofa)  {

    vector<double> id;
    vector<double> en;
    int numv;
    int max;

    string line;
    ifstream readen;  
    readen.open("atomenergy.dat");

    while (true)  {

       readen >> line;

       if ( line == "ATOMS" )  { 
                    
           readen >> numv;
           numofa.push_back(numv);
           break;

        }
    }

    id.resize(numv);  en.resize(numv);

    while ( !readen.eof() )  {

       readen >> line;

       if ( line == "c_atomenergyc[*]" )  { 

           for (int i=0; i<numv; i++)  {

                readen >> id[i];
                readen >> en[i];

           } 
        }        
    }
    
    readen.close();

    auto greatest = max_element(en.begin(), en.end());
    maxen.push_back(*greatest);
    max = distance(begin(en), greatest);

    return id[max];

}

void Monte_Carlo::MC_Perturbation(int atom, int file, int atomrem, int flag, int count, vector<double>& maxen, vector<int>& numofa, double eps, double sigma)   {
    
    Write_Lammps_In(atom,file,atomrem,flag,count,eps,sigma); /// 1) Write LAMMPS in-script 
    system("nohup mpiexec -np 6 ./lmp_mpi < mcmd.in"); /// 2) Call LAMMPS
    atomrem = Atom_Remove(maxen,numofa); /// 3) Find atom to remove
    Write_Lammps_In(atom,file,atomrem,1,count,eps,sigma); /// 4) Write LAMMPS in-script
    system("nohup mpiexec -np 6 ./lmp_mpi < mcmd.in"); /// 5) Call LAMMPS
    
}

void Monte_Carlo::MC_Iteration(int atom, int file, int atomrem, vector<double>& lx, int& count, vector<double>& maxen, vector<int>& numofa, double r_cut, int check, double eps, double sigma)  {

//    ofstream lenfile;
//    lenfile.open ("le_f" +  to_string(file) + ".dat");
//    lenfile << "atom count length maxen numofao \n";

    int flag;
    while(true)  {

//      if ( count == 1 ) { break;}
      if ( check == 1 )  { // accept the atom remove and continue

         flag = 0;
         MC_Perturbation(atom, file, atomrem, 2, count, maxen, numofa, eps, sigma);

      }
      else  { // reject the atom remove and continue

         flag++;
         MC_Perturbation(atom, file, atomrem, 0, count, maxen, numofa, eps, sigma);

      }
           
      check =  MC_Energy(lx, mce); // MC check
//      lenfile << atom << " " << count << " " << lx[count] << " " << maxen[count] << " " << numofa[count] << " \n";

      if ( lx[count] <= 2.2*r_cut || flag == 10) { break; } // minimized, stop removing 
      count++;
   }
//   lenfile.close();

}

bool Monte_Carlo::MC_Energy(vector<double>& lenx, double& mce)  {

    ifstream energyfile;  
    energyfile.open("energy.dat");

    double e1, e2, lx, de;
    double boltze, kbt= 330*0.8614e-4; // kBT @ T=330 K in eV

    energyfile >> e1; // Ei
    energyfile >> e2; // Ei+1
    energyfile >> lx; // x-dimension
    lenx.push_back(lx);

    de = e2-e1;
    if (de < 0 ) { mce = e2; return true; }
    else { 

	    boltze = exp(-de/kbt);

        srand(time(NULL));
	    double re = rand() / (RAND_MAX + 1.);

	    if (re > boltze) { mce = e2; return true; }

	    else { mce = e1; return false; }

    }
}

void Monte_Carlo::Delete_Files(int atom, int file, int count, int minmc)  {

    string delet1;
    string delet2;
    string delet;
    for (int i=0; i<count+2; i++) {

       if (i != minmc) {

            if (i == 0) {

                delet = "rm al_a" + to_string(atom) + "f" + to_string(file) + "c" + to_string(i) + ".dat";
                system(delet.c_str());

            }

            else {
                        
                delet1 = "rm al_a" + to_string(atom) + "f" + to_string(file) + "c" + to_string(i) + ".dat";
                delet2 = "rm al_a" + to_string(atom) + "f" + to_string(file) + "c" + to_string(i) + "_rem.dat";
                system(delet1.c_str());
                system(delet2.c_str());
            }
       }
     }     
}



