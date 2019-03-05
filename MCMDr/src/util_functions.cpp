#include "util_functions.h"

// Routine to split a string (useful for reading lines from files)
vector<string> split(string str, char delimiter) {
    vector<string> internal;
    stringstream ss(str); // Turn the string into a stream.
    string tok;
    while(getline(ss, tok, delimiter)) {
        internal.push_back(tok);
    }
	// cout << internal.size() << endl;
    return internal;
}

vector<string> globVector(const string& pattern){
    glob_t glob_result;
    glob(pattern.c_str(),GLOB_TILDE,NULL,&glob_result);
    vector<string> files;
    for(unsigned int i=0;i<glob_result.gl_pathc;++i){
        files.push_back(string(glob_result.gl_pathv[i]));
    }
    globfree(&glob_result);
    return files;
}

void read_dump_file(string fname, int &num_atoms, int &num_attr, vec &box_period,
		    mat &box_vecs, vector<string>& attribute_labels, mat &atm_coor,
		    mat &Atm_attr, uvec &atm_ids)
{
  ifstream data_file;
  data_file.open(fname);
  // data_file.open(fname.c_str());
  if(!data_file.is_open())
  {
    std::cerr << "Unable to open input file " << fname << "\n\n";
    exit(-1);
  }

  std::ifstream& fp = data_file;

  double origin[3]; origin[0]=0.; origin[1]=0.; origin[2]=0.;
  double supercell_edges[3][3];
  for(int c=0; c<3; c++)
    for(int d=0; d<3; d++)
      supercell_edges[c][d]=0;

  std::string full_line;
  getline(fp, full_line);
  // std::cout << full_line << std::endl;

  bool done=0;

  // Output Variables
  bool aPeriod, bPeriod, cPeriod;
  double hi_bound[3];
  int  number_of_particles;
  // std::vector<std::string> attribute_labels;
  int  xindex;
  bool scaled_coordinates;

  while(done==0)
  {
    getline(fp, full_line);
    // std::cout << full_line[0:10] << std::endl;
    std::istringstream iss(full_line);      

    if (full_line.find("ITEM: NUMBER OF ATOMS") != std::string::npos)
    {
      // std::cout << fp << std::endl;
      fp >> number_of_particles;
      // std::cout << "num particles: " << number_of_particles << std::endl;
    }
      
    else if (full_line.find("ITEM: BOX BOUNDS") != std::string::npos)          // PERIODIC IN ALL THREE DIRECTIONS
    {
      // cout << "Check: " << endl;
      // cout << (full_line.find("ITEM: BOX BOUNDS xy xz yz") != std::string::npos) << endl;
      if (full_line.find("ITEM: BOX BOUNDS xy xz yz") != std::string::npos) // TRICLINIC CRYSTAL SYSTEM
      {
        // std::cout << full_line << std::endl;
        // std::cout << full_line.substr(26,9) << std::endl;
        // X-Periodicity
        if (full_line.substr(26,2) == "pp") aPeriod = 1;
        else aPeriod = 0;
        // Y-Periodicity
        if (full_line.substr(29,2) == "pp") bPeriod = 1;
        else bPeriod = 0;
        // Z-Periodicity
        if (full_line.substr(32,2) == "pp") cPeriod = 1;
        else cPeriod = 0;

        // std::cout << "Peridicity: a, b, c: " << aPeriod << "\t" << bPeriod << "\t"<< cPeriod << "\t" << std::endl;
        double xlo_bound, xhi_bound, xy;
        double ylo_bound, yhi_bound, xz;
        double zlo_bound, zhi_bound, yz;

        fp >> xlo_bound; fp >> xhi_bound; fp >> xy;
        fp >> ylo_bound; fp >> yhi_bound; fp >> xz;
        fp >> zlo_bound; fp >> zhi_bound; fp >> yz;

        // cout << "xlo_bound: " << xlo_bound << "\t xhi_bound: " << xhi_bound << "\t xy: " << xy << endl;
        // cout << "ylo_bound: " << ylo_bound << "\t yhi_bound: " << yhi_bound << "\t xz: " << xz << endl;
        // cout << "zlo_bound: " << zlo_bound << "\t zhi_bound: " << zhi_bound << "\t yz: " << yz << endl;

        origin  [2] = zlo_bound; hi_bound[2] = zhi_bound;
        origin  [1] = ylo_bound - fmin(0.,yz); hi_bound[1] = yhi_bound - fmax(0.,yz);
        origin  [0] = xlo_bound - fmin(fmin(0.0,xy),fmin(xz,xy+xz));
        hi_bound[0] = xhi_bound - fmax(fmax(0.0,xy),fmax(xz,xy+xz));

        supercell_edges[1][0] = xy; supercell_edges[2][0] = xz; supercell_edges[2][1] = yz;
            
        for(int c=0; c<3; c++)
        {
          supercell_edges[c][c] = hi_bound[c]-origin[c];
        }
      }
      else
      {
        // std::cout << full_line << std::endl;
        // std::cout << full_line.substr(17,2) << std::endl;
        // X-Periodicity
        if (full_line.substr(17,2) == "pp") aPeriod = 1;
        else aPeriod = 0;
        // Y-Periodicity
        if (full_line.substr(20,2) == "pp") bPeriod = 1;
        else bPeriod = 0;
        // Z-Periodicity
        if (full_line.substr(23,2) == "pp") cPeriod = 1;
        else cPeriod = 0;

        for(int c=0; c<3; c++)
        {
          fp >> origin[c];
          fp >> hi_bound[c];
          supercell_edges[c][c] = hi_bound[c]-origin[c];
        }
      }
    }
    else if (full_line.find("ITEM: ATOMS ") != std::string::npos)
    {
      // DETERMINE ATOM ATTRIBUTES, INDEX OF COORDINATES, AND SCALING
      std::string entry;
      iss >> entry;   // READS IN "ITEM:"
      iss >> entry;   // READS IN "ATOMS"
      
      while(iss >> entry)                             // READ ATTRIBUTE LABELS, I.E., 'id', 'type', 'x', 'y', and 'z'
        attribute_labels.push_back(entry);
      

      xindex=-1;
      for(int c=0; c<attribute_labels.size(); c++)
      {
        if(attribute_labels[c]=="x" || attribute_labels[c]=="xs")
        {
          if(attribute_labels[c]=="x") scaled_coordinates=0;
          else                         scaled_coordinates=1;

          attribute_labels.erase(attribute_labels.begin()+c, attribute_labels.begin()+c+3);
          xindex=c;
        }
      }
      
      if(xindex==-1)
      {
        std::cout << "Insufficient xyz coordinates included in dump file.\n";
        exit(-1);
      }
      done=1;
    }
  }

  // // Print SuperCell Edges
  // for (int i=0; i<3; i++){
  //  for (int j=0; j<3; j++){
  //    std::cout << "i = " << i << "\t j = " << j << "\t :" << supercell_edges[i][j] << std::endl;
  //  }
  // }
  
  
  // for (int i=0; i<attribute_labels.size(); i++) {
  //  std::cout << attribute_labels[i] << std::endl;
  // }
  // std::cout << "Xindex: " << xindex << std::endl;

  // STORES ALL PARTICLE DATA TO particle_data[][]
  int entries = attribute_labels.size();
  vector<vector <double> > particle_data; // STORE PARTICLE DATA
  particle_data.resize(number_of_particles, std::vector<double>(entries));  // DATA FOR OUTPUT
    
  std::vector <double> xcoord; std::vector <double> ycoord; std::vector <double> zcoord;
  xcoord.resize (number_of_particles);
  ycoord.resize (number_of_particles); zcoord.resize (number_of_particles);

  num_atoms = number_of_particles;
  atm_coor = zeros<mat>(num_atoms,3);
  // atm_coor.col(0) = vec(xcoord);
  // atm_coor.col(1) = vec(ycoord);
  // atm_coor.col(2) = vec(zcoord);

  // cout << scaled_coordinates << endl;
  
  for(int c=0; c<number_of_particles; c++)
    {
      double x,y,z;

      for(int d=0; d<xindex; d++)
        fp >> particle_data[c][d];

      fp >> xcoord[c]; x = xcoord[c];
      fp >> ycoord[c]; y = ycoord[c];
      fp >> zcoord[c]; z = zcoord[c];
      
      // std::cout << "x: " << x << "\t y: " << y << "\t z: " << z << std::endl;

      for(int d=xindex+3; d<attribute_labels.size()+3; d++)
        fp >> particle_data[c][d-3];
      
      // ADJUST COORDINATES SO SYSTEM UNSCALED AND AT ORIGIN
      // if(scaled_coordinates==0)
      // 	{
      // 	  x -= origin[0]; y -= origin[1]; z -= origin[2];
      // 	}
      
      if(scaled_coordinates==1)
	{
	  double newx = supercell_edges[0][0]*x + supercell_edges[1][0]*y + supercell_edges[2][0]*z;
	  double newy = supercell_edges[0][1]*x + supercell_edges[1][1]*y + supercell_edges[2][1]*z;
	  double newz = supercell_edges[0][2]*x + supercell_edges[1][2]*y + supercell_edges[2][2]*z;
	  x=newx; y=newy; z=newz;
	}
      atm_coor(c,0) = x; atm_coor(c,1) = y; atm_coor(c,2) = z;
    }

  //Assign Outputs

  box_period(0) = aPeriod; box_period(1) = bPeriod; box_period(2) = cPeriod;

  for (int j=0; j<3; j++)
      box_vecs(j,0) = origin[j];

  for (int i=1; i<4; i++)
    for (int j=0; j<3; j++)
      box_vecs(j,i) = supercell_edges[i-1][j];

  // cout << box_vecs << endl;
  
  ///////////////////////
  num_attr = attribute_labels.size();

  Atm_attr = randu<mat>(num_atoms, num_attr);
  for (int i=0; i<attribute_labels.size(); i++)
  {
    for (int j=0; j<num_atoms; j++)
    {
      Atm_attr(j,i)=particle_data[j][i];
    }
  }
  
  
  ////////////////////////////////////////////////////////////////////////
  atm_ids = conv_to< uvec >::from(Atm_attr.col(0));
  uvec sort_inds = sort_index(atm_ids);
  // cout << "Atm_attr num rows  = " << Atm_attr.n_rows << endl;
  // cout << "Atm_IDs num rows  = " << atm_ids.n_rows << endl;
  // cout << sort_index(atm_ids) << endl;
  // cout << Atm_attr << endl;
  
  // Atm_attr.rows(sort_index(atm_ids)) = Atm_attr;
  Atm_attr = Atm_attr.rows(sort_inds);
  // atm_coor.rows(sort_index(atm_ids)) = atm_coor;
  atm_coor = atm_coor.rows(sort_inds);
  // atm_ids.rows(sort_index(atm_ids)) = atm_ids;
  atm_ids = atm_ids(sort_inds);
  // cout << atm_ids << endl;


}

mat Find_Box_Neighbors(mat nn_coors, rowvec pt, double r_cut) {

  vec diff_x = abs(nn_coors.col(0) - pt(0));
  vec diff_y = abs(nn_coors.col(1) - pt(1));
  vec diff_z = abs(nn_coors.col(2) - pt(2));

  uvec b_ids = find((diff_x <= 1.15*r_cut) && (diff_y <= 1.15*r_cut) && (diff_z <= 1.15*r_cut));

  mat Box_Coors = nn_coors.rows(b_ids);
  Box_Coors.col(0) = Box_Coors.col(0) - pt(0);
  Box_Coors.col(1) = Box_Coors.col(1) - pt(1);
  Box_Coors.col(2) = Box_Coors.col(2) - pt(2);

  return Box_Coors;

}

mat Add_ID_Type(mat Box_Coors, double r_cut) {

  vec diff = sqrt(Box_Coors.col(0)%Box_Coors.col(0)+Box_Coors.col(1)%Box_Coors.col(1)+Box_Coors.col(2)%Box_Coors.col(2));

  double gap = 0.8; // This is for introducing a gap between edge atoms and the box boundaries to eliminate overlapping
                    // I can calculate all the distances to the boundary and eliminate the ones which are closer than a 
                    // threshold distance, but it is easier to do it like this now
  double lbound = -1.15*r_cut+gap;
  double hbound = 1.15*r_cut-gap;

  uvec in_ids = find(diff <= r_cut);
  uvec out_ids = find(diff > r_cut);

  mat Coors = zeros<mat>(Box_Coors.n_rows,5);

  int k = 0;
  for (int i=0; i<Box_Coors.n_rows; i++)  {             

     if (Box_Coors(i,0) < hbound && Box_Coors(i,0) > lbound && Box_Coors(i,1) < hbound && Box_Coors(i,1) > lbound && Box_Coors(i,2) < hbound && Box_Coors(i,2) > lbound) 
     {

         Coors(k,0) = k+1;

         if ( any(i == in_ids) ) {

            Coors(k,1) = 2;

         }

         else if ( any(i == out_ids) ) {

            Coors(k,1) = 1;

         }
         
         Coors(k,2) = Box_Coors(i,0); 
         Coors(k,3) = Box_Coors(i,1); 
         Coors(k,4) = Box_Coors(i,2);
         k++;

     }
  }

  Coors.resize(k,5);
  return Coors;

}

void Write_Data_File(mat Box_Coors, double r_cut, int atom, int file, string df_name)  {

  /// Add atom IDs and types to the data, and trim the boundary atoms at the edges
  mat Coors = Add_ID_Type(Box_Coors, r_cut);

  ofstream newfile;
  newfile.setf(ios_base::fixed);  
  newfile << setprecision(5);   
  newfile.open (df_name);

  newfile << "#For atom number " << atom << " in file " << file << "\n" << endl;
  newfile << "\t" << Coors.n_rows << " atoms" << endl;
  newfile << "\t" << 0 << " bonds" << endl;
  newfile << "\t" << 0 << " angles" << endl;
  newfile << "\t" << 0 << " dihedrals" << endl;
  newfile << "\t" << 0 << " impropers \n" << endl;

  newfile << "\t" << 2 << " atom types" << endl;
  newfile << "\t" << 0 << " bond types" << endl;
  newfile << "\t" << 0 << " angle types" << endl;
  newfile << "\t" << 0 << " dihedral types" << endl;
  newfile << "\t" << 0 << " improper types \n" << endl;

  newfile << "\t" << -1.15*r_cut << " " << 1.15*r_cut << " xlo xhi" << endl;  
  newfile << "\t" << -1.15*r_cut << " " << 1.15*r_cut << " ylo yhi" << endl; 
  newfile << "\t" << -1.15*r_cut << " " << 1.15*r_cut << " zlo zhi \n" << endl;  

  newfile << "Atoms \n\n";


  Coors.raw_print(newfile);

  newfile.close();

}

void Write_ML(int atom, int file, int minmc, string fdel)  {

    /// 1) Take the selected .dat file and write the configurations, energy and stress in a new dump file using LAMMPS.

    ofstream infile;
    infile.open ("mcmdml.in");

    infile << "units	     metal" << endl;
    infile << "atom_style    atomic" << endl;
    infile << "dimension     3" << endl;
    infile << "boundary      p p p" << endl;
    infile << "timestep      0.001" << endl;
    infile << "variable      a equal " << atom << endl;
    infile << "variable      f equal " << file << endl;
    infile << "variable      c equal " << minmc << endl;

    infile << "log	         al_a$af$f_mldump.lammps" << endl;
    infile << "read_data     " << fdel << endl;
    infile << "mass          * 26.98" << endl;
    infile << "group         outside type 1" << endl;
    infile << "group         inside type 2" << endl;
    infile << "pair_style    lj/cut 8" << endl;
    infile << "pair_coeff    * *  0.5 2.755" << endl;
    infile << "neighbor      0.5 bin" << endl;
    infile << "neigh_modify  every 1 delay 0 check yes" << endl;
    infile << "fix           NPTr inside rigid single" << endl;
    infile << "fix           NPT outside npt temp 1.0 1.0 5.0 iso 0.0 0.0 1000.0" << endl;

    infile << "thermo        0" << endl;
    infile << "thermo_style  custom step temp press etotal density pxx pyy pzz pyz pxz pxy" << endl;
    infile << "run           0" << endl;

    infile << "variable      energy equal etotal" << endl;
    infile << "variable      sxx equal pxx" << endl;
    infile << "variable      syy equal pyy" << endl;
    infile << "variable      szz equal pzz" << endl;
    infile << "variable      syz equal pyz" << endl;
    infile << "variable      sxz equal pxz" << endl;
    infile << "variable      sxy equal pxy" << endl;

    infile << "print          '${energy}' file es.dat" << endl;
    infile << "print          '${sxx}' append es.dat" << endl;
    infile << "print          '${syy}' append es.dat" << endl;
    infile << "print          '${szz}' append es.dat" << endl;
    infile << "print          '${syz}' append es.dat" << endl;
    infile << "print          '${sxz}' append es.dat" << endl;
    infile << "print          '${sxy}' append es.dat" << endl;

    infile << "write_dump all custom al_a$af$f_mlip.dat id type x y z fx fy fz" << endl;

    infile.close();
	
    system("nohup mpiexec -np 4 ./lmp_mpi < mcmdml.in");    

    /// 2) Read the data from that dump file.
    int numofa;
    double xlow, xhigh, ylow, yhigh, zlow, zhigh;
    vector<vector <double> > data; 

    int l=0;
    double energy, sxx, syy, szz, syz, sxz, sxy;

    string line;
    ifstream readfile;  
    readfile.open("al_a" + to_string(atom) + "f" + to_string(file) + "_mlip.dat");

    while (true)  {

       readfile >> line;

       if ( line == "ATOMS" )  { 
                    
           readfile >> numofa;
           break;

       }
    }
    data.resize(numofa,vector<double>(8));

    while (true)  {

       readfile >> line;

       if ( line == "pp" )  { 

           readfile >> line;                        
           readfile >> line;  

		   readfile >> xlow;
		   readfile >> xhigh;
		   readfile >> ylow;
		   readfile >> yhigh;
		   readfile >> zlow;
		   readfile >> zhigh;
           break;             

        }        
    }

    while (true)  {

        readfile >> line; 

        if ( line == "fz" ) {break;}

    }

    while (true) {

       for (int j = 0; j < 8; j++) {

         readfile >> data[l][j];

       }
       l++;
       if (l == numofa) {break;}
    }

    readfile.close();

    /// 3) Calculate the new supercell vectors.
    double x = xhigh-xlow;
    double y = yhigh-ylow;
    double z = zhigh-zlow;

    /// 4) Shift the xyz coordinates with respect to the supercell vectors defined
    for (int i = 0; i < l; i++)  {
       
        data[i][2] = data[i][2] - xlow;
        data[i][3] = data[i][3] - ylow;
        data[i][4] = data[i][4] - zlow; 
       
    }

    /// 5) Read the energy and stress for the cfg. file.
    ifstream efile;  
    efile.open("es.dat");

    efile >> energy; 
    efile >> sxx; 
    efile >> syy; 
    efile >> szz; 
    efile >> syz; 
    efile >> sxz; 
    efile >> sxy; 

    efile.close();

    ///6)  Write the .cfg file.
    ofstream newfile;
    newfile.open ("mlip.cfg", ios_base::app);

    newfile << "\n BEGIN_CFG" << endl;
    newfile << " Size" << endl;
    newfile << "\t" << numofa << endl;
    newfile << " SuperCell" << endl;
    newfile << " \t\t" << x << "\t 0.0 \t 0.0" << endl;
    newfile << " \t\t" << "0.0 \t " << y << "\t 0.0" << endl;
    newfile << " \t\t" << "0.0 \t " << "0.0 \t " << z << "\n" << endl;
    newfile << "AtomData: \t id \t type \t cartes_x \t cartes_y \t cartes_z \t fx \t fy \t fz" << endl;

    for (int i = 0; i < l; i++) {

          newfile << "\t \t" << i+1 << " ";
          newfile << 0 << " ";
          newfile << data[i][2] << "\t";
          newfile << data[i][3] << "\t";
          newfile << data[i][4] << "\t";
          newfile << data[i][5] << "\t";
          newfile << data[i][6] << "\t";
          newfile << data[i][7] << "\n";

    }

    newfile << " Energy" << endl;
    newfile << "\t" << energy << endl;
    newfile << " Stress: xx \t yy \t zz \t yz \t xz \t xy" << endl;
    newfile << "\t" << sxx << "\t" << syy << "\t" << szz << "\t" << syz << "\t" << sxz << "\t" << sxy << "\n" << endl;
    newfile << " Feature   from  database:mlip.cfg" << endl;

    newfile << "END_CFG \n" << endl;
//    newfile.close();

}  



void Write_POSCAR(int atom, int file, int minmc, string fdel)  {


    /// 1) Read the data from dump file
    int numofa;
    double xlow, xhigh, ylow, yhigh, zlow, zhigh;
    vector<vector <double> > data; 

    int l=0;

    string line;
    double dum;
    ifstream readfile;  
    readfile.open(fdel);

    while (true)  {

       readfile >> line;

       if ( line == "=" )  { 
                    
           readfile >> dum;
           readfile >> numofa;
           break;

       }
    }
    data.resize(numofa,vector<double>(3));
    
    while (true)  {

       readfile >> line;

       if ( line == "types" )  { 
                    
           readfile >> xlow;
           readfile >> xhigh;
           readfile >> line;
           readfile >> line;
           readfile >> ylow;
           readfile >> yhigh;
           readfile >> line;
           readfile >> line;
           readfile >> zlow;
           readfile >> zhigh;                    
           break;

       }
    }

    while (true)  {

       readfile >> line;

       if ( line == "atomic" )  { 

           for (int i=0; i<numofa; i++)  {

               readfile >> dum;
               readfile >> dum;
               readfile >> data[i][0];
               readfile >> data[i][1];
               readfile >> data[i][2];
               readfile >> dum;
               readfile >> dum;
               readfile >> dum;

           }
           break;             
        }        
    }

    readfile.close();

    /// 2) Calculate the new supercell vectors.
    double x = xhigh-xlow;
    double y = yhigh-ylow;
    double z = zhigh-zlow;


    /// 3) Shift the xyz coordinates with respect to the supercell vectors defined
    for (int i = 0; i < numofa; i++)  {
       
        data[i][0] = data[i][0] + xhigh;
        data[i][1] = data[i][1] + yhigh;
        data[i][2] = data[i][2] + zhigh; 
       
    }


    /// 4)  Write the POSCAR file.
    ofstream newfile;
    newfile.open ("/home/tolga/Documents/github/cpp_analysis_codes-master/build/Poscars/LJ/POSCAR"+to_string(atom));

    newfile << "Al_LJ_relaxed\n" << 1 << endl;
    newfile << x << " " << 0.0000 << " " << 0.0000 << endl;
    newfile << 0.0000 << " " << y << " " << 0.0000 << endl;
    newfile << 0.0000 << " " << 0.0000 << " " << z << endl;

    newfile << numofa << "\nc" << endl;

    for (int i=0; i<numofa; i++)   {

        if (i==0)   {newfile << data[i][0] << " " << data[i][1] << " " << data[i][2] << " positions" << endl;}
        else {newfile << data[i][0] << " " << data[i][1] << " " << data[i][2] << endl;}
    }

    newfile.close(); 

}  
