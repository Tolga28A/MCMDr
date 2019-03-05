#include "sim_box.h"

void Sim_Box::set_atm_coors(mat Coors) {
  num_atoms = Coors.n_rows; 
  atm_coors = Coors;
}

void Sim_Box::set_box_props(vec bp, mat bv) {
  box_period = bp; 
  box_vecs = bv;
}

mat Sim_Box::create_img_coors(double r_cut) {
  vec nImgs;
  nImgs = calc_num_imgs(r_cut);
  int nx = nImgs(0); int ny = nImgs(1); int nz = nImgs(2);
  
  vec avec, bvec, cvec;
  avec = box_vecs.col(1); bvec = box_vecs.col(2); cvec = box_vecs.col(3);
  
  mat Img_Coors = zeros<mat>(num_atoms*(2*nx+1)*(2*ny+1)*(2*nz+1), 3);

  int ct1 = 0; int ind_start = ct1*num_atoms; int ind_end = (ct1+1)*num_atoms-1;
  Img_Coors.rows(ind_start, ind_end) = atm_coors;
  ct1 = ct1+1; ind_start = ct1*num_atoms; ind_end = (ct1+1)*num_atoms-1;

  // Periodicity is taken into account by the values of nx, ny, nz
  // If the box is not periodic, e.g. in x-direction, then nx = 0
  for (int i=-nx; i<=nx; i++) {
    for (int j=-ny; j<=ny; j++) {
      for (int k=-nz; k<=nz; k++) {
	    if (i != 0 || j != 0 || k != 0) {
	      Img_Coors.submat(ind_start, 0, ind_end, 0) = atm_coors.col(0) + i*avec(0) + j*bvec(0) + k*cvec(0);
	      Img_Coors.submat(ind_start, 1, ind_end, 1) = atm_coors.col(1) + i*avec(1) + j*bvec(1) + k*cvec(1);
	      Img_Coors.submat(ind_start, 2, ind_end, 2) = atm_coors.col(2) + i*avec(2) + j*bvec(2) + k*cvec(2);
	      ct1 = ct1+1; ind_start = ct1*num_atoms; ind_end = (ct1+1)*num_atoms-1;
	    }
      }
    }
  }

  // Only keep atoms which are within the r_cut range of the main box.
  double x_min = min(atm_coors.col(0)) - 1.15*r_cut, x_max = max(atm_coors.col(0)) + 1.15*r_cut;
  double y_min = min(atm_coors.col(1)) - 1.15*r_cut, y_max = max(atm_coors.col(1)) + 1.15*r_cut;
  double z_min = min(atm_coors.col(2)) - 1.15*r_cut, z_max = max(atm_coors.col(2)) + 1.15*r_cut;

  uvec xind1 = (Img_Coors.col(0) >= x_min); uvec xind2 = (Img_Coors.col(0) <= x_max);
  uvec yind1 = (Img_Coors.col(1) >= y_min); uvec yind2 = (Img_Coors.col(1) <= y_max);
  uvec zind1 = (Img_Coors.col(2) >= z_min); uvec zind2 = (Img_Coors.col(2) <= z_max);
  uvec cond1 = find(xind1 + xind2 + yind1 + yind2 + zind1 + zind2 == 6);
  Img_Coors = Img_Coors.rows(cond1);
  
  return Img_Coors;
}

vec Sim_Box::calc_num_imgs(double r_cut) {
  vec nimgs = randi<vec>(3);
  vec avec, bvec, cvec;
  avec = box_vecs.col(1); bvec = box_vecs.col(2); cvec = box_vecs.col(3);
  if (box_period(0)) {
    nimgs(0) = calc_nproj (avec, bvec, cvec, r_cut);
  } else {
    nimgs(0) = 0;
  }
  if (box_period(1)) {
  nimgs(1) = calc_nproj (bvec, cvec, avec, r_cut);
    } else {
    nimgs(1) = 0;
  }
  if (box_period(2)) {
  nimgs(2) = calc_nproj (cvec, avec, bvec, r_cut);
    } else {
    nimgs(2) = 0;
  }
  return nimgs;
}

int Sim_Box::calc_nproj (vec vec1, vec vec2, vec vec3, double r_cut) {
	vec aproj = cross(vec2, vec3);
	aproj = aproj/norm(aproj);
	int a_min = ceil(r_cut/abs(dot(vec1, aproj)));
	return a_min;
}
