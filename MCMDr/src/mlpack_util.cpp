#include "mlpack_util.h"

void Find_Neighbors_mlpack(mat Img_Coors, uvec atm_inds, double r_cut,
			   vector<vector<size_t> >& resultingNeighbors,
			   vector<vector<double> >& resultingDistances)
{
  // The vector-of-vector objects we will store output in.
  // std::vector<std::vector<size_t> > resultingNeighbors;
  // std::vector<std::vector<double> > resultingDistances;

  mat queryData = strans(Img_Coors.rows(atm_inds));
  mat referenceData = strans(Img_Coors);

  mlpack::range::RangeSearch<> a(referenceData);

  // The range we will use.
  mlpack::math::Range r(0.0, r_cut); // [0.0, r_cut].
  a.Search(queryData, r, resultingNeighbors, resultingDistances);
}

mat Find_Neighbors_mlpack(mat Img_Coors, rowvec pt, double r_cut)
{
  // The vector-of-vector objects we will store output in.
  std::vector<std::vector<size_t> > resultingNeighbors;
  std::vector<std::vector<double> > resultingDistances;

  colvec queryData = {pt(0), pt(1), pt(2)};
  mat referenceData = strans(Img_Coors);

  mlpack::range::RangeSearch<> a(referenceData);

  // The range we will use.
  mlpack::math::Range r(0.0, r_cut); // [0.0, r_cut].
  a.Search(queryData, r, resultingNeighbors, resultingDistances);

  uvec t1_ids = conv_to< uvec >::from(resultingNeighbors[0]);
  mat Sph_Coors = Img_Coors.rows(t1_ids);
  Sph_Coors.col(0) = Sph_Coors.col(0) - pt(0);
  Sph_Coors.col(1) = Sph_Coors.col(1) - pt(1);
  Sph_Coors.col(2) = Sph_Coors.col(2) - pt(2);
  return Sph_Coors;
}
