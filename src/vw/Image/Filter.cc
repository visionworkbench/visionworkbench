// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// Filter.cc
///
/// Instantiations of certain filter kernel functions.
///
#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <vw/Image/Filter.tcc>

namespace vw {

  template void generate_gaussian_kernel<float>( std::vector<float>& kernel, double sigma, int32 size );
  template void generate_gaussian_kernel<double>( std::vector<double>& kernel, double sigma, int32 size );

  template void generate_derivative_kernel<float>( std::vector<float>& kernel, int32 deriv, int32 size );
  template void generate_derivative_kernel<double>( std::vector<double>& kernel, int32 deriv, int32 size );

  template void generate_gaussian_derivative_kernel<float>( ImageView<float>& kernel, double sigma1, int32 deriv1, double sigma2, int32 deriv2, double angle, int32 size );
  template void generate_gaussian_derivative_kernel<double>( ImageView<double>& kernel, double sigma1, int32 deriv1, double sigma2, int32 deriv2, double angle, int32 size );

  template void generate_laplacian_of_gaussian_kernel<float>( ImageView<float>& kernel, double sigma, int32 size );
  template void generate_laplacian_of_gaussian_kernel<double>( ImageView<double>& kernel, double sigma, int32 size );

} // namespace vw
