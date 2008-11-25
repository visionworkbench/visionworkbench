// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__

// SURF.h
//
// Classes to perform SURF interest point detection. A lot of SURF
// doesn't fit into the cookie cutter model of IP.
//
// Author:
// Zachary Moratto
// Kansas State University
#include <vw/InterestPoint/SURF.h>

namespace vw {
namespace ip {
  
  /// Creates a interest scale for SURF
  SURFScaleData SURFProcessScale( vw::ImageView<double> const& integral, 
					 unsigned const& octave, unsigned const& scale, 
					 SURFParams const& params ) {
    /// Calculating how big the filter is for current octave and scale
    unsigned first_filter_size_for_octave = 9;
    for (unsigned i = 0; i < octave; ++i)
      first_filter_size_for_octave += 30*(i == 0 ? 1 : 2<<(i-1));
    unsigned filter_size = first_filter_size_for_octave;
    for (unsigned i = 0; i < scale; ++i)
      filter_size += 6*(i == 0 ? 1 : 2<<(i-1));
    unsigned sampling_step = 2<<octave;
    if (scale == 0 && octave != 0) // This for comparisons between octaves (so up'n the resolution)
      sampling_step = 2<<(octave-1);
    unsigned s_temp = floor(float(filter_size)/2.);
    unsigned starting_point = s_temp + (sampling_step - (s_temp & (sampling_step - 1)));

    vw_out(vw::DebugMessage,"interest_point")<< "Octave: " << octave 
					     << " Scale: " << scale 
					     << std::endl;
    vw_out(vw::DebugMessage,"interest_point")<< "\t> filter size:" 
					     << filter_size << std::endl;

    /// Calculating how big the interest image is going to be for our data
    unsigned data_cols = floor(float(params.cols - (starting_point<<1))/float(sampling_step)) + 1;
    unsigned data_rows = floor(float(params.rows - (starting_point<<1))/float(sampling_step)) + 1;
    
    boost::shared_ptr<vw::ImageView<float> > data_image ( new vw::ImageView<float>( data_cols, data_rows ));
    boost::shared_ptr<vw::ImageView<bool> > pol_image ( new vw::ImageView<bool>( data_cols, data_rows ));
    
    vw_out(vw::DebugMessage,"interest_point")<<"\t> interest image size:" << data_image->cols() 
					     << " " << data_image->rows() << std::endl;

    /// Filling in the interest image data by calculating the hessian
    /// for a point. Recording the determinant as the interest operator, and
    /// saving the polarity of the trace
    for (unsigned x = starting_point, ix=0; 
	 ix < data_cols; x+=sampling_step, ix++) {
      for (unsigned y = starting_point, iy=0; 
	   iy < data_rows; y+=sampling_step, iy++) {
	// Getting the determinant
	float Dxy = XYDerivative(integral,x,y,filter_size);
	float Dxx = XSecondDerivative(integral,x,y,filter_size);
	float Dyy = YSecondDerivative(integral,x,y,filter_size);
	data_image->operator()(ix,iy) = Dxx*Dyy - (0.81*Dxy*Dxy);
	pol_image->operator()(ix,iy) = 0 < (Dxx + Dyy);
      }
    }
   
    return SURFScaleData( data_image, pol_image, starting_point, sampling_step);
  }

  /// Finds the maximas in our data
  std::list<InterestPoint> SURFMaximaDetection( std::vector<SURFScaleData> const & scaleData, 
						       SURFParams const& params ) {
    std::list<InterestPoint> ip;

    // A maxima must be larger than it 8 nearest neighbor in it's own
    // layer. A maxima must also be larger that the 18 nearest points in the
    // scale above and below. The first scale for every octave above the
    // first must be over sampled to be able to perform all comparisons with
    // the scale below in the other octave.

    for (unsigned curr_oct = 0; curr_oct < params.octaves; ++curr_oct ) {
      for (unsigned curr_scale = 0; curr_scale < params.scales; ++curr_scale ) {

	// Index for the current octave and scale
	unsigned index = curr_oct*params.scales + curr_scale;
	unsigned sampling_step = 2<<curr_oct;

	if (index == 0)
	  continue;
	if (index == scaleData.size() - 1)
	  continue;

	// Calculating what the current scale is. I'm emulating the
	// numbers that were in the SURF paper
	float scale=1.2;
	for (unsigned o = 0; o < curr_oct; ++o)
	  scale+=6*pow(2,o);
	for (unsigned s = 0; s < curr_scale; ++s)
	  scale+=1.2*pow(2,curr_oct);

	// comparing across all the pixels now
	for (unsigned x = sampling_step; x < params.cols - sampling_step; 
	     x+=sampling_step) {
	  for (unsigned y = sampling_step; y < params.rows - sampling_step; 
	       y+=sampling_step) {

            #define value scaleData[index].determinant(x,y)

	    if (value < params.threshold)
	      continue;
	    
	    // From this scale
	    if (scaleData[index].isLessThan(value,x-sampling_step,y))
	      continue;
	    if (scaleData[index].isLessThan(value,x-sampling_step,y-sampling_step))
	      continue;
	    if (scaleData[index].isLessThan(value,x,y-sampling_step))
	      continue;
	    if (scaleData[index].isLessThan(value,x+sampling_step,y-sampling_step))
	      continue;
	    if (scaleData[index].isLessThan(value,x+sampling_step,y))
	      continue;
	    if (scaleData[index].isLessThan(value,x+sampling_step,y+sampling_step))
	      continue;
	    if (scaleData[index].isLessThan(value,x,y+sampling_step))
	      continue;
	    if (scaleData[index].isLessThan(value,x-sampling_step,y+sampling_step))
	      continue;

	    // From the scale below
	    if (scaleData[index-1].isLessThan(value,x-sampling_step,y))
	      continue;
	    if (scaleData[index-1].isLessThan(value,x-sampling_step,y-sampling_step))
	      continue;
	    if (scaleData[index-1].isLessThan(value,x,y-sampling_step))
	      continue;
	    if (scaleData[index-1].isLessThan(value,x+sampling_step,y-sampling_step))
	      continue;
	    if (scaleData[index-1].isLessThan(value,x+sampling_step,y))
	      continue;
	    if (scaleData[index-1].isLessThan(value,x+sampling_step,y+sampling_step))
	      continue;
	    if (scaleData[index-1].isLessThan(value,x,y+sampling_step))
	      continue;
	    if (scaleData[index-1].isLessThan(value,x-sampling_step,y+sampling_step))
	      continue;
	    if (scaleData[index-1].isLessThan(value,x,y))
	      continue;

	    // From the scale above
	    if (scaleData[index+1].isLessThan(value,x-sampling_step,y))
	      continue;
	    if (scaleData[index+1].isLessThan(value,x-sampling_step,y-sampling_step))
	      continue;
	    if (scaleData[index+1].isLessThan(value,x,y-sampling_step))
	      continue;
	    if (scaleData[index+1].isLessThan(value,x+sampling_step,y-sampling_step))
	      continue;
	    if (scaleData[index+1].isLessThan(value,x+sampling_step,y))
	      continue;
	    if (scaleData[index+1].isLessThan(value,x+sampling_step,y+sampling_step))
	      continue;
	    if (scaleData[index+1].isLessThan(value,x,y+sampling_step))
	      continue;
	    if (scaleData[index+1].isLessThan(value,x-sampling_step,y+sampling_step))
	      continue;
	    if (scaleData[index+1].isLessThan(value,x,y))
	      continue;

	    #undef value

	    // Are you still with me?
            ip.push_back(vw::ip::InterestPoint(x,y,scale,
					       scaleData[index].determinant(x,y),
					       0.0, scaleData[index].polarity(x,y),
					       index));

	  }
	}

      }
    }

    return ip;
  }

  /// Write SURF debug interest images
  /// This method is for debug purposes only. It will write out the
  /// interest images that are used for interest detect for all octaves and
  /// scales
  void SURFWriteDebugImages( std::vector<SURFScaleData> const& scaleData,
			     SURFParams const& params ) {
    if (scaleData.size() == 0)
      return;
    for ( unsigned octave_num = 0; octave_num < params.octaves; ++octave_num) {
	for ( unsigned scale_num = 0; scale_num < params.scales; ++scale_num) {
	  unsigned index = octave_num*params.scales + scale_num;
	  unsigned sampling_step = 2<<octave_num;
	  if (scale_num == 0 && octave_num != 0)
	    sampling_step = 2<<(octave_num - 1);
	  unsigned data_cols = floor(float(params.cols)/float(sampling_step)) + 1;
	  unsigned data_rows = floor(float(params.rows)/float(sampling_step)) + 1;

	  vw::ImageView<float> debug( data_cols, data_rows);
	  for ( unsigned sx = 0, dx = 0; sx < params.cols; sx+=sampling_step,dx++)
	    for ( unsigned sy = 0, dy = 0; sy < params.rows; sy+=sampling_step, dy++)
	      debug(dx,dy) = scaleData[index].determinant(sx,sy);
	
	  // Normalizing
	  float min = debug(0,0);
	  for (vw::ImageView<float>::iterator iter = debug.begin();
	     iter != debug.end(); ++iter)
	    if (*iter < min) min = *iter;
	  std::cout << "\t> min: " << min;
	  debug -= min;
	  float max = debug(0,0);
	  for (vw::ImageView<float>::iterator iter = debug.begin();
	       iter != debug.end(); ++iter)
	    if (*iter > max) max = *iter;
	  std::cout << " max: " << max << std::endl;
	  debug /= max;

	  // Name
	  std::ostringstream ostr;
	  ostr << "dbg_iimage_" << octave_num << "_" << scale_num << ".png";
	  write_image( ostr.str(), debug );
	  
	}
      }
  }

  /// SURF Subpixel refinement
  /// This gets floating point precision for the interest points, will
  /// also throw at some points if they move too far.
  void SURFSubpixelRefinement( std::list<vw::ip::InterestPoint>& ip,
			       std::vector<SURFScaleData> const& scaleData,
			       SURFParams const& params ) {
    // Calculating scaleSizes before hand, this is important for when
    // it comes time to interpolate what the actual scale is. SURF has lotsa
    // weird scale values. Here's what to expect:
    // (1.2, 2.4, 3.6, 4.8, 7.2, 9.6, 12.0, 14.4, 19.2 ... and so on.
    std::vector<float> scaleSize(params.octaves*params.scales);
    scaleSize[0] = 1.2;
    for (unsigned curr_oct = 0; curr_oct < params.octaves - 1; ++curr_oct )
      scaleSize[(curr_oct+1)*params.scales] = scaleSize[curr_oct*params.scales] + 6*pow(2,curr_oct);
    for (unsigned curr_oct = 0; curr_oct < params.octaves; ++curr_oct ) {
      for (unsigned curr_scale = 1; curr_scale < params.scales; ++curr_scale ) {
	unsigned index = curr_oct*params.scales + curr_scale;
	scaleSize[index]=scaleSize[index-1]+1.2*pow(2,curr_oct);
      }
    }

    // Performing refinement on all pixels
    for (std::list<vw::ip::InterestPoint>::iterator point = ip.begin();
	 point != ip.end(); ++point) {

      signed dx = 0, dy = 0;
      vw::Vector3 result;
      unsigned iter;

      for (iter = 0; iter < 5; ++iter) {
	// Apply changes
	(*point).ix += dx;
	(*point).iy += dy;

	// Calculate what the sample step is (I'm going to keep this in
	// terms of the main image)
	unsigned octave = (unsigned)floor(double((*point).index)/double(params.scales));
	unsigned sampling_step = 2<<octave;
	unsigned filter_size = 9;
	for (unsigned i = 0; i < octave; ++i)
	  filter_size += 30*(i == 0 ? 1 : 2<<(i-1));
	unsigned s_temp = floor(float(filter_size)/2.);
	signed starting_point = s_temp + (sampling_step - (s_temp & (sampling_step - 1)));

	if ((*point).ix < starting_point || (*point).ix > int(params.cols) - starting_point) {
	  vw_out(DebugMessage, "interest_point") << "SURFSubpixel Refinement: ix is out of bounds.\n";
	  continue;
	}
	if ((*point).iy < starting_point || (*point).iy > int(params.rows) - starting_point) {
	  vw_out(DebugMessage, "interest_point") << "SURFSubpixel Refinement: iy is out of bounds.\n";
	  continue;
	}

	// Finding gradients
	vw::Vector3 B = SURFGradient3D(scaleData, (*point).ix, (*point).iy,
				       (*point).index, sampling_step );

	// Finding second derivatives
	vw::Matrix3x3 A = SURFHessian3D(scaleData, (*point).ix, (*point).iy,
					(*point).index, sampling_step );

	result = vw::math::solve(A,B);

	dx = dy = 0;
	if ( (result(0) > float(sampling_step)/2 || result(0) < -float(sampling_step)/2) && ( (*point).ix+round(result(0))) < params.cols - starting_point
	     && ((*point).ix+round(result(0))) > starting_point)
	  dx = result(0) > 0 ? sampling_step : -sampling_step;

        if ( (result(1) > float(sampling_step)/2 || result(1) < -float(sampling_step)/2) && ( (*point).iy+round(result(1))) < params.rows - starting_point 
	     && ((*point).iy+round(result(1))) > starting_point)
	  dy = result(1) > 0 ? sampling_step : -sampling_step;
      
	if (dx == 0 && dy == 0) break;
	
      } //End of iter

      // Deciding if to remove for too many iterations
      if (iter >= 4 ) {
	vw_out(vw::DebugMessage,"interest_point")<<"removing point\tx:" << (*point).ix << "\ty:" << (*point).iy << "\ts:" << (*point).scale << std::endl;
	point = ip.erase(point);
	point--;
	continue;
      }

      // Calculating the refinement
      {
	(*point).x = float((*point).ix) + result(0);
	(*point).y = float((*point).iy) + result(1);
	// Not trying to make you think or anything
	(*point).scale += (result(2)>0 ? result(2) : -result(2))*(scaleSize[result(2) > 0 ? (*point).index+1 : (*point).index-1]-scaleSize[(*point).index]);
      }
    } // end of for
  }

  // gradient 3D
  // - scaleData = Contains the interest data
  // - ix        = x location to evaluate at
  // - iy        = y location to evaluate at
  // - index     = index to evaluate at, representative of octave & scale
  // - step      = the step size that is for the interest data
  // This calculates the gradient, hmm... it negative
  Vector3& SURFGradient3D( std::vector<SURFScaleData> const& scaleData,
			       int const& ix, int const& iy, 
			       unsigned const& index, unsigned const& step ) {

    vw::Vector3 B;
    B(0) = -0.5*(scaleData[index].determinant(ix+step,iy) -
		 scaleData[index].determinant(ix-step,iy));
    B(1) = -0.5*(scaleData[index].determinant(ix,iy+step) -
		 scaleData[index].determinant(ix,iy-step));
    B(2) = -0.5*(scaleData[index+1].determinant(ix,iy) -
		 scaleData[index-1].determinant(ix,iy) );
    
    return B;
  }

  
  // hessian 3D
  // - scaleData = contains the interest data
  // - ix        = x location to evaluate at
  // - iy        = y location to evaluate at
  // - index     = index to evaluate at, representative of octave & scale
  // - step      = the step size that is for the interest data
  // This calculates the hessian matrix
  Matrix3x3& SURFHessian3D( std::vector<SURFScaleData> const& scaleData, 
			    int const& ix, int const& iy, 
			    unsigned const& index, unsigned const& step ) {

    vw::Matrix3x3 A;
    // Dxx
    A(0,0) = scaleData[index].determinant(ix+step,iy) +
      scaleData[index].determinant(ix-step,iy) -
      2.0*scaleData[index].determinant(ix,iy);
    // Dyy
    A(1,1) = scaleData[index].determinant(ix,iy+step) +
      scaleData[index].determinant(ix,iy-step) -
      2.0*scaleData[index].determinant(ix,iy);
    // Dzz
    A(2,2) = scaleData[index+1].determinant(ix,iy) +
      scaleData[index-1].determinant(ix,iy) -
      2.0*scaleData[index].determinant(ix,iy);
    
    // keep'n it invertable yo :)
    A(0,0) += 1e-29;
    A(1,1) += 1e-29;
    A(2,2) += 1e-29;

    // Dxy
    A(0,1) = A(1,0) = 0.25*(scaleData[index].determinant(ix+step,iy+step) +
			    scaleData[index].determinant(ix-step,iy-step) -
			    scaleData[index].determinant(ix+step,iy-step) -
			    scaleData[index].determinant(ix-step,iy+step));
	
    // Dxs
    A(0,2) = A(2,0) = 0.25*(scaleData[index+1].determinant(ix+step,iy) +
			    scaleData[index-1].determinant(ix-step,iy) -
			    scaleData[index-1].determinant(ix+step,iy) -
			    scaleData[index+1].determinant(ix-step,iy) );
    // Dys
    A(1,2) = A(2,1) = 0.25*(scaleData[index+1].determinant(ix,iy+step) +
			    scaleData[index-1].determinant(ix,iy-step) -
			    scaleData[index-1].determinant(ix,iy+step) -
			    scaleData[index+1].determinant(ix,iy-step) );

    return A;
  }

  // SURF Orientation
  // - integral  = Integral used for calculations
  // - ix        = x location to evaluate at
  // - iy        = y location to evaluate at
  // - scale     = scale at which the feature was found like 7.2

  // This calculates the orientationin radians according to the SURF
  // paper. Well, kinda, this implementation is a little lazy.
  float SURFOrientation( vw::ImageView<double> const& integral, 
			 int const& ix, int const& iy,
			 float const& scale ) {
    int iscale = round(scale);

    std::vector<float> h_response(49);
    std::vector<float> v_response(49);
    std::vector<float> angle(49);
    float ori = 0;

    // Calculating all the responses
    for (signed idx = 0, i = ix-3*iscale; 
	 i <=(ix+3*iscale); i+= iscale ) {
      for (signed j = iy-3*iscale; j<=(iy+3*iscale); j+=iscale) {
	// Calculating the haar wavelet response
	if ( (i - 2*iscale < 0 || j - 2*iscale < 0) || 
	     (i + 2*iscale >= integral.cols() || j + 2*iscale >= integral.rows() ) ) {
	  continue;
	}
	h_response[idx] = HHaarWavelet( integral, i, j, iscale*4 );
	v_response[idx] = VHaarWavelet( integral, i, j, iscale*4 );
	angle[idx] = atan2(v_response[idx],h_response[idx]);
	idx++;
      }
    }

    // Fitting a slice to find the main reponse
    const float pi_6 = 3.14159/6.0;
    const float two_pi = 3.14159*2.0;
    float sumx, sumy, mod, greatest_mod;
    greatest_mod = 0;
    for (float a = 0; a < two_pi; a += 0.01) {
      sumx = sumy = 0;
      for (int idx = 0; idx < 49; idx++ ) {
	// Is it in my slice
	if ( (angle[idx] > a - pi_6 && angle[idx] < a + pi_6 ) ||
	     (angle[idx] + two_pi > a - pi_6 && angle[idx] + two_pi < a + pi_6 ) ||
	     (angle[idx] - two_pi > a - pi_6 && angle[idx] - two_pi < a + pi_6 ) ) {
	  
	  float diff_a = fabs( angle[idx] - a );
	  if ( diff_a > 3.14159 )
	    diff_a = fabs( diff_a - two_pi );

	  float weight = -6*diff_a/3.14159 + 1;

	  sumx += weight*h_response[idx];
	  sumy += weight*v_response[idx];
	}
      }
      
      mod = sumx*sumx + sumy*sumy;
      if ( mod > greatest_mod ) {
	greatest_mod = mod;
	ori = a;
      }
    }

    return ori;
  }


}} // namespace vw::ip
