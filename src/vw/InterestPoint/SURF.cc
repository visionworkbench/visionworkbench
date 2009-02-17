// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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
#include <vw/Image/ImageMath.h>
#include <vw/Math/LinearAlgebra.h>

namespace vw {
namespace ip {
  
  /// Creates a interest scale for SURF
  SURFScaleData SURFProcessScale( vw::ImageView<double> const& integral, 
				  unsigned const& filter_size, 
				  unsigned const& sampling_step,
				  SURFParams const& params ) {
    /// Calculating how big the filter is for current octave and scale
    unsigned s_temp = floor(float(filter_size)/2.);
    unsigned starting_point = s_temp + (sampling_step - (s_temp & (sampling_step - 1)));

    vw_out(vw::DebugMessage,"interest_point")<< "Building Interest Image: ";
    vw_out(vw::DebugMessage,"interest_point")<< "\t> filter size: " << filter_size 
					     << "\n\t> sampling step: " << sampling_step
					     << std::endl;

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
   
    return SURFScaleData( data_image, pol_image, filter_size, sampling_step);
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
      for (unsigned curr_scale = 1; curr_scale < (params.scales - 1 ); ++curr_scale ) {

	// Index for the current octave and scale
	unsigned index0 = params.convToIndex( curr_oct, curr_scale );
	unsigned indexDown = params.convToIndex( curr_oct, curr_scale - 1 );
	unsigned indexUp = params.convToIndex( curr_oct, curr_scale + 1 );
	unsigned sampling_step = 2<<curr_oct;

	float scale = scaleData[index0].filter_size();
	
	// comparing across all the pixels now
	for (unsigned x = sampling_step; x < params.cols - sampling_step; 
	     x+=sampling_step) {
	  for (unsigned y = sampling_step; y < params.rows - sampling_step; 
	       y+=sampling_step) {

            #define value scaleData[index0].determinant(x,y)

	    if (value < params.threshold)
	      continue;
	    
	    // From this scale
	    if (scaleData[index0].isLessThan(value,x-sampling_step,y))
	      continue;
	    if (scaleData[index0].isLessThan(value,x-sampling_step,y-sampling_step))
	      continue;
	    if (scaleData[index0].isLessThan(value,x,y-sampling_step))
	      continue;
	    if (scaleData[index0].isLessThan(value,x+sampling_step,y-sampling_step))
	      continue;
	    if (scaleData[index0].isLessThan(value,x+sampling_step,y))
	      continue;
	    if (scaleData[index0].isLessThan(value,x+sampling_step,y+sampling_step))
	      continue;
	    if (scaleData[index0].isLessThan(value,x,y+sampling_step))
	      continue;
	    if (scaleData[index0].isLessThan(value,x-sampling_step,y+sampling_step))
	      continue;

	    // From the scale below
	    if (scaleData[indexDown].isLessThan(value,x-sampling_step,y))
	      continue;
	    if (scaleData[indexDown].isLessThan(value,x-sampling_step,y-sampling_step))
	      continue;
	    if (scaleData[indexDown].isLessThan(value,x,y-sampling_step))
	      continue;
	    if (scaleData[indexDown].isLessThan(value,x+sampling_step,y-sampling_step))
	      continue;
	    if (scaleData[indexDown].isLessThan(value,x+sampling_step,y))
	      continue;
	    if (scaleData[indexDown].isLessThan(value,x+sampling_step,y+sampling_step))
	      continue;
	    if (scaleData[indexDown].isLessThan(value,x,y+sampling_step))
	      continue;
	    if (scaleData[indexDown].isLessThan(value,x-sampling_step,y+sampling_step))
	      continue;
	    if (scaleData[indexDown].isLessThan(value,x,y))
	      continue;

	    // From the scale above
	    if (scaleData[indexUp].isLessThan(value,x-sampling_step,y))
	      continue;
	    if (scaleData[indexUp].isLessThan(value,x-sampling_step,y-sampling_step))
	      continue;
	    if (scaleData[indexUp].isLessThan(value,x,y-sampling_step))
	      continue;
	    if (scaleData[indexUp].isLessThan(value,x+sampling_step,y-sampling_step))
	      continue;
	    if (scaleData[indexUp].isLessThan(value,x+sampling_step,y))
	      continue;
	    if (scaleData[indexUp].isLessThan(value,x+sampling_step,y+sampling_step))
	      continue;
	    if (scaleData[indexUp].isLessThan(value,x,y+sampling_step))
	      continue;
	    if (scaleData[indexUp].isLessThan(value,x-sampling_step,y+sampling_step))
	      continue;
	    if (scaleData[indexUp].isLessThan(value,x,y))
	      continue;

	    #undef value

	    // Are you still with me?
            ip.push_back(vw::ip::InterestPoint(x,y,scale,
					       scaleData[index0].determinant(x,y),
					       0.0, scaleData[index0].polarity(x,y),
					       curr_oct, curr_scale));

	    // Note at this point the Interest Point's scale actually
	    // corresponds to the filter size. (it's not that
	    // important)

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


    for ( unsigned octave_num = 0; octave_num < params.octaves; ++octave_num) {

      for ( unsigned scale_num = 0; scale_num < params.scales; ++scale_num) {
	unsigned index = params.convToIndex( octave_num, scale_num );
	unsigned sampling_step = 2 << octave_num;
	
	unsigned data_cols = floor(float(params.cols)/float(sampling_step)) + 1;
	unsigned data_rows = floor(float(params.rows)/float(sampling_step)) + 1;
	
	vw::ImageView<float> debug( data_cols, data_rows);
	for ( int sx = 0, dx = 0; sx < params.cols; sx+=sampling_step,dx++)
	  for ( int sy = 0, dy = 0; sy < params.rows; sy+=sampling_step, dy++)
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
	
      }  // End scale iteration
    } // End octave iteration
  }

  /// SURF Subpixel refinement
  /// This gets floating point precision for the interest points, will
  /// also throw at some points if they move too far.
  void SURFSubpixelRefinement( std::list<vw::ip::InterestPoint>& ip,
			       std::vector<SURFScaleData> const& scaleData,
			       SURFParams const& params ) {

    // Performing refinement on all pixels
    for (std::list<vw::ip::InterestPoint>::iterator point = ip.begin();
	 point != ip.end(); ++point) {

      signed dx = 0, dy = 0, ds = 0;
      vw::Vector3 result;
      unsigned iter;
      unsigned sampling_step = 2<<(*point).octave;
      
      for (iter = 0; iter < 5; ++iter) {
	
	// Apply changes
	(*point).ix += dx;
	(*point).iy += dy;
	(*point).scale_lvl += ds;
	dx = dy = ds = 0;

	unsigned index0 = params.convToIndex( (*point).octave,
					      (*point).scale_lvl );

	// Has this point left the boundaries?  Note: isValid does
	//    double duty. It can determine if the point has left the
	//    edge of interest image
	if ( !scaleData[index0].isValid( (*point).ix, (*point).iy ) ||
	     (*point).scale_lvl < 1 ||
	     (*point).scale_lvl > params.scales - 2 ) {

	  // Axe this point then
	  iter = 5;
	  break;
	}

	// Finding gradients
	vw::Vector3 B = SURFGradient3D( scaleData, (*point), params );

	// Finding second derivatives
	vw::Matrix3x3 A = SURFHessian3D( scaleData, (*point), params );

	result = vw::math::solve(A,B);

	if (result(0) > 0.5 || result(0) < -0.5)
	  dx = result(0) > 0 ? sampling_step : -sampling_step;
        if (result(1) > 0.5 || result(1) < -0.5)
	  dy = result(1) > 0 ? sampling_step : -sampling_step;
	if (result(2) > 0.5 || result(2) < -0.5 )
	  ds = result(2) > 0 ? 1 : - 1;
      
	if (dx == 0 && dy == 0 && ds == 0) break;

      } //End of iter

      // Deciding if to remove for too many iterations
      if (iter > 4 ) {
	point = ip.erase(point);
	point--;
	continue;
      }
      
      // Calculating the refinement
      {
	(*point).x = float((*point).ix) + result(0)*sampling_step;
	(*point).y = float((*point).iy) + result(1)*sampling_step;
	(*point).ix += round(result(0)*sampling_step);
	(*point).iy += round(result(1)*sampling_step);

	// Scale Interpolation
	//   At this time IP's scale becomes actually scale size (which is the gaussian variance of the filter)
	unsigned index0 = params.convToIndex( (*point).octave, (*point).scale_lvl );
	unsigned index1 = params.convToIndex( (*point).octave,
					      (*point).scale_lvl + ( result(2)>0 ? 1 : -1 ) );

	(*point).scale = params.calcScale( float(scaleData[index0].filter_size()) + 
					   fabs( result(2) )*( float(scaleData[index1].filter_size()) -
							       float(scaleData[index0].filter_size()) ) );
					   
      }

    } // end of for
  }

  // gradient 3D
  // - scaleData = Contains the interest data
  // - ix        = x location to evaluate at
  // - iy        = y location to evaluate at
  // - index     = index to evaluate at, representative of octave & scale
  // - step      = the step size that is for the interest data
  // This calculates the gradient, hmm... its negative
  Vector3 SURFGradient3D( std::vector<SURFScaleData> const& scaleData,
			  vw::ip::InterestPoint const& ip, 
			  SURFParams const& params ) {

    int step = 0x2<<ip.octave;
    int index0 = params.convToIndex( ip.octave, ip.scale_lvl );
    int indexDown = params.convToIndex( ip.octave, ip.scale_lvl - 1);
    int indexUp = params.convToIndex( ip.octave, ip.scale_lvl + 1);

    vw::Vector3 B;
    B(0) = -0.5*(scaleData[index0].determinant(ip.ix+step,ip.iy) -
		 scaleData[index0].determinant(ip.ix-step,ip.iy));
    B(1) = -0.5*(scaleData[index0].determinant(ip.ix,ip.iy+step) -
		 scaleData[index0].determinant(ip.ix,ip.iy-step));
    B(2) = -0.5*(scaleData[indexUp].determinant(ip.ix,ip.iy) -
		 scaleData[indexDown].determinant(ip.ix,ip.iy) );
    
    return B;
  }

  
  // hessian 3D
  // - scaleData = contains the interest data
  // - ix        = x location to evaluate at
  // - iy        = y location to evaluate at
  // - index     = index to evaluate at, representative of octave & scale
  // - step      = the step size that is for the interest data
  // This calculates the hessian matrix
  Matrix3x3 SURFHessian3D( std::vector<SURFScaleData> const& scaleData, 
			   vw::ip::InterestPoint const& ip, 
			   SURFParams const& params ) {

    int step = 0x2<<ip.octave;
    int index0 = params.convToIndex( ip.octave, ip.scale_lvl );
    int indexDown = params.convToIndex( ip.octave, ip.scale_lvl - 1);
    int indexUp = params.convToIndex( ip.octave, ip.scale_lvl + 1 );

    vw::Matrix3x3 A;
    // Dxx
    A(0,0) = scaleData[index0].determinant(ip.ix+step, ip.iy) +
      scaleData[index0].determinant(ip.ix-step,ip.iy) -
      2.0*scaleData[index0].determinant(ip.ix,ip.iy);
    // Dyy
    A(1,1) = scaleData[index0].determinant(ip.ix,ip.iy+step) +
      scaleData[index0].determinant(ip.ix,ip.iy-step) -
      2.0*scaleData[index0].determinant(ip.ix,ip.iy);
    // Dzz
    A(2,2) = scaleData[indexUp].determinant(ip.ix,ip.iy) +
      scaleData[indexDown].determinant(ip.ix,ip.iy) -
      2.0*scaleData[index0].determinant(ip.ix,ip.iy);
    
    // Keep'n it invertable
    A(0,0) += 1e-30;
    A(1,1) += 1e-30;
    A(2,2) += 1e-30;

    // Dxy
    A(0,1) = A(1,0) = 0.25*(scaleData[index0].determinant(ip.ix+step,ip.iy+step) +
			    scaleData[index0].determinant(ip.ix-step,ip.iy-step) -
			    scaleData[index0].determinant(ip.ix+step,ip.iy-step) -
			    scaleData[index0].determinant(ip.ix-step,ip.iy+step));
	
    // Dxs
    A(0,2) = A(2,0) = 0.25*(scaleData[indexUp].determinant(ip.ix+step,ip.iy) +
			    scaleData[indexDown].determinant(ip.ix-step,ip.iy) -
			    scaleData[indexDown].determinant(ip.ix+step,ip.iy) -
			    scaleData[indexUp].determinant(ip.ix-step,ip.iy) );
    // Dys
    A(1,2) = A(2,1) = 0.25*(scaleData[indexUp].determinant(ip.ix,ip.iy+step) +
			    scaleData[indexDown].determinant(ip.ix,ip.iy-step) -
			    scaleData[indexDown].determinant(ip.ix,ip.iy+step) -
			    scaleData[indexUp].determinant(ip.ix,ip.iy-step) );

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
			 float const& x, float const& y,
			 float const& scale) {

    std::vector<float> h_response(169);
    std::vector<float> v_response(169);
    std::vector<float> angle(169);
    int measures = 0;

    for ( int i = -6; i <= 6; i++ ) {
      for ( int j = -6; j <= 6; j++ ) {
	// Is this this within the radius of 6
	if ( (i*i + j*j) > 36 )
	  continue;

	// Is this still on the image?
	if ( (floor(x + i*scale - round(scale*4)/2) < 0 ) || 
	     (floor(y + j*scale - round(scale*4)/2) < 0 ) ||
	     (ceil(x + i*scale + round(scale*4)/2) + 1 >= integral.cols()) || 
	     (ceil(y + j*scale + round(scale*4)/2) + 1 >= integral.rows()) )
	  continue;

	float distance_2 = i*i + j*j;
	float weight = exp(-distance_2/8)/5.0133;

      	h_response[measures] =
	  weight*HHaarWavelet( integral, 
			       int(round(x + i*scale)), 
			       int(round(y + j*scale)), 
			       round(scale*4) );
	v_response[measures] =
	  weight*VHaarWavelet( integral, 
			       int(round(x + i*scale)), 
			       int(round(y + j*scale)), 
			       round(scale*4) );
	angle[measures] = atan2( v_response[measures],
				 h_response[measures] );
	measures++;
      }
    }

    // Fitting a slice to find the response
    const float pi_6 = 3.14159/6.0;
    const float two_pi = 3.14159*2.0;
    float sumx, sumy;
    float mod, second_mod=0, greatest_mod=0;
    float second_ori=0, greatest_ori=0;
    for ( float a = 0; a < two_pi; a+= 0.5) {
      sumx = sumy = 0;
      for ( int idx = 0; idx < measures; idx++ ) {
	// Is it in my slice
	if ( (angle[idx] > a - pi_6 && angle[idx] < a + pi_6 ) ||
	     (angle[idx] + two_pi > a - pi_6 && angle[idx] + two_pi < a + pi_6 ) ||
	     (angle[idx] - two_pi > a - pi_6 && angle[idx] - two_pi < a + pi_6 ) ) {

	  sumx += h_response[idx];
	  sumy += v_response[idx];
	}
      }

      mod = sumx*sumx + sumy*sumy;
      if ( mod > greatest_mod ) {
	// Storing second greatest
	second_mod = greatest_mod;
	second_ori = greatest_ori;

	greatest_mod = mod;
	greatest_ori = atan2(sumy,sumx);
      }
    }

    return greatest_ori;
  }

  // SURF Descriptor
  Vector<float> SURFDescriptor( ImageView<double> const& integral,
				Matrix<float,20,20> const& gaussian,
				InterestPoint const& ip,
				bool extended = false ) {

    Vector<float> result(64);
    if ( extended )
      result.set_size(128);
    Matrix<float,20,20> h_response;
    Matrix<float,20,20> v_response;
    float scaling = 1.0/ip.scale;

    // Building Transforms
    TransformRef txform( compose(ResampleTransform(scaling,scaling),
				 RotateTransform(-ip.orientation),
				 TranslateTransform(-ip.x,-ip.y) ) );
    TransformRef rotate_only( RotateTransform(-ip.orientation ) );

    // Wrapping integral
    InterpolationView<EdgeExtensionView<ImageView<double>, ConstantEdgeExtension>, BilinearInterpolation> wrapped_integral = interpolate( integral, BilinearInterpolation() );

    // Building responses
    for ( char x = 0; x < 20; x++ ) {
      for ( char y = 0; y < 20; y++ ) {

	// Sampling Point's location
	Vector2 location = txform.reverse( Vector2(float(x) - 9.5,
						   float(y) - 9.5) );

	// Is this response within the image?
	if ( location.x() + ip.scale+1 < integral.cols() &&
	     location.x() - ip.scale >= 0 &&
	     location.y() + ip.scale+1 < integral.rows() &&
	     location.y() - ip.scale >= 0 ) {

	  Vector2 response;
	  
	  response[0] = HHaarWavelet( wrapped_integral,
				      location.x(), location.y(),
				      2*ip.scale );
	  response[1] = VHaarWavelet( wrapped_integral,
				      location.x(), location.y(),
				      2*ip.scale );

	  // Rotating and weighting the response
	  response = rotate_only.forward(response);
	  response *= gaussian(x,y);
	  h_response(x,y) = response.x();
	  v_response(x,y) = response.y();
	  
	} else {
	  // well, boo..
	  h_response(x,y) = 0;
	  v_response(x,y) = 0;
	}
      }
    }

    // Building the descriptor
    for ( char x = 0; x < 4; x++ ) {
      for ( char y = 0; y < 4; y++ ) { 
	int dest;
	if (extended)
	  dest = 32*x + 8*y;
	else
	  dest = 16*x + 4*y;

	// Summing responses
	for ( char ix = 0; ix < 5; ix++ ) {
	  for ( char iy = 0; iy < 5; iy++ ) {
	    char sx = x*5+ix;
	    char sy = y*5+iy;
	    
	    if (!extended) {
	      // SURF 64
	      // Dx
	      result(dest+0) += h_response(sx, sy);
	      // |Dx|
	      result(dest+1) += fabs(h_response(sx, sy));
	      // Dy
	      result(dest+2) += v_response(sx, sy);
	      // |Dy|
	      result(dest+3) += fabs(v_response(sx, sy));
	    } else {
	      // SURF 128
	      if ( v_response(sx, sy) >= 0 ) {
		// Dx
		result(dest+0) += h_response(sx, sy);
		// |Dx|
		result(dest+1) += fabs(h_response(sx, sy));
	      } else {
		// Dx
		result(dest+2) += h_response(sx, sy);
		// |Dx|
		result(dest+3) += fabs(h_response(sx, sy));
	      }
	      if ( h_response(sx, sy) >= 0 ) {
		// Dy
		result(dest+4) += v_response(sx, sy);
		// |Dy|
		result(dest+5) += fabs(v_response(sx, sy));
	      } else {
		// Dy
		result(dest+6) += v_response(sx, sy);
		// |Dy|
		result(dest+7) += fabs(v_response(sx, sy));
	      }
	    }
	  }
	}
      }
    }

    return normalize(result);

  }

  // MSURF Descriptor
  Vector<float> MSURFDescriptor( ImageView<double> const& integral,
				 Matrix<float,4,4> const& overall_gaus,
				 Matrix<float,9,9> const& sub_region_gaus,
				 InterestPoint const& ip,
				 bool extended ) {

    Vector<float> result(64);
    if ( extended )
      result.set_size(128);
    Matrix<float,24,24> h_response;
    Matrix<float,24,24> v_response;
    float scaling = 1.0/ip.scale;

    // Building Transforms
    TransformRef txform( compose(ResampleTransform(scaling,scaling),
				 RotateTransform(-ip.orientation),
				 TranslateTransform(-ip.x,-ip.y) ) );
    TransformRef rotate_only( RotateTransform(-ip.orientation ) );

    // Wrapping integral
    InterpolationView<EdgeExtensionView<ImageView<double>, ConstantEdgeExtension>, BilinearInterpolation> wrapped_integral = interpolate( integral, BilinearInterpolation() );

    // Building responses
    for ( char x = 0; x < 24; x++ ) {
      for ( char y = 0; y < 24; y++ ) {

	// Sampling Point's Location
	Vector2 location = txform.reverse( Vector2(float(x) - 11.5,
						   float(y) - 11.5) );

	// Is this response within the image?
	if ( location.x() + ip.scale+1 < integral.cols() &&
	     location.x() - ip.scale >= 0 &&
	     location.y() + ip.scale+1 < integral.rows() &&
	     location.y() - ip.scale >= 0 ) {

	  Vector2 response;
	  
	  response[0] = HHaarWavelet( wrapped_integral,
				      location.x(), location.y(),
				      2*ip.scale );
	  response[1] = VHaarWavelet( wrapped_integral,
				      location.x(), location.y(),
				      2*ip.scale );

	  // Rotating the response
	  response = rotate_only.forward(response);
	  h_response(x,y) = response.x();
	  v_response(x,y) = response.y();
	  
	} else {
	  // well, boo..
	  h_response(x,y) = 0;
	  v_response(x,y) = 0;
	}
      }
    }

    // Building the descriptor
    for ( char x = 0; x < 4; x++ ) {
      for ( char y = 0; y < 4; y++ ) {
	int dest;
	if (extended)
	  dest = 32*x+8*y;
	else
	  dest = 16*x+4*y;
	
	// Summing responses
	for ( char ix = 0; ix < 9; ix++ ) {
	  for ( char iy = 0; iy < 9; iy++ ) {
	    char sx = x*5+ix;
	    char sy = y*5+iy;
	    float mp = overall_gaus(x,y) * sub_region_gaus(ix,iy);

	    if (!extended) {
	      // MSURF 64
	      // Dx
	      result(dest+0) += mp * h_response(sx, sy);
	      // |Dx|
	      result(dest+1) += mp * fabs(h_response(sx, sy));
	      // Dy
	      result(dest+2) += mp * v_response(sx, sy);
	      // |Dy|
	      result(dest+3) += mp * fabs(v_response(sx, sy));
	    } else {
	      // MSURF 128
	      if ( v_response(sx, sy) >= 0 ) {
		// Dx
		result(dest+0) += mp * h_response(sx, sy);
		// |Dx|
		result(dest+1) += mp * fabs(h_response(sx, sy));
	      } else {
		// Dx
		result(dest+2) += mp * h_response(sx, sy);
		// |Dx|
		result(dest+3) += mp * fabs(h_response(sx, sy));
	      }
	      if ( h_response(sx, sy) >= 0 ) {
		// Dy
		result(dest+4) += mp * v_response(sx, sy);
		// |Dy|
		result(dest+5) += mp * fabs(v_response(sx, sy));
	      } else {
		// Dy
		result(dest+6) += mp * v_response(sx, sy);
		// |Dy|
		result(dest+7) += mp * fabs(v_response(sx, sy));
	      }
	    }
	  }
	}
      }
    }

    return normalize(result);
  }
}} // namespace vw::ip
