// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <cxxtest/TestSuite.h>
#include <vw/InterestPoint/IntegralImage.h>
#include <vw/InterestPoint/SURF.h>
#include <vw/Image.h>
#include <vw/FileIO.h>

using namespace vw;
using namespace vw::ip;

class TestSURF : public CxxTest::TestSuite
{
 public:
  
  void test_subpixel_refinement() {

    // Params
    SURFParams params( .001, 1, 3 );
    params.cols = 40;
    params.rows = 40;

    std::list<Vector2> test_points;
    test_points.push_back(Vector2(22.25,15.625));
    test_points.push_back(Vector2(22.83,15.625));
    test_points.push_back(Vector2(22.83,15.17));
    test_points.push_back(Vector2(22.25,15.2));
    test_points.push_back(Vector2(13.1,10.6));
    test_points.push_back(Vector2(12.3,11.7));
    test_points.push_back(Vector2(24.6,26.7));
    test_points.push_back(Vector2(28.4,32.6));

    for (std::list<Vector2>::iterator cur = test_points.begin();
	 cur != test_points.end();
	 cur++ ) {

      // Building fictitious data
      std::vector<SURFScaleData> interest_scales;
      interest_scales.resize(3);
      Vector2 center = (*cur);
      for (unsigned i = 0; i < 3; i++ ) {
	
	// Parameters
	int sampling_step = 2;
	int filter_size = 3+i*3;
	unsigned s_temp = floor(float(filter_size)/2.);
	unsigned starting_point = s_temp + (sampling_step - (s_temp & (sampling_step - 1 )));
	
	unsigned data_cols = floor(float(params.cols - (starting_point<<1))/float(sampling_step))+1;
	unsigned data_rows = floor(float(params.rows - (starting_point<<1))/float(sampling_step))+1;
	
	boost::shared_ptr<vw::ImageView<float> > data_image ( new vw::ImageView<float>( data_cols, data_rows ));
	boost::shared_ptr<vw::ImageView<bool> > pol_image ( new vw::ImageView<bool>( data_cols, data_rows ));

	// Filling in the data_image with a circle gradient
	for (unsigned x = starting_point, ix = 0;
	     ix < data_cols; x+=sampling_step, ix++) {
	  for (unsigned y = starting_point, iy = 0;
	       iy < data_rows; y+=sampling_step, iy++) {
	    Vector2 loc(x,y);
	    loc -= center;
	    data_image->operator()(ix,iy) = ( i==1?20:15 ) - norm_2(loc);
	    pol_image->operator()(ix,iy) = false;
	  }
	}
	
	interest_scales[i] = SURFScaleData( data_image, pol_image, 
					    filter_size, sampling_step );
	
      }
      
      // Creating Interest Points at various locations to the center.
      InterestPointList ip;
      ip.push_back( InterestPoint(22,16,params.calcScale(6),
				  interest_scales[1].determinant(22,16),
				  0.0, interest_scales[1].polarity(22,16),
				  0,1) );
      ip.push_back( InterestPoint(18,16,params.calcScale(6),
				  interest_scales[1].determinant(18,16),
				  0.0, interest_scales[1].polarity(18,16),
				  0,1) );
      ip.push_back( InterestPoint(22,20,params.calcScale(6),
				  interest_scales[1].determinant(22,20),
				  0.0, interest_scales[1].polarity(22,20),
				  0,1) );
      ip.push_back( InterestPoint(30,30,params.calcScale(6),
				  interest_scales[1].determinant(30,30),
				  0.0, interest_scales[1].polarity(30,30),
				  0,1) );
      ip.push_back( InterestPoint(24,24,params.calcScale(6),
				  interest_scales[1].determinant(24,24),
				  0.0, interest_scales[1].polarity(24,24),
				  0,1) );
      ip.push_back( InterestPoint(12,12,params.calcScale(6),
				  interest_scales[1].determinant(12,12),
				  0.0, interest_scales[1].polarity(12,12),
				  0,1) );
      
      // Performing SubpixelRefinement
      SURFSubpixelRefinement( ip, interest_scales, params );
      
      // Checking data
      for (InterestPointList::iterator p = ip.begin();
	   p != ip.end(); ++p ) {
	TS_ASSERT_DELTA( (*p).x, center.x(), .2 );
	TS_ASSERT_DELTA( (*p).y, center.y(), .2 );
	TS_ASSERT_DELTA( (*p).scale, params.calcScale( 6 ), .2 );
      }
      
    }
  }

  void test_orientation() {
    ImageView<float> gradient;
    read_image( gradient, "noisy_gradient_60.png");
    ImageView<double> integral = IntegralImage( gradient );

    InterestPointList ip;
    ip.push_back( InterestPoint(50,50,2.4,0.0) );
    ip.push_back( InterestPoint(50,50,4.8,0.0) );
    ip.push_back( InterestPoint(25,25,2.4,0.0) );
    ip.push_back( InterestPoint(25,25,4.8,0.0) );
    ip.push_back( InterestPoint(25,50,2.4,0.0) );
    ip.push_back( InterestPoint(25,50,4.8,0.0) );
    ip.push_back( InterestPoint(50,25,2.4,0.0) );
    ip.push_back( InterestPoint(50,25,4.8,0.0) );
    ip.push_back( InterestPoint(75,75,4.8,0.0) );
    ip.push_back( InterestPoint(75,25,9.6,0.0) );
    ip.push_back( InterestPoint(10,90,9.6,0.0) );
    ip.push_back( InterestPoint(20,80,9.6,0.0) );

    float sixty_degrees = 60*3.14159/180;

    for ( InterestPointList::iterator point = ip.begin();
	  point != ip.end(); ++point ) {
      (*point).orientation = SURFOrientation( integral,
					      (*point).x, (*point).y,
					      (*point).scale );
      TS_ASSERT_DELTA((*point).orientation, sixty_degrees, .1 );
    }

  }
};
