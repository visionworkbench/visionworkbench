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
#ifndef __VW_INTERESTPOINT_SURF_H__
#define __VW_INTERESTPOINT_SURF_H__

#include <vw/InterestPoint/IntegralImage.h>
#include <vw/InterestPoint/InterestData.h>
#include <vw/InterestPoint/Detector.h>
#include <vw/InterestPoint/Descriptor.h>
#include <vw/Math.h>
#include <vw/Image.h>
#include <vw/Core/ThreadPool.h>

#include <vector>
#include <list>

namespace vw {
namespace ip {

  /// Red Herring
  class SURFInterestOperator {
    
    static const float DEFAULT_INTEREST_THRESHOLD = 0.03;
    float m_threshold;
  public:
    SURFInterestOperator() : m_threshold(DEFAULT_INTEREST_THRESHOLD) {}
    SURFInterestOperator( float threshold ) : m_threshold(threshold) {}
   
    /// DO NOT USE, this is only for compatibility with IP code
    template <class ViewT>
    inline ImageViewRef< typename ViewT::pixel_type>
    operator() (ImageViewBase<ViewT> const& source, float scale= 1.0) const {
      vw_out(DebugMessage, "interest_point") << "DO NOT USE SURFInterestOperator it does nothing.\n";
    }

    template <class DataT>
    inline bool threshold (InterestPoint const& pt, DataT const& data) const {
      return false;
    }

    template <class DataT>
    inline bool threshold (DataT const& data) const {
      return data > m_threshold;
    }

    float threshold ( void ) const {
      return m_threshold;
    }
  };

  /// SURFParams
  /// Just a fancier way to share data with the functions
  class SURFParams {
  public:
    SURFParams( float thres = .1, unsigned oct = 4, unsigned scal = 4 ) {
      cols = rows = 0;    // Remember to set this later
      threshold = thres;
      octaves = oct;
      scales = scal;
      lobe_size = 3;
      
      // Building indexes
      indexes = new unsigned[octaves * scales];
      unsigned count = 0;
      for (unsigned o = 0; o < octaves; o++ ) {
	for (unsigned s = 0; s < scales; s++ ) {
	  if ( o == 0 ) {
	    indexes[o*scales + s] = count;
	    count++;
	  } else if ( s == 0 ) {
	    indexes[o*scales] = indexes[o*scales - 3];
	  } else if ( s == 1 ) {
	    indexes[o*scales+1] = indexes[o*scales-1];
	  } else {
	    indexes[o*scales + s] = count;
	    count++;
	  }
	}
      }

    }
    template <class ViewT>
    void setForImage( ImageViewBase<ViewT> const& source ) {
      ImageView<float> image = pixel_cast<PixelGray<float> >(source);

      cols = image.cols();
      rows = image.rows();
    }

    template <class T>
    float calcScale( T const& filter_size ) const {
      return 1.2*float(filter_size)/9.0;
    }

    unsigned convToIndex( unsigned const& octave, unsigned const& scale ) const {
      unsigned idx = octave*scales + scale;
      VW_DEBUG_ASSERT( idx < 16, 
		       vw::ArgumentErr() << "convToIndex: Incorrect octave and scale given" );
      return indexes[idx];
    }
    int cols;
    int rows;
    float threshold;
    unsigned octaves;
    unsigned scales;
    unsigned lobe_size;
  private:
    // This is used for converting octave and scale to an index of the
    // interest data. It is fixed for 4 octaves and 4 scales right now
    unsigned *indexes;
  };

  /// SURFScaleData
  /// This is a wrapper for the interest data to make it easier to
  /// compare different scales to each other
  class SURFScaleData {
  public:
    SURFScaleData( void ) {
    }

    // Fancy constructor
    SURFScaleData( boost::shared_ptr<vw::ImageView<float> > data,
		   boost::shared_ptr<vw::ImageView<bool> > polarity,
		   unsigned const& filter_size, unsigned const& sampStep) :
    m_data(data), m_polarity(polarity), m_filter_size(filter_size), m_sampling_step(sampStep)
    { 
      // calculating sampling_step_2, (finds x of 2^x = sampling_step)
      unsigned temp = m_sampling_step;
      m_sampling_step_2 = 0;
      while ( (temp & 0x1) == 0 && temp != 0 ) {
	m_sampling_step_2++;
	temp>>=1;
      }

      // Calculating the starting offset
      unsigned s_temp = floor(float(m_filter_size)/2.);
      m_starting_offset = s_temp + (m_sampling_step - (s_temp & ( m_sampling_step - 1)));
    }

    float determinant( signed const& x, signed const& y) const {
      vw::Vector3i cov = convertCoord( x, y);
      if ( cov[2] == 0 ){
	return 0.0;
      }
      return m_data->operator()(cov[0],cov[1]);
    }

    bool polarity( signed const& x, signed const& y) const {
      vw::Vector3i cov = convertCoord( x, y);
      if ( cov[2] == 0 )
	return false;
      return m_polarity->operator()(cov[0],cov[1]);
    }

    bool isValid( signed const& x, signed const& y) const {
      vw::Vector3i cov = convertCoord( x, y);
      if ( cov[2] == 0 )
	return false;
      return true;
    }
    
    // isLessThan
    // - value = value to compare interest image against
    // - x,y   = point to compare value against 
    //           (in the coordinates of the orginal image)
    bool isLessThan( float const& value,
		     int const& x, int const& y) const {
      vw::Vector3i cov = convertCoord( x, y);
      if ( cov[2] == 0 )
	return false;
      return value < m_data->operator()(cov[0],cov[1]);
    }

    unsigned filter_size( void ) const {
      return m_filter_size;
    }
  private:
    unsigned m_filter_size;
    unsigned m_starting_offset;
    unsigned m_sampling_step;
    unsigned m_sampling_step_2;
    boost::shared_ptr<vw::ImageView<float> > m_data;
    boost::shared_ptr<vw::ImageView<bool> > m_polarity;

    vw::Vector3i convertCoord( int const& x, int const& y) const {
      if ( (x & 1) || (y & 1) ) // is it odd?
	return vw::Vector3i(0,0,0);
      int tx = x, ty = y;
      tx -= m_starting_offset;
      ty -= m_starting_offset;
      if (tx < 0 || ty < 0) // Have we left our space?
	return vw::Vector3i(0,0,0);
      if (tx & (m_sampling_step - 1) || ty & (m_sampling_step - 1) ) // mod not equal zero
	return vw::Vector3i(0,0,0);
      tx = tx >> m_sampling_step_2;
      ty = ty >> m_sampling_step_2;
      if ( tx >= m_data->cols() || ty >= m_data->rows() )
	return vw::Vector3i(0,0,0);
      return vw::Vector3i(tx,ty,1);
    }
  };

  /// Creates a interest scale for SURF
  SURFScaleData SURFProcessScale( vw::ImageView<double> const& integral, 
				  unsigned const& filter_size,
				  unsigned const& sampling_step,
				  SURFParams const& params );

  /// Finds the maximas in our data
  std::list<InterestPoint> SURFMaximaDetection( std::vector<SURFScaleData> const & scaleData, 
						SURFParams const& params );
  
  /// Writes SURF debug interest images
  void SURFWriteDebugImages ( std::vector<SURFScaleData> const& scaleData,
			      SURFParams const& params );

  /// SURF Subpixel refinement
  /// This gets floating point precision for the interest points, will
  /// also throw at some points if they move too far.
  void SURFSubpixelRefinement( std::list<vw::ip::InterestPoint>& ip,
			       std::vector<SURFScaleData> const& scaleData,
			       SURFParams const& params );

  // SURF Gradient 3D
  // This calculates the gradient, hmm... it negative
  Vector3 SURFGradient3D( std::vector<SURFScaleData> const& scaleData,
			  vw::ip::InterestPoint const& ip,
			  SURFParams const& params );

  // SURF Hessian 3D
  // This calculates the hessian matrix for sub pixel correlation
  Matrix3x3 SURFHessian3D( std::vector<SURFScaleData> const& scaleData, 
			   vw::ip::InterestPoint const& ip, 
			   SURFParams const& params );

  // SURF Orientation
  // This calculates the orientation of a feature
  float SURFOrientation( vw::ImageView<double> const&, int const&, int const&,
			 float const& );  

  // SURF Descriptor
  // This calculated the descriptor of a feature
  Vector<float> SURFDescriptor( vw::ImageView<double> const&,
				Matrix<float,20,20> const&,
				vw::ip::InterestPoint const& ip,
				bool extended );

  template <class ViewT>
  class SURFInterestScaleTask : public Task {
    // Disable copyable semantics
    SURFInterestScaleTask(SURFInterestScaleTask& copy) {}
    void operator=(SURFInterestScaleTask& copy) {}

    ViewT m_integral;
    SURFScaleData& m_data;
    Mutex& m_mutex;
    unsigned m_filter_size;
    unsigned m_sampling_step;
    SURFParams m_params;
    int m_id;

  public:
    SURFInterestScaleTask(ViewT const& integral, unsigned const filter_size,
			  unsigned sampling_step, SURFParams const params,
			  SURFScaleData& scaleData, Mutex &mutex, int id) :
    m_integral(integral), m_filter_size(filter_size), m_sampling_step( sampling_step ), m_params(params), m_data(scaleData), m_mutex(mutex), m_id(id) {
    }

    // This actually does the work
    void operator()() {

      m_data = SURFProcessScale( m_integral,
				 m_filter_size,
				 m_sampling_step,
				 m_params );

    }

    SURFScaleData results(){ return m_data; }
  };

  template <class ViewT>
  class SURFOrientationCalcTask : public Task {
    // Disable copyable semantics
    SURFOrientationCalcTask(SURFOrientationCalcTask& copy) {}
    void operator=(SURFOrientationCalcTask& copy) {}

    ViewT m_integral;
    InterestPointList& m_individual_list;
    InterestPointList& m_main_list;
    Mutex& m_mutex;
    int m_id;

  public:
    SURFOrientationCalcTask(ViewT const& integral, InterestPointList& my_ip_list, 
			    InterestPointList& main_ip_list, Mutex &mutex, int id) :
    m_integral(integral), m_individual_list(my_ip_list), m_main_list(main_ip_list), m_mutex(mutex), m_id(id) {
    }

    void operator()(){

      // Processing my individual list of interest points
      for (InterestPointList::iterator point = m_individual_list.begin();
	   point != m_individual_list.end(); ++point) {
	(*point).orientation = SURFOrientation( m_integral,
						(*point).ix, (*point).iy,
						(*point).scale );
      }

      // Injecting my finished points into the global pile
      {
	Mutex::Lock lock(m_mutex);
	m_main_list.splice(m_main_list.end(), m_individual_list);
      } // End mutex lock's scope

    }
  };

  template <class ViewT>
  class SURFDescriptorTask : public Task {
    // Disable copyable semantics
    SURFDescriptorTask( SURFDescriptorTask& copy ) {}
    void operator=(SURFDescriptorTask& copy ) {}

    ViewT m_integral;
    Matrix<float,20,20> m_gauss;
    bool m_extended;
    InterestPointList& m_individual_list;
    InterestPointList& m_main_list;
    Mutex& m_mutex;
    int m_id;
    
  public:
    SURFDescriptorTask( ViewT const& integral,  
			Matrix<float,20,20>& gauss,
			bool extended,
			InterestPointList& my_ip_list,
			InterestPointList& main_ip_list, 
			Mutex &mutex, int id) :
    m_integral(integral), m_gauss(gauss), m_individual_list(my_ip_list), 
      m_extended(extended), m_main_list(main_ip_list), m_mutex(mutex), m_id(id) {
    }

    void operator()(){
      // Processing my individual list of interest points
      for (InterestPointList::iterator point = m_individual_list.begin();
	   point != m_individual_list.end(); ++point) {
	(*point).descriptor = SURFDescriptor( m_integral, m_gauss,
					      (*point), m_extended );
					      
      }

      // Injecting my finished points into the global pile
      {
	Mutex::Lock lock(m_mutex);
	m_main_list.splice(m_main_list.end(), m_individual_list);
      }
    }
  };

  /// This class actually performs all the work. The SURF operator is
  /// actually just a red herring that is required to match the IP module's
  /// frame.
  template <class InterestT>
  class FH9InterestPointDetector : public InterestDetectorBase<FH9InterestPointDetector<InterestT> >
  {

    FH9InterestPointDetector(FH9InterestPointDetector<InterestT> const& copy) {}
    ImageView<double> m_integral;

  public:
    static const int IP_DEFAULT_SCALES = 4; // Don't play with
    static const int IP_DEFAULT_OCTAVES = 4;

    /// Setting max_points = 0 will disable interest point culling.
    /// Otherwies, the max_points most "interesting" points are
    /// returned.
    FH9InterestPointDetector(ImageView<double> const& integral, int max_points = 0, int num_threads = 1)
      : m_max_points(max_points), m_num_threads(num_threads), m_integral(integral) {
      m_params = SURFParams();
    }

    FH9InterestPointDetector(InterestT const& interest, ImageView<double> const& integral, int max_points = 0, int num_threads = 1 ) 
      : m_max_points(max_points), m_num_threads(num_threads), m_integral(integral) {
      
      m_params = SURFParams( interest.threshold(), IP_DEFAULT_OCTAVES, IP_DEFAULT_SCALES );
    }

    FH9InterestPointDetector(InterestT const& interest, ImageView<double> const& integral, int scales, int octaves, int max_points, int num_threads = 1 )
      : m_max_points(max_points), m_num_threads(num_threads), m_integral(integral) {
      
      m_params = SURFParams( interest.threshold(), octaves, scales );
    }

    /// Process Image ////////////////////////////
    /// Detect interest points in the source image.
    template <class ViewT>
    InterestPointList process_image(ImageViewBase<ViewT> const& image) {

      // Timing
      Timer *total = new Timer("Total elapsed time", DebugMessage, "interest_point");

      // Set SURF Params;
      m_params.setForImage( image );
  
      // Create Interest Scales
      vw_out(InfoMessage,"interest_point") << "\tCalculating Interest Data.\n";
      std::vector<SURFScaleData> interest_scales;
      {
	unsigned filter_size = 3;
	int id_count = 0;
	FifoWorkQueue queue( m_num_threads );
	Mutex mutex;

	// Determing how many interest scales are required
	{
	  int size = 4;
	  size += m_params.octaves*2;
	  interest_scales.resize(size);
	}

	for ( unsigned octave = 0; octave < m_params.octaves; octave++ ) {
	  unsigned sampling_step = 0x2 << octave;
	  unsigned iterations = 2;
	  if ( octave == 0 )
	    iterations = 4;
       
	  for ( unsigned i = 0; i < iterations; i++ ) {
	    filter_size += 6 * ( octave == 0 ? 1 : 0x2 << ( octave -1 ) );
	    
	    if ( filter_size + 2*sampling_step > m_integral.cols() ||
		 filter_size + 2*sampling_step > m_integral.rows() ) {
	      m_params.octaves--;
	      break;
	    }
	    
	    if ( m_num_threads < 2 ) {
	      // Single Threaded Version
	      interest_scales[id_count] = SURFProcessScale( m_integral,
							    filter_size,
							    sampling_step,
							    m_params );
	    } else {
	      // Multithread Version
	      typedef SURFInterestScaleTask<ImageView<double> > task_type;
	      boost::shared_ptr<task_type> task ( new task_type( m_integral, filter_size,
								 sampling_step, m_params, 
								 interest_scales[id_count],
								 mutex, id_count ) );
	      queue.add_task( task );
	    }

	    id_count++;
	  }
	}

	if ( m_num_threads > 1 ) {

	  queue.join_all();
	}
      }
      
      // Write Debug Images ?
      // SURFWriteDebugImages( interest_scales, m_params );
      
      // Locate Maximas
      vw_out(InfoMessage,"interest_point") << "\tPerforming non-maximal suppression.\n";
      std::list<InterestPoint> ip = SURFMaximaDetection( interest_scales,
							 m_params );

      vw_out(InfoMessage,"interest_point") << "\t\toutput points: " << ip.size() << std::endl;
      
      // Subpixel refinement
      vw_out(InfoMessage,"interest_point") << "\tPerforming sub-pixel refinement.\n";
      SURFSubpixelRefinement( ip, interest_scales, m_params );
      
      vw_out(InfoMessage,"interest_point") << "\t\toutput points: " << ip.size() << std::endl;
      
      // Cull
      if ( m_max_points != 0 ) {
	vw_out(InfoMessage,"interest_point") << "\tCulling data to a " << m_max_points << " interest points.\n";
	int orginal_num_points = ip.size();
	ip.sort();
	if ( m_max_points < ip.size() )
	  ip.resize( m_max_points );
	vw_out(DebugMessage,"interest_point") << "\tdone (removed " << orginal_num_points - ip.size() << " interest points, " << ip.size() << " remaining.).\n";
      }
      
      // Assign orientations
      vw_out(InfoMessage,"interest_point") << "\tCalculating orientation.\n";
      {
	if ( m_num_threads < 2 ) {
	  // Single threaded edition
	  for (std::list<InterestPoint>::iterator point = ip.begin();
	       point != ip.end(); ++point) {
	    (*point).orientation = SURFOrientation( m_integral,
						    (*point).ix, (*point).iy,
						    (*point).scale );
	  }
	} else {
	  // Multithreaded

	  // Breaking the points in lengths that can be processed individually.
	  unsigned count = 0;
	  std::vector<InterestPointList> ip_separated(m_num_threads);
	  for (InterestPointList::iterator point = ip.begin();
	       point != ip.end(); ++point) {
	    
	    ip_separated[count].push_back( (*point) );

	    count++;
	    count %= m_num_threads;
	  }
	  ip.clear();
	  
	  // Giving data to threads
	  FifoWorkQueue queue( m_num_threads );
	  Mutex mutex;
	  typedef SURFOrientationCalcTask<ImageView<double> > task_type;
	  for ( count = 0; count < m_num_threads; count++ ) {
	    boost::shared_ptr<task_type> task ( new task_type( m_integral, ip_separated[count],
							       ip, mutex, count ) );

	    queue.add_task( task );
	  }

	  queue.join_all();
	}
      }
      
      // End Timing
      delete total;

      return ip;
    }


  protected:
    int m_num_threads;
    int m_max_points;
    SURFParams m_params;
  };

  // SURF Descriptor 64/128 bit
  struct SURFDescriptorGenerator {

    bool m_extended;
    int m_num_threads;
    Matrix<float,20,20> gaussian_weight;

  public:

    // Constructor
    SURFDescriptorGenerator( bool extended=false, int num_threads = 1 ) : 
    m_extended(extended), m_num_threads(num_threads) {
      for (char x = 0; x < 20; x++ ) {
	for (char y = 0; y < 20; y++ ) {
	  float dist = (float(x) - 9.5)*(float(x) - 9.5) + (float(y) - 9.5)*(float(y) - 9.5);
	  
	  gaussian_weight(x,y) = exp( -dist/21.78 );
	}
      }
    }

    // Given an image and a list of interest points, set the
    // descriptor field of the interest points. Sorry that this
    // doesn't fit into the Descriptor Generator Base.
    void operator()( ImageView<double> const& integral, InterestPointList& points ) {
      
      // Timing
      Timer *total = new Timer("Total elapsed time", DebugMessage, "interest_point");

      if ( m_num_threads < 2 ) {
	// Single Threaded
	for (InterestPointList::iterator i = points.begin(); 
	     i != points.end(); ++i ) {
	  (*i).descriptor = SURFDescriptor( integral,
					    gaussian_weight,
					    (*i), m_extended );
	}
      } else {
	// Multithreaded

	// Breaking the points in lengths that can be processed individually.
	unsigned count = 0;
	std::vector<InterestPointList> ip_separated(m_num_threads);
	for (InterestPointList::iterator i = points.begin();
	     i != points.end(); ++i ) {
	  
	  ip_separated[count].push_back( (*i) );

	  count++;
	  count %= m_num_threads;
	}
	points.clear();

	// Giving data to threads
	FifoWorkQueue queue( m_num_threads );
	Mutex mutex;
	typedef SURFDescriptorTask<ImageView<double> > task_type;
	for ( count = 0; count < m_num_threads; count++ ) {
	  boost::shared_ptr<task_type> task ( new task_type( integral,
							     gaussian_weight,
							     m_extended,
							     ip_separated[count],
							     points,
							     mutex, count ) );
	  queue.add_task( task );
	}
	
	queue.join_all();
      }

      // End Timing
      delete total;
    }

  };

}} // namespace vw::ip

#endif // __VW_INTERESTPOINT_SURF_H__
