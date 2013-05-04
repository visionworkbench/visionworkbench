// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


/// \file BoxFilter.h
///
/// Classes that evaluate a filter designed by summing boxes that
/// using an integral image. Immediate use is in OBALoG. Users could
/// also apply this to CenSuRE (<sp?) or SURF.
///
#ifndef __VW_INTERESTPOINT_BOX_FILTER_H__
#define __VW_INTERESTPOINT_BOX_FILTER_H__

#include <vector>

#include <vw/Image/ImageView.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/EdgeExtension.h>

namespace vw {
namespace ip {

  // Defining Box Filter Object
  // _____________________________________________________________
  struct SumBox {
    Vector2i start; // Indexed from center of filter
    Vector2i size;
    float weight;
  };
  typedef std::vector<SumBox> BoxFilter;

  // Functor for evaluating a Box Filter at origin
  // _____________________________________________________________
  template <class IntegralAccessT>
  typename IntegralAccessT::pixel_type
  inline apply_box_filter_at_point( IntegralAccessT const& src,
                                    BoxFilter const& box ) {
    typename IntegralAccessT::pixel_type result = 0;
    for ( size_t b = 0; b < box.size(); b++ ) {
      IntegralAccessT instance = src;
      typename IntegralAccessT::pixel_type box_sum = 0;
      instance.advance( box[b].start[0],
                        box[b].start[1] );
      box_sum += *instance;
      instance.advance( box[b].size[0], 0 );
      box_sum -= *instance;
      instance.advance( -box[b].size[0], box[b].size[1] );
      box_sum -= *instance;
      instance.advance( box[b].size[0], 0 );
      box_sum += *instance;
      result += box[b].weight * box_sum;
    }
    return result;
  }

  // BoxFilterView
  // _____________________________________________________________
  template <class IntegralT>
  class BoxFilterView : public ImageViewBase<BoxFilterView<IntegralT> > {
  private:
    IntegralT m_integral;
    BoxFilter const& m_filter;
    int m_max_filter_size;
    int m_pixel_buffer;

  public:
    /// The pixel type of the image view.
    typedef typename IntegralT::pixel_type pixel_type;
    /// Return type.
    typedef pixel_type result_type;
    /// Accessor Type.
    typedef ProceduralPixelAccessor<BoxFilterView<IntegralT> > pixel_accessor;

    /// Construct BoxFilterView
    BoxFilterView( IntegralT const& image, BoxFilter const& box ) :
    m_integral( image ), m_filter( box ) {
      m_max_filter_size = 0;
      for ( size_t b = 0; b < m_filter.size(); b++ ) {
        if ( m_max_filter_size < m_filter[b].size[0] )
          m_max_filter_size = m_filter[b].size[0];
        if ( m_max_filter_size < m_filter[b].size[1] )
          m_max_filter_size = m_filter[b].size[1];
      }
      m_pixel_buffer = m_max_filter_size >> 1;
    }

    /// Standard access
    inline int32 cols() const { return m_integral.cols()-1; }
    inline int32 rows() const { return m_integral.rows()-1; }
    inline int32 planes() const { return m_integral.planes(); }
    inline pixel_accessor origin() const { return pixel_accessor( *this ); }

    /// Evaluate
    inline result_type operator()( int32 x, int32 y, int32 p=0 ) const {
      if ( x < m_pixel_buffer || x >= m_integral.cols()-m_pixel_buffer-1 ||
           y < m_pixel_buffer || y >= m_integral.rows()-m_pixel_buffer-1 )
        return result_type();
      return apply_box_filter_at_point( m_integral.origin().advance(x,y,p),
                                        m_filter );
    }

    /// Rasterize function
    typedef BoxFilterView<typename IntegralT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
      return prerasterize_type( m_integral.prerasterize(bbox), m_filter );
    }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
      vw::rasterize( prerasterize(bbox), dest, bbox);
    }
  };

  // Convenience Wrappers
  // _____________________________________________________
  template <class ImageT>
  BoxFilterView<ImageT> box_filter( ImageViewBase<ImageT> const& integral,
                                    BoxFilter const& box ) {
    return BoxFilterView<ImageT>( integral.impl(), box );
  }

}}

#endif//__VW_INTERESTPOINT_BOX_FILTER_H__

