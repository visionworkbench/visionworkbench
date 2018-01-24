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


/// \file Algorithms2.h
///
/// Some very simple code which needs to be tightened. It is used only
/// in very few places and is not very generic, hence not in the main
/// Algorithms.h.
///
#ifndef __VW_IMAGE_ALGORITHMS2_H__
#define __VW_IMAGE_ALGORITHMS2_H__

#include <vw/Image/Algorithms.h>

namespace vw {

// Round a double value to int, and clamp to the min and max of this type of int.
template<class IntType>
IntType round_and_clamp(double val){

  // Round. Be careful to work with doubles before clamping, to not overflow.
  val = round(val);

  // Clamp. Note that this won't work right unless IntType is an int type.
  if (val < double(std::numeric_limits<IntType>::min()) )
    val = double(std::numeric_limits<IntType>::min());
  if (val > double(std::numeric_limits<IntType>::max()) )
    val = double(std::numeric_limits<IntType>::max());

  return static_cast<IntType>(val);
}

// Pick the first channel of an image. Round its values, then clamp to the bounds
// for the given output type and cast to this output type.
template <class ImageT, class OutputPixelT>
class RoundAndClamp: public ImageViewBase< RoundAndClamp<ImageT, OutputPixelT> >{
  ImageT const& m_img;

public:
  RoundAndClamp( ImageT const& img): m_img(img){}

  typedef typename ImageT::pixel_type InputPixelT;
  typedef OutputPixelT pixel_type;
  typedef OutputPixelT result_type;
  typedef ProceduralPixelAccessor<RoundAndClamp> pixel_accessor;

  inline int32 cols() const { return m_img.cols(); }
  inline int32 rows() const { return m_img.rows(); }
  inline int32 planes() const { return 1; }

  inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

  inline result_type operator()( double/*i*/, double/*j*/, int32/*p*/ = 0 ) const {
    vw_throw(NoImplErr() << "RoundAndClamp::operator()(...) is not implemented");
    return result_type();
  }

  typedef CropView<ImageView<result_type> > prerasterize_type;
  inline prerasterize_type prerasterize(BBox2i const& bbox) const {

    ImageView<result_type> tile(bbox.width(), bbox.height());
    for (int col = bbox.min().x(); col < bbox.max().x(); col++){
      for (int row = bbox.min().y(); row < bbox.max().y(); row++){

        // I could not figure out in reasonable time how to select a channel
        // from a pixel which can be single or compound. Hence use the
        // functionality for selecting a channel from an image.
        ImageView<InputPixelT> A(1, 1);
        A(0, 0) = m_img(col, row);

        // First cast to double, as the values of A could be out of range
        // for the OutputPixelT data type.
        ImageView<double> B = select_channel(A, 0);

        // Now round, clamp, and cast. 
        tile(col - bbox.min().x(), row - bbox.min().y() ) = round_and_clamp<OutputPixelT>(B(0, 0));
      }
    }
    
    return prerasterize_type(tile, -bbox.min().x(), -bbox.min().y(),
                             cols(), rows() );
  }

  template <class DestT>
  inline void rasterize(DestT const& dest, BBox2i bbox) const {
    vw::rasterize(prerasterize(bbox), dest, bbox);
  }
};
template <class ImageT, class OutputPixelT>
RoundAndClamp<ImageT, OutputPixelT> round_and_clamp(ImageT const& img){
  return RoundAndClamp<ImageT, OutputPixelT>(img);
}

} // namespace vw

#endif // __VW_IMAGE_ALGORITHMS2_H__
