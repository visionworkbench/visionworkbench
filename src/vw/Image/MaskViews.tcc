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

namespace vw {

  // *******************************************************************
  /// create_mask( view, value )
  ///
  /// Given a view with pixels of type PixelT and a pixel value to
  /// consider as the "no data" or masked value, returns a view with
  /// pixels that are of the PixelMask<PixelT>, with the appropriate
  /// pixels masked.
  /// - Should safely accept inputs which already contain a mask, in 
  ///   which case the input mask is ignored.
  
  // *******************************************************************
  /// apply_mask( view, value )
  ///
  /// Given a view with pixels of the type PixelMask<T>, this view
  /// returns an image with pixels of type T where any pixel that was
  /// marked as "invalid" in the mask is replaced with the constant
  /// pixel value passed in as value.  The value is T() by default.
  ///
  // *******************************************************************
  /// copy_mask(view, mask)
  ///
  /// Copies a mask from one image to another.
  ///

  // *******************************************************************
  /// mask_to_alpha(view)
  ///
  /// Converts a mask channel to an alpha channel, generating an image that
  /// is transparent wherever the data is masked.
  ///

  // *******************************************************************
  /// alpha_to_mask(view)
  ///
  /// Converts an channel to a mask channel, generating an image that
  /// is transparent wherever the data has an alpha value of 0.
  ///

  // *******************************************************************
  /// EdgeMaskView
  ///
  /// Create an image with zero-valued (i.e. default constructor)
  /// pixels around the edges masked out.

  //******************************************************************
  /// invert_mask(view)
  ///
  /// Given a view with pixels of type PixelMask<T>, this will toggle
  /// all valids to invalid, and invalid to valids.
  template <class PixelT>
  class InvertPixelMask : public ReturnFixedType<PixelT> {
  public:
    inline PixelT operator()( PixelT value ) const {
      toggle(value);
      return value;
    }
  };

  //*****************************************************************
  /// validate_mask(view)
  ///
  /// Given an image of PixelMasks, this will make all pixels valid
  template <class PixelT>
  class ValidatePixelMask : public ReturnFixedType<PixelT> {
  public:
    inline PixelT operator()( PixelT value ) const {
      validate(value);
      return value;
    }
  };

  //*****************************************************************
  /// invalidate_mask(view)
  ///
  /// Given an image of PixelMasks, this will make all pixels invalid
  template <class PixelT>
  class InvalidatePixelMask : public ReturnFixedType<PixelT> {
  public:
    inline PixelT operator()( PixelT value ) const {
      invalidate(value);
      return value;
    }
  };

  //******************************************************************
  /// union_mask(view, mask)
  ///
  /// Unions 'mask' w/ view. View's data is returned
  ///
  template <class PixelT>
  class UnionPixelMask : public ReturnFixedType<typename MaskedPixelType<PixelT>::type> {
    typedef typename MaskedPixelType<PixelT>::type return_type;
  public:
    template <class MaskedPixelT>
    inline return_type operator()( PixelT const& value, MaskedPixelT const& mask ) const {
      return_type result = value;
      if ( is_valid(value) || is_valid(mask) )
        validate(result);
      else
        invalidate(result);
      return result;
    }
  };

  //******************************************************************
  /// intersect_mask(view, mask)
  ///
  /// Intersects 'mask' w/ view. View's data is returned
  ///
  template <class PixelT>
  class IntersectPixelMask : public ReturnFixedType<typename MaskedPixelType<PixelT>::type> {
    typedef typename MaskedPixelType<PixelT>::type return_type;
  public:
    template <class MaskedPixelT>
    inline return_type operator()( PixelT const& value, MaskedPixelT const& mask ) const {
      return_type result = value;
      if ( is_valid(value) && is_valid(mask) )
        validate(result);
      else
        invalidate(result);
      return result;
    }
  };

} // namespace vw

