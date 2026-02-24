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
  // fill()
  // *******************************************************************

  /// Fill an image with a constant pixel value
  template <class ImageT, class ValueT>
  void fill( ImageViewBase<ImageT> const& image, ValueT value_ ) {
    int32 planes=image.impl().planes(), rows=image.impl().rows(), cols=image.impl().cols();
    typename ImageT::pixel_type value( value_ );
    typename ImageT::pixel_accessor plane = image.impl().origin();
    for( int32 p=planes; p; --p ) {
      typename ImageT::pixel_accessor row = plane;
      for( int32 r=rows; r; --r ) {
        typename ImageT::pixel_accessor col = row;
        for( int32 c=cols; c; --c ) {
          *col = value;
          col.next_col();
        }
        row.next_row();
      }
      plane.next_plane();
    }
  }

  // *******************************************************************
  // nonzero_data_bounding_box()
  // *******************************************************************

  template <class ViewT>
  BBox2i nonzero_data_bounding_box( ImageViewBase<ViewT> const& image_ ) {
    const typename ViewT::pixel_type zero = typename ViewT::pixel_type();
    ViewT const& image = static_cast<ViewT const&>(image_);
    int32 x=0, y=0, cols=0, rows=0;
    int32 i, j, icols = image.cols(), irows = image.rows();
    for( j=0; j<irows; ++j ) {
      for( i=0; i<icols; ++i ) {
        if( image(i,j) != zero ) break;
      }
      if( i != icols ) break;
    }
    if( j != irows ) {
      y = j;
      for( j=irows-1; j; --j ) {
        for( i=0; i<icols; ++i ) {
          if( image(i,j) != zero ) break;
        }
        if( i != icols ) break;
      }
      rows = j - y + 1;
      for( i=0; i<icols; ++i ) {
        for( j=y; j<y+rows; ++j ) {
          if( image(i,j) != zero ) break;
        }
        if( j != y+rows ) break;
      }
      x = i;
      for( i=icols-1; i; --i ) {
        for( j=y; j<y+rows; ++j ) {
          if( image(i,j) != zero ) break;
        }
        if( j != y+rows ) break;
      }
      cols = i - x + 1;
    }
    return BBox2i( x, y, cols, rows );
  }
}
