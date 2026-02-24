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
  // grassfire()
  // *******************************************************************

  // Computes the 4-connected grassfire image of an image.
  // (Specifically, computes the Manhattan distance from each pixel to
  // the nearest pixel with zero value, assuming the borders of the
  // image are zero.)
  template <class SourceT, class OutputT>
  void grassfire( ImageViewBase<SourceT> const& src, ImageView<OutputT>& dst,
                  bool ignore_borders ) {
    int32 cols = src.impl().cols(), rows = src.impl().rows();
    dst.set_size( cols, rows );

    // Bug fix: This code crashes if images are too small
    if (cols <= 1 || rows <= 1) {
      fill(dst, 0);
      return;
    }

    typedef typename SourceT::pixel_accessor src_accessor;
    typedef typename ImageView<OutputT>::pixel_accessor dst_accessor;

    src_accessor srow = src.impl().origin();
    dst_accessor drow = dst.origin();
    const typename SourceT::pixel_type zero = typename SourceT::pixel_type();

    if (!ignore_borders) {

      { // First row
        src_accessor scol = srow;
        dst_accessor dcol = drow;
        for( int32 col=cols; col; --col ) {
          // Pixels on the edge can only be distance one or zero.
          *dcol = ((*scol)==zero)?0:1;
          scol.next_col();
          dcol.next_col();
        }
        srow.next_row();
        drow.next_row();
      }

      // Raster through the image down to the second last row.
      for( int32 row=rows-2; row; --row ) {
        src_accessor scol = srow;
        dst_accessor dcol = drow;
        // The first pixel in the row is one or zero.
        *dcol = ((*scol)==zero)?0:1;
        scol.next_col();
        dcol.next_col();
        // Loop through columns from col to the second last col
        for( int32 col=cols-2; col; --col ) {
          if( (*scol)==zero ) // Input pixel is zero, so is the output pixel.
            (*dcol)=0;
          else { // Output count is 1 + min of the left and upper values.
            dst_accessor s1 = dcol, s2 = dcol;
            (*dcol) = 1 + std::min( *(s1.prev_col()), *(s2.prev_row()) );
          }
          scol.next_col();
          dcol.next_col();
        }
        // The last pixel in the row is again one or zero.
        *dcol = ((*scol)==zero)?0:1;
        srow.next_row();
        drow.next_row();
      }
      { // Last row
        src_accessor scol = srow;
        dst_accessor dcol = drow;
        for( int32 col=cols; col; --col ) {
          // These locations on the edge can only be one or zero.
          *dcol = ((*scol)==zero)?0:1;
          scol.next_col();
          dcol.next_col();
        }
      }

      // Raster through the image in reverse order.
      drow.advance(cols-2,-1);
      for( int32 row=rows-2; row; --row ) {
        dst_accessor dcol = drow;
        for( int32 col=cols-2; col; --col ) {
          if( (*dcol)!=0 ) {
            // If min(lower, right)+1 is less than the current value, update.
            dst_accessor s1 = dcol, s2 = dcol;
            int32 m = std::min( *(s1.next_col()), *(s2.next_row()) );
            if( m < *dcol )
              *dcol = m + 1;
          }
          dcol.prev_col();
        } // End col loop
        drow.prev_row();
      } // End row loop

    } else { // ignore_borders

      // The largest possible distance in the image with no zeros.
      OutputT init_val = cols + rows;

      // Raster through the image
      for( int32 row=0; row<rows; ++row ) {
        src_accessor scol = srow;
        dst_accessor dcol = drow;

        // Loop through columns
        for( int32 col=0; col<cols; ++col ) {
          if( (*scol)==zero ) // Input pixel is zero, so is the output pixel.
            (*dcol) = 0;
          else { // Output count is 1 + min of the left and upper values.

            dst_accessor s1 = dcol, s2 = dcol;
            if (row > 0) {
              if (col > 0)
                (*dcol) = std::min( *(s1.prev_col()), *(s2.prev_row()) ) + 1;
              else
                (*dcol) = std::min(init_val, *(s2.prev_row())+1);
            } else {
              if (col > 0)
                (*dcol) = std::min(init_val, *(s1.prev_col())+1);
              else
                (*dcol) = init_val;
            }
          }

          scol.next_col();
          dcol.next_col();
        } // End loop through columns
        srow.next_row();
        drow.next_row();
      }// End loop through rows

      // Raster through the image in reverse order.
      drow.advance(cols-1,-1);
      for( int32 row=rows; row; --row ) {
        dst_accessor dcol = drow;
        for( int32 col=cols; col; --col ) {

          if( (*dcol)!=0 ) {
            // If min(lower, right)+1 is less than the current value, update.
            dst_accessor s1 = dcol, s2 = dcol;
            OutputT m = *dcol;
            if (col < cols)
              m = std::min( m, *(s1.next_col())+1 );
            if (row < rows)
              m = std::min( m, *(s2.next_row())+1 );

            *dcol = m;
            // When ignoring borders we need to enforce a cap on the output values.
            if ( *dcol > init_val )
              *dcol = init_val;
          }

          dcol.prev_col();
        } // End col loop
        drow.prev_row();
      } // End row loop
    } // End border ignore case

  }


  // *******************************************************************
  // centerline_weights()
  // *******************************************************************


// A function that compute weights (positive in the image and zero
// outside) based on finding where each image line data values start
// and end and the centerline. The same thing is repeated for columns.
// Then two such functions are multiplied. This works better than
// grassfire for images with simple boundary and without holes.
// If fill_holes is true, this will return positive weights inside holes.
template<class ImageT>
void centerline_weights(ImageT const& img, ImageView<double> & weights,
                         BBox2i roi, bool fill_holes, bool use_min_weight){
  int numRows = img.rows();
  int numCols = img.cols();

  // Arrays to be returned out of this function
  std::vector<double> hCenterLine  (numRows, 0);
  std::vector<double> hMaxDistArray(numRows, 0);
  std::vector<double> vCenterLine  (numCols, 0);
  std::vector<double> vMaxDistArray(numCols, 0);

  std::vector<int> minValInRow(numRows, 0);
  std::vector<int> maxValInRow(numRows, 0);
  std::vector<int> minValInCol(numCols, 0);
  std::vector<int> maxValInCol(numCols, 0);

  for (int k = 0; k < numRows; k++){
    minValInRow[k] = numCols;
    maxValInRow[k] = 0;
  }
  for (int col = 0; col < numCols; col++){
    minValInCol[col] = numRows;
    maxValInCol[col] = 0;
  }

  // Note that we do just a single pass through the image to compute
  // both the horizontal and vertical min/max values.
  for (int row = 0 ; row < numRows; row++) {
    for (int col = 0; col < numCols; col++) {

      if ( !is_valid(img(col,row)) )
        continue;

      // Record the first and last valid column in each row
      if (col < minValInRow[row]) minValInRow[row] = col;
      if (col > maxValInRow[row]) maxValInRow[row] = col;

      // Record the first and last valid row in each column
      if (row < minValInCol[col]) minValInCol[col] = row;
      if (row > maxValInCol[col]) maxValInCol[col] = row;   
    }
  }

  // For each row, record central column and the column width
  for (int row = 0; row < numRows; row++) {
    hCenterLine   [row] = (minValInRow[row] + maxValInRow[row])/2.0;
    hMaxDistArray [row] =  maxValInRow[row] - minValInRow[row];
    if (hMaxDistArray[row] < 0){
      hMaxDistArray[row]=0;
    }
  }

  // For each row, record central column and the column width
  for (int col = 0 ; col < numCols; col++) {
    vCenterLine   [col] = (minValInCol[col] + maxValInCol[col])/2.0;
    vMaxDistArray [col] =  maxValInCol[col] - minValInCol[col];
    if (vMaxDistArray[col] < 0){
      vMaxDistArray[col]=0;
    }
  }

  BBox2i output_bbox = roi;
  if (roi.empty())
    output_bbox = bounding_box(img);

  // Compute the weighting for each pixel in the image
  weights.set_size(output_bbox.width(), output_bbox.height());
  fill(weights, 0);

  for (int row = output_bbox.min().y(); row < output_bbox.max().y(); row++){
    for (int col = output_bbox.min().x(); col < output_bbox.max().x(); col++){
      Vector2 pix(col, row);
      double new_weight = 0; // Invalid pixels usually get zero weight
      if (is_valid(img(col,row)) || fill_holes) {
        double weight_h = compute_line_weights(pix, true,  hCenterLine, hMaxDistArray);
        double weight_v = compute_line_weights(pix, false, vCenterLine, vMaxDistArray);
        if (use_min_weight)
          new_weight = std::min(weight_h, weight_v);
        else
          new_weight = weight_h*weight_v;
      }
      weights(col-output_bbox.min().x(), row-output_bbox.min().y()) = new_weight;

    }
  }

} // End function centerline_weights



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
