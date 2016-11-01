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


/// \file AlgorithmFunctions.h
///
/// Basic algorithms operating on images. This is for functions that are non lazy.
///
#ifndef __VW_IMAGE_ALGORITHM_FUNCTIONS_H__
#define __VW_IMAGE_ALGORITHM_FUNCTIONS_H__

#include <vw/Image/ImageView.h>

/// Used in blobindex
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

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
  void grassfire( ImageViewBase<SourceT> const& src, ImageView<OutputT>& dst ) {
    int32 cols = src.impl().cols(), rows = src.impl().rows();
    dst.set_size( cols, rows );

    // Bug fix: This code crashes if images are too small
    if (cols <= 1 || rows <= 1) {
      for (int col = 0; col < cols; col++) {
	for (int row = 0; row < rows; row++) {
	  dst(col, row) = 0;
	}
      }
      return;
    }
    
    typedef typename SourceT::pixel_accessor src_accessor;
    typedef typename ImageView<OutputT>::pixel_accessor dst_accessor;

    src_accessor srow = src.impl().origin();
    dst_accessor drow = dst.origin();
    const typename SourceT::pixel_type zero = typename SourceT::pixel_type();

    { // First row
      src_accessor scol = srow;
      dst_accessor dcol = drow;
      for( int32 col=cols; col; --col ) {
        *dcol = ((*scol)==zero)?0:1;
        scol.next_col();
        dcol.next_col();
      }
      srow.next_row();
      drow.next_row();
    }
    for( int32 row=rows-2; row; --row ) {
      src_accessor scol = srow;
      dst_accessor dcol = drow;
      *dcol = ((*scol)==zero)?0:1;
      scol.next_col();
      dcol.next_col();
      for( int32 col=cols-2; col; --col ) {
        if( (*scol)==zero ) (*dcol)=0;
        else {
          dst_accessor s1 = dcol, s2 = dcol;
          (*dcol) = 1 + std::min( *(s1.prev_col()), *(s2.prev_row()) );
        }
        scol.next_col();
        dcol.next_col();
      }
      *dcol = ((*scol)==zero)?0:1;
      srow.next_row();
      drow.next_row();
    }
    { // Last row
      src_accessor scol = srow;
      dst_accessor dcol = drow;
      for( int32 col=cols; col; --col ) {
        *dcol = ((*scol)==zero)?0:1;
        scol.next_col();
        dcol.next_col();
      }
    }
    drow.advance(cols-2,-1);
    for( int32 row=rows-2; row; --row ) {
      dst_accessor dcol = drow;
      for( int32 col=cols-2; col; --col ) {
        if( (*dcol)!=0 ) {
          dst_accessor s1 = dcol, s2 = dcol;
          int32 m = std::min( *(s1.next_col()), *(s2.next_row()) );
          if( m < *dcol ) *dcol = m + 1;
        }
        dcol.prev_col();
      }
      drow.prev_row();
    }
  }

  // Without destination given, return in a newly-created ImageView<int32>
  template <class SourceT>
  ImageView<int32> grassfire( ImageViewBase<SourceT> const& src ) {
    int32 cols = src.impl().cols(), rows = src.impl().rows();
    ImageView<int32> result( cols, rows );
    grassfire( src, result );
    return result;
  }

  // *******************************************************************
  // bounding_box()
  // *******************************************************************

  template <class ViewT>
  BBox2i bounding_box( ImageViewBase<ViewT> const& image_ ) {
    return BBox2i( 0, 0, image_.impl().cols(), image_.impl().rows() );
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


  // *******************************************************************
  // is_opaque()
  // *******************************************************************

  template <class ImageT>
  bool is_opaque_helper( ImageT const& image, true_type ) {
    for( int32 y=0; y<image.rows(); ++y )
      for( int32 x=0; x<image.cols(); ++x )
        if( ! (is_opaque( image(x,y) ) ) )
          return false;
    return true;
  }

  template <class ImageT>
  bool is_opaque_helper( ImageT const& /*image*/, false_type ) {
    return true;
  }

  /// Returns true if the given image is entirely opaque, or false if
  /// it is at least partially transparent.
  template <class ImageT>
  bool is_opaque( ImageViewBase<ImageT> const& image ) {
    return is_opaque_helper( image.impl(), typename PixelHasAlpha<typename ImageT::pixel_type>::type() );
  }


  // *******************************************************************
  // is_transparent()
  // *******************************************************************

  template <class ImageT>
  bool is_transparent_helper( ImageT const& image, true_type ) {
    for( int32 y=0; y<image.rows(); ++y )
      for( int32 x=0; x<image.cols(); ++x )
        if( ! is_transparent(image(x,y)) ) return false;
    return true;
  }

  template <class ImageT>
  bool is_transparent_helper( ImageT const& /*image*/, false_type ) {
    return false;
  }

  /// Returns true if the given image is entirely transparent, or false if
  /// it is opaque or only partially transparent.
  template <class ImageT>
  bool is_transparent( ImageViewBase<ImageT> const& image ) {
    return is_transparent_helper( image.impl(), typename PixelHasAlpha<typename ImageT::pixel_type>::type() );
  }

  // *******************************************************************
  // image_blocks()
  // *******************************************************************

  /// A utility routine that, given an image, returns a vector of
  /// bounding boxes for sub-regions of the image of the specified
  /// size.  Note that bounding boxes along the right and bottom edges
  /// of the image will not have the specified dimension unless the
  /// image width and height are perfectly divisible by the bounding
  /// box width and height, respectively. This routine is useful if you
  /// want to apply an operation to a large image one region at a time.
  /// It will operate on any object that has cols() and rows() methods.
  inline std::vector<BBox2i>
  image_blocks(BBox2i const& object, int32 block_width, int32 block_height) {
    std::vector<BBox2i> bboxes;

    int32 j_offset = 0;
    while ( j_offset < object.height() ) {
      int32 j_dim = (object.height() - j_offset) < block_height ? (object.height() - j_offset) : block_height;
      int32 i_offset = 0;
      while ( i_offset < object.width() ) {
        int32 i_dim = (object.width() - i_offset) < block_width ? (object.width() - i_offset) : block_width;
        bboxes.push_back(BBox2i(i_offset + object.min().x(),
                                j_offset + object.min().y(),
                                i_dim,j_dim));
        i_offset += i_dim;
      }
      j_offset += j_dim;
    }
    return bboxes;
  }

  template <class T>
  inline std::vector<BBox2i>
  image_blocks(ImageViewBase<T> const& view, int32 block_width, int32 block_height ) {
    return image_blocks( bounding_box(view.impl()), block_width, block_height );
  }
  
  
/* Replaced by the BlobIndexThreaded class
  // ********************************************************************
  // blob_index()
  // ********************************************************************

  /// A utility that numbers off blobs of valid pixels. So it's
  /// important that you use pixel mask, or the entire image will be
  /// numbered one. This is for the most part a clone of the Matlab function bwlabel.
  template <class SourceT>
  void blob_index( ImageViewBase<SourceT> const& src,
                   ImageView<uint32>           & dst ) {

    if ( src.impl().planes() > 1 )
      vw_throw( NoImplErr() << "Blob index currently only works with 2D images." );
    dst.set_size( src.impl().cols(),
                  src.impl().rows() );
    fill(dst,0);

    // Initialize Graph used to pair blobs
    typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS> Graph;
    Graph connections;
    uint32 m_blob_count=1;

    { // Initial Pass
      typename SourceT::pixel_accessor s_acc = src.impl().origin();
      typename SourceT::pixel_accessor p_s_acc = src.impl().origin(); // previous
      typename ImageView<uint32>::pixel_accessor d_acc = dst.origin();
      typename ImageView<uint32>::pixel_accessor p_d_acc = dst.origin(); // previous

      // Top Corner
      if ( is_valid(*s_acc) ) {
        *d_acc = m_blob_count;
        m_blob_count++;
      }

      // Top Row
      s_acc.next_col();
      d_acc.next_col();
      for ( int32 i = 1; i < dst.cols(); i++ ) {
        if ( is_valid(*s_acc) ) {
          if ( is_valid(*p_s_acc) ) {
            *d_acc = *p_d_acc;
          } else {
            *d_acc = m_blob_count;
            m_blob_count++;
          }
          }
        s_acc.next_col();
        d_acc.next_col();
        p_s_acc.next_col();
        p_d_acc.next_col();
      }
    }

    { // Everything else (9 connected)
      typename SourceT::pixel_accessor s_acc_row = src.impl().origin();
      typename ImageView<uint32>::pixel_accessor d_acc_row = dst.origin();
      s_acc_row.advance(0,1);
      d_acc_row.advance(0,1);

      for (int j = dst.rows()-1; j; --j ) { // Not for indexing
        typename SourceT::pixel_accessor s_acc = s_acc_row;
        typename SourceT::pixel_accessor p_s_acc = s_acc_row;
        typename ImageView<uint32>::pixel_accessor d_acc = d_acc_row;
        typename ImageView<uint32>::pixel_accessor p_d_acc = d_acc_row;

        // Process
        for ( int i = dst.cols(); i; --i ) {
          if ( is_valid(*s_acc) ) {
            if ( i != dst.cols() ) {
              // Left
              p_s_acc.advance(-1,0);
              p_d_acc.advance(-1,0);
              if ( is_valid(*p_s_acc) ) {
                if ( (*d_acc != 0) && (*d_acc != *p_d_acc) ) {
                  boost::add_edge(*p_d_acc,*d_acc,connections);
                } else
                  *d_acc = *p_d_acc;
              }
              // Upper Left
              p_s_acc.advance(0,-1);
              p_d_acc.advance(0,-1);
              if ( is_valid(*p_s_acc) ) {
                if ( (*d_acc != 0) && (*d_acc != *p_d_acc) ) {
                  boost::add_edge(*p_d_acc,*d_acc,connections);
                } else
                  *d_acc = *p_d_acc;
              }
            } else {
              p_s_acc.advance(-1,-1);
              p_d_acc.advance(-1,-1);
            }
            // Upper
            p_s_acc.advance(1,0);
            p_d_acc.advance(1,0);
            if ( is_valid(*p_s_acc) ) {
              if ( (*d_acc != 0) && (*d_acc != *p_d_acc) ) {
                boost::add_edge(*p_d_acc,*d_acc,connections);
              } else
                *d_acc = *p_d_acc;
            }
            // Upper Right
            p_s_acc.advance(1,0);
            p_d_acc.advance(1,0);
            if ( i != 1 )
              if ( is_valid(*p_s_acc) ) {
                if ( (*d_acc != 0) && (*d_acc != *p_d_acc) ) {
                  boost::add_edge(*p_d_acc,*d_acc,connections);
                } else
                  *d_acc = *p_d_acc;
              }
            // Setting if not
            p_s_acc.advance(-1,1);
            p_d_acc.advance(-1,1);
            if ( *d_acc == 0 ) {
              *d_acc = m_blob_count;
              m_blob_count++;
            }
          }
          s_acc.next_col();
          p_s_acc.next_col();
          d_acc.next_col();
          p_d_acc.next_col();
        } // end row process

        s_acc_row.next_row();
        d_acc_row.next_row();
      }
    }

    // Making sure connections has vertices for all indexes made
    add_edge(m_blob_count-1,m_blob_count-1,connections);
    std::vector<uint32> component(boost::num_vertices(connections));
    m_blob_count = boost::connected_components(connections, &component[0])-1;

    { // Update index map to optimal numbering
      ImageView<uint32>::pixel_accessor p_d_acc = dst.origin(); // previous

      for ( int32 r = 0; r < dst.rows(); r++ ) {
        ImageView<uint32>::pixel_accessor d_acc = p_d_acc;
        for ( int32 c = 0; c < dst.cols(); c++ ) {
          if ( (*d_acc) != 0 ) {

            *d_acc = component[*d_acc];
          }
          d_acc.next_col();
        }
        p_d_acc.next_row();
      }
    }
  }

  // Simple interface
  template <class SourceT>
  ImageView<uint32> blob_index( ImageViewBase<SourceT> const& src ) {
    ImageView<uint32> result( src.impl().cols(), src.impl().rows() );
    blob_index( src, result );
    return result;
  }
*/

}

#endif//__VW_IMAGE_ALGORITHMS_FUNCTIONS_H__
