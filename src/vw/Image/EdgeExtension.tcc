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
  // The edge extension types
  // *******************************************************************

  // Here we define the edge extension modes as polymorphic functors.
  // You can extend the list of supported edge extension modes by
  // creating a new functor type derived from EdgeExtensionBase and
  // implementing the function call operator, as shown in the
  // examples below.

  /// A special type providing no edge extension.  This is mainly
  /// used as a signal to certain types to adopt fundamentally
  /// different behavior.
  struct NoEdgeExtension : EdgeExtensionBase {
    template <class ViewT>
    inline typename ViewT::pixel_type operator()( const ViewT &view, int32 i, int32 j, int32 p=0) const {
      return view(i,j,p);
    }
    template <class ViewT>
    inline BBox2i source_bbox( ViewT const& /*view*/, BBox2i const& bbox ) const {
      return bbox;
    }
  };

  /// An edge extention type that extends the image with zeroes in
  /// all directions.
  struct ZeroEdgeExtension : EdgeExtensionBase {
    template <class ViewT>
    inline typename ViewT::pixel_type operator()( const ViewT &view, int32 i, int32 j, int32 p=0) const {
      if( i>=0 && j>=0 && i<view.cols() && j<view.rows() )
        return view(i,j,p);
      else
        return typename ViewT::pixel_type();
    }
    template <class ViewT>
    inline BBox2i source_bbox( ViewT const& view, BBox2i const& bbox ) const {
      BBox2i result = bbox;
      result.crop( BBox2i( 0, 0, view.cols(), view.rows() ) );
      return result;
    }
  };

  /// An edge extention type that extends the image with a
  /// user-supplied value in all directions.
  template <class PixelT>
  struct ValueEdgeExtension : EdgeExtensionBase {
    PixelT m_pix;
    ValueEdgeExtension(PixelT pix) : m_pix(pix) {}

    template <class ViewT>
    inline PixelT operator()( const ViewT &view, int32 i, int32 j, int32 p=0) const {
      if( i>=0 && j>=0 && i<view.cols() && j<view.rows() )
        return view(i,j,p);
      else
        return m_pix;
    }

    template <class ViewT>
    inline BBox2i source_bbox( ViewT const& view, BBox2i const& bbox ) const {
      BBox2i result = bbox;
      result.crop( BBox2i( 0, 0, view.cols(), view.rows() ) );
      return result;
    }
  };

  /// An edge extension type that extends the image using constant
  /// functions in all directions.  In other words, it returns the
  /// nearest valid pixel.
  struct ConstantEdgeExtension : EdgeExtensionBase {
    template <class ViewT>
    inline typename ViewT::pixel_type operator()( const ViewT &view, int32 i, int32 j, int32 p=0) const {
      return view((i<0) ? 0 : (i>=view.cols()) ? (view.cols()-1) : i,
                  (j<0) ? 0 : (j>=view.rows()) ? (view.rows()-1) : j, p);
    }
    template <class ViewT>
    inline BBox2i source_bbox( ViewT const& view, BBox2i const& bbox ) const {
      BBox2i result = bbox;
      if( bbox.min().x() < 0 ) result.min().x() = 0;
      else if( bbox.min().x() >= view.cols() ) result.min().x() = view.cols()-1;
      if( bbox.min().y() < 0 ) result.min().y() = 0;
      else if( bbox.min().y() >= view.rows() ) result.min().y() = view.rows()-1;
      if( bbox.max().x() > view.cols() ) result.max().x() = view.cols();
      else if( bbox.max().x() <= 0 ) result.max().x() = 1;
      if( bbox.max().y() > view.rows() ) result.max().y() = view.rows();
      else if( bbox.max().y() <= 0 ) result.max().y() = 1;
      return result;
    }
  };

  /// A periodic edge extension type.
  struct PeriodicEdgeExtension : EdgeExtensionBase {
    template <class ViewT>
    inline typename ViewT::pixel_type operator()( const ViewT &view, int32 i, int32 j, int32 p=0) const {
      int32 d_i=i, d_j=j;
      d_i %= int(view.cols());
      if( d_i < 0 ) d_i += view.cols();
      d_j %= int(view.rows());
      if( d_j < 0 ) d_j += view.rows();
      return view(d_i,d_j,p);
    }
    template <class ViewT>
    inline BBox2i source_bbox( ViewT const& view, BBox2i const& bbox ) const {
      BBox2i result;
      if( bbox.width() >= view.cols() ) {
        result.min().x() = 0;
        result.max().x() = view.cols();
      }
      else {
        result.min().x() = mod(bbox.min().x(), view.cols());
        result.max().x() = mod(bbox.max().x()-1, view.cols())+1;
        if( result.min().x() >= result.max().x() ) {
          result.min().x() = 0;
          result.max().x() = view.cols();
        }
      }
      if( bbox.height() >= view.rows() ) {
        result.min().y() = 0;
        result.max().y() = view.rows();
      }
      else {
        result.min().y() = mod(bbox.min().y(), view.rows());
        result.max().y() = mod(bbox.max().y()-1, view.rows())+1;
        if( result.min().y() >= result.max().y() ) {
          result.min().y() = 0;
          result.max().y() = view.rows();
        }
      }
      return result;
    }
  };

  /// A cylindrical edge extension type: periodic in the x axis, constant in the y axis.
  struct CylindricalEdgeExtension : EdgeExtensionBase {
    template <class ViewT>
    inline typename ViewT::pixel_type operator()( const ViewT &view, int32 i, int32 j, int32 p=0) const {
      int32 d_i=i;
      d_i %= int(view.cols());
      if( d_i < 0 ) d_i += view.cols();
      int32 d_j = (j<0) ? 0 : (j>=view.rows()) ? (view.rows()-1) : j;
      return view(d_i,d_j,p);
    }
    template <class ViewT>
    inline BBox2i source_bbox( ViewT const& view, BBox2i const& bbox ) const {
      BBox2i result;
      if( bbox.width() >= view.cols() ) {
        result.min().x() = 0;
        result.max().x() = view.cols();
      }
      else {
        result.min().x() = mod(bbox.min().x(), view.cols());
        result.max().x() = mod(bbox.max().x()-1, view.cols())+1;
        if( result.min().x() >= result.max().x() ) {
          result.min().x() = 0;
          result.max().x() = view.cols();
        }
      }
      if( bbox.min().y() < 0 ) result.min().y() = 0;
      else if( bbox.min().y() >= view.rows() ) result.min().y() = view.rows()-1;
      else result.min().y() = bbox.min().y();
      if( bbox.max().y() > view.rows() ) result.max().y() = view.rows();
      else if( bbox.max().y() <= 0 ) result.max().y() = 1;
      else result.max().y() = bbox.max().y();
      return result;
    }
  };

  /// A reflection edge extension type.
  struct ReflectEdgeExtension : EdgeExtensionBase {
    template <class ViewT>
    inline typename ViewT::pixel_type operator()( const ViewT &view, int32 i, int32 j, int32 p=0) const {
      int32 d_i=i, d_j=j;
      if( d_i < 0 ) d_i = -d_i;
      int32 vcm1 = view.cols() - 1;
      d_i %= 2*vcm1;
      if( d_i > vcm1 ) d_i = 2*vcm1 - d_i;
      if( d_j<0 ) d_j=-d_j;
      int32 vrm1 = view.rows() - 1;
      d_j %= 2*vrm1;
      if( d_j > vrm1 ) d_j = 2*vrm1 - d_j;
      return view(d_i,d_j,p);
    }
    template <class ViewT>
    inline BBox2i source_bbox( ViewT const& view, BBox2i const& bbox ) const {
      BBox2i result;
      if( bbox.width() >= 2*view.cols()-2 ) {
        result.min().x() = 0;
        result.max().x() = view.cols();
      }
      else {
        result.min().x() = mod(bbox.min().x(), 2*view.cols()-2);
        result.max().x() = mod(bbox.max().x()-1, 2*view.cols()-2)+1;
        if( result.min().x() >= result.max().x() ) {
          result.max().x() = (std::min)((std::max)(result.max().x(),2*view.cols()-1-result.min().x()),view.cols());
          result.min().x() = 0;
        }
        else if( result.min().x() >= view.cols() ) {
          int tmp = 2*view.cols()-1-result.min().x();
          result.min().x() = 2*view.cols()-1-result.max().x();
          result.max().x() = tmp;
        }
        else if( result.max().x() > view.cols() ) {
          result.min().x() = (std::min)(result.min().x(),2*view.cols()-1-result.max().x());
          result.max().x() = view.cols();
        }
      }
      if( bbox.height() >= 2*view.rows()-2 ) {
        result.min().y() = 0;
        result.max().y() = view.rows();
      }
      else {
        result.min().y() = mod(bbox.min().y(), 2*view.rows()-2);
        result.max().y() = mod(bbox.max().y()-1, 2*view.rows()-2)+1;
        if( result.min().y() >= result.max().y() ) {
          result.max().y() = (std::min)((std::max)(result.max().y(),2*view.rows()-1-result.min().y()),view.rows());
          result.min().y() = 0;
        }
        else if( result.min().y() >= view.rows() ) {
          int tmp = 2*view.rows()-1-result.min().y();
          result.min().y() = 2*view.rows()-1-result.max().y();
          result.max().y() = tmp;
        }
        else if( result.max().y() > view.rows() ) {
          result.min().y() = (std::min)(result.min().y(),2*view.rows()-1-result.max().y());
          result.max().y() = view.rows();
        }
      }
      return result;
    }
  };

  /// A linear extrapolation edge extension type.
  struct LinearEdgeExtension : EdgeExtensionBase {
    template <class ViewT>
    inline typename ViewT::pixel_type operator()( const ViewT &view, int32 i, int32 j, int32 p=0) const {
      int32 vcm1 = view.cols() - 1;
      int32 vrm1 = view.rows() - 1;
      if( i < 0 ) {
        if( j < 0 ) return view(0,0,p) - i*(view(0,0,p)-view(1,0,p)) - j*(view(0,0,p)-view(0,1,p));
        else if( j > vrm1 ) return view(0,vrm1,p) - i*(view(0,vrm1,p)-view(1,vrm1,p)) + (j-vrm1)*(view(0,vrm1,p)-view(0,vrm1-1,p));
        else return view(0,j,p) - i*(view(0,j,p)-view(1,j,p));
      }
      else if( i > vcm1 ) {
        if( j < 0 ) return view(vcm1,0,p) + (i-vcm1)*(view(vcm1,0,p)-view(vcm1-1,0,p)) - j*(view(vcm1,0,p)-view(vcm1,1,p));
        else if( j > vrm1 ) return view(vcm1,vrm1,p) + (i-vcm1)*(view(vcm1,vrm1,p)-view(vcm1-1,vrm1,p)) + (j-vrm1)*(view(vcm1,vrm1,p)-view(vcm1,vrm1-1,p));
        else return view(vcm1,j,p) + (i-vcm1)*(view(vcm1,j,p)-view(vcm1-1,j,p));
      }
      else {
        if( j < 0 ) return view(i,0,p) - j*(view(i,0,p)-view(i,1,p));
        else if( j > vrm1 ) return view(i,vrm1,p) + (j-vrm1)*(view(i,vrm1,p)-view(i,vrm1-1,p));
        else return view(i,j,p);
      }
    }
    template <class ViewT>
    inline BBox2i source_bbox( ViewT const& view, BBox2i const& bbox ) const {
      BBox2i result = bbox;
      if( bbox.min().x() < 0 ) {
        result.min().x() = 0;
        result.max().x() = (std::max)( result.max().x(), 2 );
      }
      if( bbox.max().x() > view.cols() ) {
        result.max().x() = view.cols();
        result.min().x() = (std::min)( result.min().x(), view.cols()-2 );
      }
      if( bbox.min().y() < 0 ) {
        result.min().y() = 0;
        result.max().y() = (std::max)( result.max().y(), 2 );
      }
      if( bbox.max().y() > view.rows() ) {
        result.max().y() = view.rows();
        result.min().y() = (std::min)( result.min().y(), view.rows()-2 );
      }
      return result;
    }
  };

} // namespace vw

