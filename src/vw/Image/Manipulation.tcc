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
// copy()
// *******************************************************************

// Class definition
template <class ImageT>
class CopyView : public ImageViewBase<CopyView<ImageT> >
{
private:
  ImageView<typename ImageT::pixel_type> m_child;
public:
  typedef typename ImageView<typename ImageT::pixel_type>::pixel_type pixel_type;
  typedef pixel_type const& result_type;
  typedef typename ImageView<typename ImageT::pixel_type>::pixel_accessor pixel_accessor;

  CopyView( ImageT const& image ) : m_child(image.cols(),image.rows(),image.planes()) {
    image.rasterize( m_child, BBox2i(0,0,image.cols(),image.rows()) );
  }

  inline int32 cols() const { return m_child.cols(); }
  inline int32 rows() const { return m_child.rows(); }
  inline int32 planes() const { return m_child.planes(); }

  inline pixel_accessor origin() const { return m_child.origin(); }
  inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { return m_child(i,j,p); }

  typedef CopyView prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& /*bbox*/ ) const { return *this; }
  template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { m_child.rasterize( dest, bbox ); }
};


// *******************************************************************
// Transpose
// *******************************************************************

// Specialized pixel accessor
template <class ChildT>
class TransposePixelAccessor
{
private:
  ChildT m_child;
public:
  typedef typename ChildT::pixel_type pixel_type;
  typedef typename ChildT::result_type result_type;
  typedef typename ChildT::offset_type offset_type;
  TransposePixelAccessor( ChildT const& acc ) : m_child(acc) {}

  inline TransposePixelAccessor& next_col() { m_child.next_row(); return *this; }
  inline TransposePixelAccessor& prev_col() { m_child.prev_row(); return *this; }
  inline TransposePixelAccessor& next_row() { m_child.next_col(); return *this; }
  inline TransposePixelAccessor& prev_row() { m_child.prev_col(); return *this; }
  inline TransposePixelAccessor& next_plane() { m_child.next_plane(); return *this; }
  inline TransposePixelAccessor& prev_plane() { m_child.prev_plane(); return *this; }
  inline TransposePixelAccessor& advance( offset_type di, offset_type dj, ssize_t dp=0 ) { m_child.advance(dj,di,dp); return *this; }

  inline result_type operator*() const { return *m_child; }
};

// Class definition
template <class ImageT>
class TransposeView : public ImageViewBase<TransposeView<ImageT> >
{
private:
  ImageT m_child;
public:

  typedef typename ImageT::pixel_type pixel_type;
  typedef typename ImageT::result_type result_type;
  typedef TransposePixelAccessor<typename ImageT::pixel_accessor> pixel_accessor;

  TransposeView( ImageT const& image ) : m_child(image) {}

  inline int32 cols() const { return m_child.rows(); }
  inline int32 rows() const { return m_child.cols(); }
  inline int32 planes() const { return m_child.planes(); }

  inline pixel_accessor origin() const { return m_child.origin(); }

  inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { return m_child(j,i,p); }

  template <class ViewT>
  TransposeView const& operator=( ImageViewBase<ViewT> const& view ) const {
    view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
    return *this;
  }

  template <class ViewT>
  TransposeView& operator=( ImageViewBase<ViewT> const& view ) {
    *const_cast<const TransposeView*>(this) = view.impl();
    return *this;
  }

  ImageT const& child() const {
    return m_child;
  }

  /// \cond INTERNAL
  typedef TransposeView<typename ImageT::prerasterize_type> prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
    BBox2i child_bbox( bbox.min().y(), bbox.min().x(), bbox.height(), bbox.width() );
    return prerasterize_type( m_child.prerasterize(child_bbox) );
  }
  template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
  /// \endcond
};


// *******************************************************************
// Rotate180
// *******************************************************************

// Specialized pixel accessor
template <class ChildT>
class Rotate180PixelAccessor
{
private:
  ChildT m_child;
public:
  typedef typename ChildT::pixel_type  pixel_type;
  typedef typename ChildT::result_type result_type;
  typedef typename ChildT::offset_type offset_type;
  Rotate180PixelAccessor( ChildT const& image ) : m_child(image) {}

  inline Rotate180PixelAccessor& next_col  () { m_child.prev_col  (); return *this; }
  inline Rotate180PixelAccessor& prev_col  () { m_child.next_col  (); return *this; }
  inline Rotate180PixelAccessor& next_row  () { m_child.prev_row  (); return *this; }
  inline Rotate180PixelAccessor& prev_row  () { m_child.next_row  (); return *this; }
  inline Rotate180PixelAccessor& next_plane() { m_child.next_plane(); return *this; }
  inline Rotate180PixelAccessor& prev_plane() { m_child.prev_plane(); return *this; }
  inline Rotate180PixelAccessor& advance( offset_type di, offset_type dj, ssize_t dp=0 ) { m_child.advance(-di,-dj,dp); return *this; }

  inline result_type operator*() const { return *m_child; }
};

// Image View Class
template <class ImageT>
class Rotate180View : public ImageViewBase<Rotate180View<ImageT> >
{
private:
  ImageT m_child;
public:

  typedef typename ImageT::pixel_type   pixel_type;
  typedef typename ImageT::result_type result_type;
  typedef Rotate180PixelAccessor<typename ImageT::pixel_accessor> pixel_accessor;

  Rotate180View( ImageT const& image ) : m_child(image) {}

  inline int32 cols  () const { return m_child.cols  (); }
  inline int32 rows  () const { return m_child.rows  (); }
  inline int32 planes() const { return m_child.planes(); }

  inline pixel_accessor origin() const { return m_child.origin().advance(cols()-1,rows()-1); }

  inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { return m_child(cols()-1-i,rows()-1-j,p); }

  template <class ViewT>
  Rotate180View const& operator=( ImageViewBase<ViewT> const& view ) const {
    view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
    return *this;
  }

  template <class ViewT>
  Rotate180View& operator=( ImageViewBase<ViewT> const& view ) {
    *const_cast<const Rotate180View*>(this) = view.impl();
    return *this;
  }

  ImageT const& child() const {
    return m_child;
  }

  /// \cond INTERNAL
  typedef Rotate180View<typename ImageT::prerasterize_type> prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
    BBox2i child_bbox( cols()-bbox.max().x(), rows()-bbox.max().y(), bbox.width(), bbox.height() );
    return prerasterize_type( m_child.prerasterize(child_bbox) );
  }
  template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { 
    vw::rasterize( prerasterize(bbox), dest, bbox ); 
  }
  /// \endcond
};


// *******************************************************************
// Rotate90CW
// *******************************************************************

// Specialized pixel accessor
template <class ChildT>
class Rotate90CWPixelAccessor
{
private:
  ChildT m_child;
public:
  typedef typename ChildT::pixel_type pixel_type;
  typedef typename ChildT::result_type result_type;
  typedef typename ChildT::offset_type offset_type;

  Rotate90CWPixelAccessor( ChildT const& acc ) : m_child(acc) {}
  inline Rotate90CWPixelAccessor& next_col() { m_child.prev_row(); return *this; }
  inline Rotate90CWPixelAccessor& prev_col() { m_child.next_row(); return *this; }
  inline Rotate90CWPixelAccessor& next_row() { m_child.next_col(); return *this; }
  inline Rotate90CWPixelAccessor& prev_row() { m_child.prev_col(); return *this; }
  inline Rotate90CWPixelAccessor& next_plane() { m_child.next_plane(); return *this; }
  inline Rotate90CWPixelAccessor& prev_plane() { m_child.prev_plane(); return *this; }
  inline Rotate90CWPixelAccessor& advance( offset_type di, offset_type dj, ssize_t dp=0 ) { m_child.advance(dj,-di,dp); return *this; }

  inline result_type operator*() const { return *m_child; }
};

// Class definition
template <class ImageT>
class Rotate90CWView : public ImageViewBase<Rotate90CWView<ImageT> >
{
private:
  ImageT m_child;
public:

  typedef typename ImageT::pixel_type pixel_type;
  typedef typename ImageT::result_type result_type;
  typedef Rotate90CWPixelAccessor<typename ImageT::pixel_accessor> pixel_accessor;

  Rotate90CWView( ImageT const& image ) : m_child(image) {}

  inline int32 cols() const { return m_child.rows(); }
  inline int32 rows() const { return m_child.cols(); }
  inline int32 planes() const { return m_child.planes(); }

  inline pixel_accessor origin() const { return m_child.origin().advance(0,cols()-1); }

  inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { return m_child(j,cols()-1-i,p); }

  template <class ViewT>
  Rotate90CWView const& operator=( ImageViewBase<ViewT> const& view ) const {
    view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
    return *this;
  }

  template <class ViewT>
  Rotate90CWView& operator=( ImageViewBase<ViewT> const& view ) {
    *const_cast<const Rotate90CWView*>(this) = view.impl();
    return *this;
  }

  ImageT const& child() const {
    return m_child;
  }

  /// \cond INTERNAL
  typedef Rotate90CWView<typename ImageT::prerasterize_type> prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
    BBox2i child_bbox( bbox.min().y(), cols()-bbox.max().x(), bbox.height(), bbox.width() );
    return prerasterize_type( m_child.prerasterize(child_bbox) );
  }
  template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
  /// \endcond
};


// *******************************************************************
// Rotate90CCW
// *******************************************************************

// Specialized pixel accessor
template <class ChildT>
class Rotate90CCWPixelAccessor
{
private:
  ChildT m_child;
public:
  typedef typename ChildT::pixel_type pixel_type;
  typedef typename ChildT::result_type result_type;
  typedef typename ChildT::offset_type offset_type;
  Rotate90CCWPixelAccessor( ChildT const& acc ) : m_child(acc) {}

  inline Rotate90CCWPixelAccessor& next_col() { m_child.next_row(); return *this; }
  inline Rotate90CCWPixelAccessor& prev_col() { m_child.prev_row(); return *this; }
  inline Rotate90CCWPixelAccessor& next_row() { m_child.prev_col(); return *this; }
  inline Rotate90CCWPixelAccessor& prev_row() { m_child.next_col(); return *this; }
  inline Rotate90CCWPixelAccessor& next_plane() { m_child.next_plane(); return *this; }
  inline Rotate90CCWPixelAccessor& prev_plane() { m_child.prev_plane(); return *this; }
  inline Rotate90CCWPixelAccessor& advance( offset_type di, offset_type dj, ssize_t dp=0 ) { m_child.advance(-dj,di,dp); return *this; }

  inline result_type operator*() const { return *m_child; }
};

// Class definition
template <class ImageT>
class Rotate90CCWView : public ImageViewBase<Rotate90CCWView<ImageT> >
{
private:
  ImageT m_child;
public:

  typedef typename ImageT::pixel_type pixel_type;
  typedef typename ImageT::result_type result_type;
  typedef Rotate90CCWPixelAccessor<typename ImageT::pixel_accessor> pixel_accessor;

  Rotate90CCWView( ImageT const& image ) : m_child(image) {}

  inline int32 cols() const { return m_child.rows(); }
  inline int32 rows() const { return m_child.cols(); }
  inline int32 planes() const { return m_child.planes(); }

  inline pixel_accessor origin() const { return m_child.origin().advance(rows()-1,0); }

  inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { return m_child(rows()-1-j,i,p); }

  template <class ViewT>
  Rotate90CCWView const& operator=( ImageViewBase<ViewT> const& view ) const {
    view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
    return *this;
  }

  template <class ViewT>
  Rotate90CCWView& operator=( ImageViewBase<ViewT> const& view ) {
    *const_cast<const Rotate90CCWView*>(this) = view.impl();
    return *this;
  }

  ImageT const& child() const {
    return m_child;
  }

  /// \cond INTERNAL
  typedef Rotate90CCWView<typename ImageT::prerasterize_type> prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
    BBox2i child_bbox( rows()-bbox.max().y(), bbox.min().x(), bbox.height(), bbox.width() );
    return prerasterize_type( m_child.prerasterize(child_bbox) );
  }
  template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
  /// \endcond
};


// *******************************************************************
// FlipVertical
// *******************************************************************

// Specialized pixel accessor
template <class ChildT>
class FlipVerticalPixelAccessor
{
private:
  ChildT m_child;
public:
  typedef typename ChildT::pixel_type pixel_type;
  typedef typename ChildT::result_type result_type;
  typedef typename ChildT::offset_type offset_type;
  FlipVerticalPixelAccessor( ChildT const& acc ) : m_child(acc) {}

  inline FlipVerticalPixelAccessor& next_col() { m_child.next_col(); return *this; }
  inline FlipVerticalPixelAccessor& prev_col() { m_child.prev_col(); return *this; }
  inline FlipVerticalPixelAccessor& next_row() { m_child.prev_row(); return *this; }
  inline FlipVerticalPixelAccessor& prev_row() { m_child.next_row(); return *this; }
  inline FlipVerticalPixelAccessor& next_plane() { m_child.next_plane(); return *this; }
  inline FlipVerticalPixelAccessor& prev_plane() { m_child.prev_plane(); return *this; }
  inline FlipVerticalPixelAccessor& advance( offset_type di, offset_type dj, ssize_t dp=0 ) { m_child.advance(di,-dj,dp); return *this; }

  inline result_type operator*() const { return *m_child; }
};

// Class definition
template <class ImageT>
class FlipVerticalView : public ImageViewBase<FlipVerticalView<ImageT> >
{
private:
  ImageT m_child;
public:

  typedef typename ImageT::pixel_type pixel_type;
  typedef typename ImageT::result_type result_type;
  typedef FlipVerticalPixelAccessor<typename ImageT::pixel_accessor> pixel_accessor;

  FlipVerticalView( ImageT const& image ) : m_child(image) {}

  inline int32 cols() const { return m_child.cols(); }
  inline int32 rows() const { return m_child.rows(); }
  inline int32 planes() const { return m_child.planes(); }

  inline pixel_accessor origin() const { return m_child.origin().advance(0,rows()-1); }

  inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { return m_child(i,rows()-1-j,p); }

  ImageT const& child() const {
    return m_child;
  }

  template <class ViewT>
  FlipVerticalView const& operator=( ImageViewBase<ViewT> const& view ) const {
    view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
    return *this;
  }

  template <class ViewT>
  FlipVerticalView& operator=( ImageViewBase<ViewT> const& view ) {
    *const_cast<const FlipVerticalView*>(this) = view.impl();
    return *this;
  }

  /// \cond INTERNAL
  typedef FlipVerticalView<typename ImageT::prerasterize_type> prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
    BBox2i child_bbox( bbox.min().x(), rows()-bbox.max().y(), bbox.width(), bbox.height() );
    return prerasterize_type( m_child.prerasterize(child_bbox) );
  }
  template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
  /// \endcond
};


// *******************************************************************
// FlipHorizontal
// *******************************************************************

// Specialized pixel accessor
template <class ChildT>
class FlipHorizontalPixelAccessor
{
private:
  ChildT m_child;
public:
  typedef typename ChildT::pixel_type pixel_type;
  typedef typename ChildT::result_type result_type;
  typedef typename ChildT::offset_type offset_type;
  FlipHorizontalPixelAccessor( ChildT const& acc ) : m_child(acc) {}

  inline FlipHorizontalPixelAccessor& next_col() { m_child.prev_col(); return *this; }
  inline FlipHorizontalPixelAccessor& prev_col() { m_child.next_col(); return *this; }
  inline FlipHorizontalPixelAccessor& next_row() { m_child.next_row(); return *this; }
  inline FlipHorizontalPixelAccessor& prev_row() { m_child.prev_row(); return *this; }
  inline FlipHorizontalPixelAccessor& next_plane() { m_child.next_plane(); return *this; }
  inline FlipHorizontalPixelAccessor& prev_plane() { m_child.prev_plane(); return *this; }
  inline FlipHorizontalPixelAccessor& advance( offset_type di, offset_type dj, ssize_t dp=0 ) { m_child.advance(-di,dj,dp); return *this; }

  inline result_type operator*() const { return *m_child; }
};

// Class definition
template <class ImageT>
class FlipHorizontalView : public ImageViewBase<FlipHorizontalView<ImageT> >
{
private:
  ImageT m_child;
public:

  typedef typename ImageT::pixel_type pixel_type;
  typedef typename ImageT::result_type result_type;
  typedef FlipHorizontalPixelAccessor<typename ImageT::pixel_accessor> pixel_accessor;

  FlipHorizontalView( ImageT const& image ) : m_child(image) {}

  inline int32 cols() const { return m_child.cols(); }
  inline int32 rows() const { return m_child.rows(); }
  inline int32 planes() const { return m_child.planes(); }

  inline pixel_accessor origin() const { return m_child.origin().advance(cols()-1,0); }

  inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { return m_child(cols()-1-i,j,p); }

  ImageT const& child() const {
    return m_child;
  }

  template <class ViewT>
  FlipHorizontalView const& operator=( ImageViewBase<ViewT> const& view ) const {
    view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
    return *this;
  }

  template <class ViewT>
  FlipHorizontalView& operator=( ImageViewBase<ViewT> const& view ) {
    *const_cast<const FlipHorizontalView*>(this) = view.impl();
    return *this;
  }

  /// \cond INTERNAL
  typedef FlipHorizontalView<typename ImageT::prerasterize_type> prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
    BBox2i child_bbox( cols()-bbox.max().x(), bbox.min().y(), bbox.width(), bbox.height() );
    return prerasterize_type( m_child.prerasterize(child_bbox) );
  }
  template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
  /// \endcond
};

// *******************************************************************
// crop()
// *******************************************************************

// Class definition
template <class ImageT>
class CropView : public ImageViewBase< CropView<ImageT> > {
private:
  typedef typename boost::mpl::if_<IsFloatingPointIndexable<ImageT>, double, int32>::type offset_type;

  ImageT m_child;
  offset_type m_ci, m_cj;
  int32 m_di, m_dj; // Cropped width and height

public:
  typedef typename ImageT::pixel_type pixel_type;
  typedef typename ImageT::result_type result_type;
  typedef typename ImageT::pixel_accessor pixel_accessor;

  CropView( ImageT const& image, offset_type const upper_left_i, offset_type const upper_left_j, int32 const width, int32 const height ) :
    m_child(image), m_ci(upper_left_i), m_cj(upper_left_j), m_di(width), m_dj(height) {}

  template<class RealT>
  CropView( ImageT const& image, BBox<RealT,2> const& bbox) :
    m_child(image),
    m_ci((offset_type)(bbox.min()[0])),
    m_cj((offset_type)(bbox.min()[1])),
    m_di(int32(.5+(bbox.width()))),
    m_dj(int32(.5+(bbox.height()))) {}

  inline int32 cols() const { return m_di; }
  inline int32 rows() const { return m_dj; }
  inline int32 planes() const { return m_child.planes(); }

  inline pixel_accessor origin() const { return m_child.origin().advance(m_ci, m_cj); }

  inline result_type operator()( offset_type i, offset_type j, int32 p=0 ) const { return m_child(m_ci + i, m_cj + j, p); }

  CropView const& operator=( CropView const& view ) const {
    view.rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
    return *this;
  }

  template <class ViewT>
  CropView const& operator=( ImageViewBase<ViewT> const& view ) const {
    view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
    return *this;
  }

  ImageT const& child() const {
    return m_child;
  }

  /// \cond INTERNAL
  typedef CropView<typename ImageT::prerasterize_type> prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
    return prerasterize_type( m_child.prerasterize(bbox+Vector2i(m_ci,m_cj)), m_ci, m_cj, m_di, m_dj );
  }
  template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
    // FIXME Warning: This does not respect floating-point offsets!
    vw::rasterize( prerasterize(bbox), dest, bbox );
  }
  /// \endcond
};

// *******************************************************************
// subsample()
// *******************************************************************

// Specialized image accessor
template <class ChildT>
class SubsamplePixelAccessor {
  ChildT m_child;
  int32 m_xdelta, m_ydelta;
public:
  typedef typename ChildT::pixel_type pixel_type;
  typedef typename ChildT::result_type result_type;
  typedef typename ChildT::offset_type offset_type;
  SubsamplePixelAccessor( ChildT const& acc , int32 xdelta, int32 ydelta) : m_child(acc), m_xdelta(xdelta), m_ydelta(ydelta) {}

  inline SubsamplePixelAccessor& next_col() { m_child.advance(  m_xdelta, 0 ); return *this; }
  inline SubsamplePixelAccessor& prev_col() { m_child.advance( -m_xdelta, 0 ); return *this; }
  inline SubsamplePixelAccessor& next_row() { m_child.advance( 0,  m_ydelta ); return *this; }
  inline SubsamplePixelAccessor& prev_row() { m_child.advance( 0, -m_ydelta ); return *this; }
  inline SubsamplePixelAccessor& next_plane() { m_child.next_plane(); return *this; }
  inline SubsamplePixelAccessor& prev_plane() { m_child.prev_plane(); return *this; }
  inline SubsamplePixelAccessor& advance( offset_type di, offset_type dj, ssize_t dp=0 ) {
    m_child.advance((offset_type)m_xdelta*di,(offset_type)m_ydelta*dj,dp);
    return *this;
  }

  inline result_type operator*() const { return *m_child; }
};

// Class definition
template <class ImageT>
class SubsampleView : public ImageViewBase<SubsampleView<ImageT> > {
  ImageT m_child;
  int32 m_xdelta, m_ydelta;
public:
  typedef typename ImageT::pixel_type pixel_type;
  typedef typename ImageT::result_type result_type;
  typedef SubsamplePixelAccessor<typename ImageT::pixel_accessor> pixel_accessor;

  SubsampleView( ImageT const& image, int32 subsampling_factor ) :
    m_child(image), m_xdelta(subsampling_factor), m_ydelta(subsampling_factor) {
    VW_ASSERT( m_xdelta > 0 && m_ydelta > 0,
               ArgumentErr() << "SubsampleView: Arguments must be greater than zero." );
  }
  SubsampleView( ImageT const& image, int32 xfactor, int32 yfactor ) :
    m_child(image), m_xdelta(xfactor), m_ydelta(yfactor) {
    VW_ASSERT( m_xdelta > 0 && m_ydelta > 0,
               ArgumentErr() << "SubsampleView: Arguments must be greater than zero." );
  }

  inline int32 cols() const { return 1 + (m_child.cols()-1)/m_xdelta; }
  inline int32 rows() const { return 1 + (m_child.rows()-1)/m_ydelta; }
  inline int32 planes() const { return m_child.planes(); }

  inline pixel_accessor origin() const { return pixel_accessor(m_child.origin(), m_xdelta, m_ydelta); }
  inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { return m_child(m_xdelta*i,m_ydelta*j,p); }

  ImageT const& child() const { return m_child; }

  /// \cond INTERNAL
  // This complicated prerasterize call is to reduce the overall
  // memory requirement in the event the user is greatly reducing
  // the input.
  typedef CropView<ImageView<pixel_type> > prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
    ImageView<pixel_type> buffer( bbox.width(), bbox.height() );

    // The size of the chunks that we will be rastering one at a time.
    Vector2i sub_region_size( bbox.width () / m_xdelta,
                              bbox.height() / m_ydelta );
    if ( sub_region_size.x() < 2 )
      sub_region_size.x() = 2;
    if ( sub_region_size.y() < 2 )
      sub_region_size.y() = 2;

    typedef std::vector<BBox2i> ContainerT;
    ContainerT bboxes = subdivide_bbox( bbox, sub_region_size.x(), sub_region_size.y() );

    typedef SubsampleView<typename ImageT::prerasterize_type> input_pre_type;

    for ( ContainerT::const_iterator b = bboxes.begin(); b != bboxes.end(); ++b ) {
      // The math for the bbox may seem weird. It's not just the
      // size of the bbox scaled up. It's the minimum width we need
      // to achieve the output samples we want. I used this method
      // to avoid having to put conditionals in the math .. as
      // otherwise this code would need a BBox crop to the input's size.
      vw::rasterize(
        input_pre_type(
          m_child.prerasterize( BBox2i(m_xdelta*(*b).min().x(),
                                       m_ydelta*(*b).min().y(),
                                       m_xdelta*((*b).width () - 1) + 1,
                                       m_xdelta*((*b).height() - 1) + 1 ) ),
          m_xdelta, m_ydelta ),
        crop(buffer, *b - bbox.min()), *b );
    }

    return crop( buffer, -bbox.min().x(), -bbox.min().y(), cols(), rows() );
  }
  template <class DestT>
  inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
    vw::rasterize( prerasterize(bbox), dest, bbox );
  }
  /// \endcond
};


// *******************************************************************
// select_col()
// *******************************************************************

/// Return a single column from an image
/// \see vw::select_col
template <class ImageT>
class SelectColView : public ImageViewBase<SelectColView<ImageT> > {
  ImageT m_child;
  int32 m_col;
public:
  typedef typename ImageT::pixel_type pixel_type;
  typedef typename ImageT::result_type result_type;
  typedef typename ImageT::pixel_accessor pixel_accessor;

  SelectColView( ImageT const& image, int32 col ) : m_child(image), m_col(col) {}

  int32 cols() const { return 1; }
  inline int32 rows() const { return m_child.rows(); }
  inline int32 planes() const { return m_child.planes(); }

  inline pixel_accessor origin() const { return m_child.origin().advance(m_col,0,0); }
  inline result_type operator()( int32 /*i*/, int32 j, int32 p=0) const { return m_child(m_col,j,p); }

  SelectColView const& operator=( SelectColView const& view ) const {
    view.rasterize( *this, BBox2i(0,0,view.cols(),view.rows()) );
    return *this;
  }

  template <class ViewT>
  SelectColView const& operator=( ImageViewBase<ViewT> const& view ) const {
    view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
    return *this;
  }

  /// \cond INTERNAL
  typedef SelectColView<typename ImageT::prerasterize_type> prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const { return prerasterize_type( m_child.prerasterize(bbox), m_col ); }
  template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
  /// \endcond
};

// *******************************************************************
// select_row()
// *******************************************************************

/// Return a single row from an image
/// \see vw::select_row
template <class ImageT>
class SelectRowView : public ImageViewBase<SelectRowView<ImageT> > {
  ImageT m_child;
  int32 m_row;
public:
  typedef typename ImageT::pixel_type pixel_type;
  typedef typename ImageT::result_type result_type;
  typedef typename ImageT::pixel_accessor pixel_accessor;

  SelectRowView( ImageT const& image, int32 row ) : m_child(image), m_row(row) {}

  inline int32 cols() const { return m_child.cols(); }
  inline int32 rows() const { return 1; }
  inline int32 planes() const { return m_child.planes(); }

  inline pixel_accessor origin() const { return m_child.origin().advance(0,m_row,0); }

  inline result_type operator()( int32 i, int32 /*j*/, int32 p=0) const { return m_child(i,m_row,p); }

  SelectRowView const& operator=( SelectRowView const& view ) const {
    view.rasterize( *this, BBox2i(0,0,view.cols(),view.rows()) );
    return *this;
  }

  template <class ViewT>
  SelectRowView const& operator=( ImageViewBase<ViewT> const& view ) const {
    view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
    return *this;
  }

  typedef SelectRowView<typename ImageT::prerasterize_type> prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const { return prerasterize_type( m_child.prerasterize(bbox), m_row ); }
  template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
};

// *******************************************************************
// select_plane()
// *******************************************************************

/// Return a single plane from a multi-plane image
/// \see vw::select_plane
template <class ImageT>
class SelectPlaneView : public ImageViewBase<SelectPlaneView<ImageT> > {
  ImageT m_child;
  int32 m_plane;
public:
  typedef typename ImageT::pixel_type pixel_type;
  typedef typename ImageT::result_type result_type;
  typedef typename ImageT::pixel_accessor pixel_accessor;

  SelectPlaneView( ImageT const& image, int32 plane ) : m_child(image), m_plane(plane) {}

  inline int32 cols() const { return m_child.cols(); }
  inline int32 rows() const { return m_child.rows(); }
  inline int32 planes() const { return 1; }

  inline pixel_accessor origin() const { return m_child.origin().advance(0,0,m_plane); }
  inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { return m_child(i,j,m_plane+p); }

  template <class ViewT>
  SelectPlaneView const& operator=( ImageViewBase<ViewT> const& view ) const {
    view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
    return *this;
  }

  template <class ViewT>
  SelectPlaneView& operator=( ImageViewBase<ViewT> const& view ) {
    *const_cast<const SelectPlaneView*>(this) = view.impl();
    return *this;
  }

  /// \cond INTERNAL
  typedef SelectPlaneView<typename ImageT::prerasterize_type> prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const { return prerasterize_type( m_child.prerasterize(bbox), m_plane ); }
  template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
  /// \endcond
};


// *******************************************************************
// select_channel()
// *******************************************************************

/// A channel selecting functor, used by \ref select_channel().
template <class ImageT>
struct SelectChannelFunctor {
  int32 m_channel;
public:
  SelectChannelFunctor( int32 channel ) : m_channel(channel) {}

  // Computes an appropriate reference-to-channel type.
  typedef typename CompoundChannelType<typename ImageT::pixel_type>::type   channel_type;
  typedef typename CopyCVR<typename ImageT::result_type,channel_type>::type result_type;

  result_type operator()( typename ImageT::result_type pixel ) const {
    return compound_select_channel<result_type>(pixel,m_channel);
  }
};

/// A convenience function to select the alpha channel of an image.
template <class ImageT>
UnaryPerPixelView<ImageT,SelectChannelFunctor<ImageT> >
select_alpha_channel( ImageViewBase<ImageT>& image ) {
  // FIXME: This should be a static assertion
  if( ! PixelHasAlpha<typename ImageT::pixel_type>::value )
    vw_throw( ArgumentErr() << "Image has no alpha channel in call to select_alpha_channel()" );
  return select_channel( image, PixelNumChannels<typename ImageT::pixel_type>::value - 1 );
}

/// A convenience function to select the alpha channel of an image
/// (const overload).
template <class ImageT>
UnaryPerPixelView<ImageT,SelectChannelFunctor<const ImageT> >
select_alpha_channel( ImageViewBase<ImageT> const& image ) {
  // FIXME: This should be a static assertion
  if( ! PixelHasAlpha<typename ImageT::pixel_type>::value )
    vw_throw( ArgumentErr() << "Image has no alpha channel in call to select_alpha_channel()" );
  return select_channel( image, PixelNumChannels<typename ImageT::pixel_type>::value - 1 );
}


// *******************************************************************
// channels_to_planes()
// *******************************************************************

/// A channels-to-planes pixel accessor adaptor.
///
/// This is a special wrapper pixel accessor type, used by
/// \ref vw::ChannelsToPlanesView, that treats the channels in a
/// multi-channel image as planes.
template <class ChildT>
class ChannelsToPlanesAccessor
{
private:
  ChildT m_child;
  int32 m_channel;
public:
  typedef typename CompoundChannelType<typename ChildT::pixel_type>::type pixel_type;
  typedef typename CopyCVR<typename ChildT::result_type, pixel_type>::type result_type;
  typedef typename ChildT::offset_type offset_type;

  ChannelsToPlanesAccessor( ChildT const& acc ) : m_child(acc), m_channel(0) {}
  inline ChannelsToPlanesAccessor& next_col() { m_child.next_col(); return *this; }
  inline ChannelsToPlanesAccessor& prev_col() { m_child.prev_col(); return *this; }
  inline ChannelsToPlanesAccessor& next_row() { m_child.next_row(); return *this; }
  inline ChannelsToPlanesAccessor& prev_row() { m_child.prev_row(); return *this; }
  inline ChannelsToPlanesAccessor& next_plane() { ++m_channel; return *this; }
  inline ChannelsToPlanesAccessor& prev_plane() { --m_channel; return *this; }
  inline ChannelsToPlanesAccessor& advance( offset_type di, offset_type dj, ssize_t dp=0 ) { m_child.advance(di,dj); m_channel+=dp; return *this; }

  inline result_type operator*() const { return compound_select_channel<result_type>(*m_child,m_channel); }
};

/// A view that turns a one plane, multi-channel view into a mulit-plane, one channel view.
/// \see vw::channels_to_planes
template <class ImageT>
class ChannelsToPlanesView : public ImageViewBase<ChannelsToPlanesView<ImageT> > {
  ImageT m_child;
public:

  typedef typename CompoundChannelType<typename ImageT::pixel_type>::type pixel_type;
  typedef typename CopyCVR<typename ImageT::result_type, pixel_type>::type result_type;

  typedef typename boost::mpl::if_< IsCompound<typename ImageT::pixel_type>,
                                    ChannelsToPlanesAccessor<typename ImageT::pixel_accessor>,
                                    typename ImageT::pixel_accessor >::type pixel_accessor;

  ChannelsToPlanesView( ImageT const& image ) : m_child(image) {
    VW_ASSERT( m_child.planes()==1 , ArgumentErr() << "ChannelsToPlanesView: The image must be single plane.");
  }

  inline int32 cols() const { return m_child.cols(); }
  inline int32 rows() const { return m_child.rows(); }
  inline int32 planes() const { return m_child.channels(); }

  inline pixel_accessor origin() const { return pixel_accessor(m_child.origin()); }

  inline result_type operator()( int32 i, int32 j, int32 p=0 ) const {
    return compound_select_channel<result_type>(m_child(i,j),p);
  }

  template <class ViewT>
  ChannelsToPlanesView const& operator=( ImageViewBase<ViewT> const& view ) const {
    view.impl().rasterize( *this, BBox2i(0,0,view.impl().cols(),view.impl().rows()) );
    return *this;
  }

  template <class ViewT>
  ChannelsToPlanesView& operator=( ImageViewBase<ViewT> const& view ) {
    *const_cast<const ChannelsToPlanesView*>(this) = view.impl();
    return *this;
  }

  ImageT const& child() const {
    return m_child;
  }

  /// \cond INTERNAL
  typedef ChannelsToPlanesView<typename ImageT::prerasterize_type> prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const { return prerasterize_type( m_child.prerasterize(bbox) ); }
  template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
  /// \endcond
};

// *******************************************************************
// planes_to_channels()
// *******************************************************************

/// A view that turns a multi-plane, single-channel view into a
/// one-plane, multi-channel view.
/// \see vw::planes_to_channels
template <class PixelT, class ImageT>
class PlanesToChannelsView : public ImageViewBase<PlanesToChannelsView<PixelT,ImageT> > {
  ImageT m_child;
public:

  typedef PixelT pixel_type;
  typedef PixelT result_type;
  typedef ProceduralPixelAccessor<PlanesToChannelsView> pixel_accessor;

  PlanesToChannelsView( ImageT const& image ) : m_child(image) {
    VW_ASSERT( m_child.channels()==1 &&
               boost::numeric_cast<size_t>(m_child.planes())==CompoundNumChannels<PixelT>::value,
               ArgumentErr() << "PlanesToChannelsView: The image must be multi-plane, single-channel.");
  }

  inline int32 cols() const { return m_child.cols(); }
  inline int32 rows() const { return m_child.rows(); }
  inline int32 planes() const { return 1; }

  inline pixel_accessor origin() const { return pixel_accessor(*this); }

  inline result_type operator()( int32 i, int32 j, int32 /*p*/=0 ) const {
    result_type result;
    typedef typename CompoundChannelType<result_type>::type channel_type;
    for ( size_t c=0; c<CompoundNumChannels<PixelT>::value; ++c )
      compound_select_channel<channel_type&>(result,c) = m_child(i,j,c);
    return result;
  }

  ImageT const& child() const {
    return m_child;
  }

  /// \cond INTERNAL
  typedef PlanesToChannelsView<PixelT, typename ImageT::prerasterize_type> prerasterize_type;
  inline prerasterize_type prerasterize( BBox2i const& bbox ) const { return prerasterize_type( m_child.prerasterize(bbox) ); }
  template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
  /// \endcond
};


// *******************************************************************
// pixel_cast()
// *******************************************************************

/// A pixel casting functor, used by \ref pixel_cast().
template <class PixelT>
struct PixelCastFunctor : ReturnFixedType<PixelT> {
  template <class ArgT>
  inline PixelT operator()( ArgT pixel ) const {
    return pixel_cast<PixelT>(pixel);
  }
};



// *******************************************************************
// pixel_cast_rescale()
// *******************************************************************

/// A pixel casting functor, used by \ref pixel_cast_rescale().
template <class PixelT>
struct PixelCastRescaleFunctor : ReturnFixedType<PixelT> {
  template <class ArgT>
  inline PixelT operator()( ArgT pixel ) const {
    return pixel_cast_rescale<PixelT>(pixel);
  }
};


// *******************************************************************
// channel_cast()
// *******************************************************************

/// A pixel channel casting functor, used by \ref channel_cast().
template <class ChannelT>
struct PixelChannelCastFunctor : UnaryReturnBinaryTemplateBind2nd<CompoundChannelCast,ChannelT> {
  template <class ArgT>
  inline typename CompoundChannelCast<ArgT,ChannelT>::type operator()( ArgT const& pixel ) const {
    return channel_cast<ChannelT>(pixel);
  }
};



// *******************************************************************
// channel_cast_round_and_clamp()
// *******************************************************************

/// A pixel channel casting, rounding and clamping functor, used by
/// \ref channel_cast_round_and_clamp().
/// This uses the function with the same name from vw/Image/PixelTypeInfo.h.
template <class ChannelT>
struct PixelChannelCastRoundClampFunctor : UnaryReturnBinaryTemplateBind2nd<CompoundChannelCast,ChannelT> {
  template <class ArgT>
  inline typename CompoundChannelCast<ArgT,ChannelT>::type operator()( ArgT const& pixel ) const {
    return channel_cast_round_and_clamp<ChannelT>(pixel);
  }
};


// *******************************************************************
// channel_cast_rescale()
// *******************************************************************

/// A pixel channel casting and rescaling functor, used by
/// \ref channel_cast_rescale().
template <class ChannelT>
struct PixelChannelCastRescaleFunctor : UnaryReturnBinaryTemplateBind2nd<CompoundChannelCast,ChannelT> {
  template <class ArgT>
  inline typename CompoundChannelCast<ArgT,ChannelT>::type operator()( ArgT const& pixel ) const {
    return channel_cast_rescale<ChannelT>(pixel);
  }
};

// *******************************************************************
// weighted_rgb_to_gray()
// *******************************************************************
/* TODO: Why does this not compile here?
/// A weighted rgb-to-gray pixel conversion functor.
class WeightedRGBToGrayFunctor {
  double m_rw, m_gw, m_bw;
public:
  template <class ArgsT> struct result {};
  template <class FuncT, class ChannelT> struct result<FuncT(PixelRGB<ChannelT>)> { typedef PixelGray<ChannelT> type; };
  template <class FuncT, class ChannelT> struct result<FuncT(PixelRGBA<ChannelT>)> { typedef PixelGrayA<ChannelT> type; };
  WeightedRGBToGrayFunctor( double rw, double gw, double bw ) : m_rw(rw), m_gw(gw), m_bw(bw) {}
  template <class ChannelT> inline PixelGrayA<ChannelT> operator()( PixelRGBA<ChannelT> const& rgb ) const {
    return weighted_rgb_to_gray( rgb, m_rw, m_gw, m_bw );
  }
  template <class ChannelT> inline PixelGray<ChannelT> operator()( PixelRGB<ChannelT> const& rgb ) const {
    return weighted_rgb_to_gray( rgb, m_rw, m_gw, m_bw );
  }
};
*/
} // namespace vw

