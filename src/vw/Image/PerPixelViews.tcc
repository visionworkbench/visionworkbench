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
  // BinaryPerPixelView
  // *******************************************************************

  // Specialized Accessor
  template <class Image1IterT, class Image2IterT, class FuncT>
  class BinaryPerPixelAccessor {
    Image1IterT m_iter1;
    Image2IterT m_iter2;
    FuncT const& m_func;
  public:
    typedef typename boost::result_of<FuncT(typename Image1IterT::pixel_type,typename Image2IterT::pixel_type)>::type result_type;
    typedef typename boost::remove_cv<typename boost::remove_reference<result_type>::type>::type pixel_type;
    typedef typename boost::mpl::if_<boost::is_same<typename Image1IterT::offset_type, typename Image2IterT::offset_type>,
                                     typename Image1IterT::offset_type, int32>::type offset_type;

    BinaryPerPixelAccessor( Image1IterT const& iter1, Image2IterT const& iter2, FuncT const& func ) : m_iter1(iter1), m_iter2(iter2), m_func(func) {}
    inline BinaryPerPixelAccessor& next_col  () { m_iter1.next_col();   m_iter2.next_col();   return *this; }
    inline BinaryPerPixelAccessor& prev_col  () { m_iter1.prev_col();   m_iter2.prev_col();   return *this; }
    inline BinaryPerPixelAccessor& next_row  () { m_iter1.next_row();   m_iter2.next_row();   return *this; }
    inline BinaryPerPixelAccessor& prev_row  () { m_iter1.prev_row();   m_iter2.prev_row();   return *this; }
    inline BinaryPerPixelAccessor& next_plane() { m_iter1.next_plane(); m_iter2.next_plane(); return *this; }
    inline BinaryPerPixelAccessor& prev_plane() { m_iter1.prev_plane(); m_iter2.prev_plane(); return *this; }
    inline BinaryPerPixelAccessor& advance( offset_type di, offset_type dj, ssize_t dp=0 )
      { m_iter1.advance(di,dj,dp); m_iter2.advance(di,dj,dp); return *this; }
    inline result_type operator*() const { return m_func(*m_iter1,*m_iter2); }
  };

  // Image View Class Definition
  template <class Image1T, class Image2T, class FuncT>
  class BinaryPerPixelView : public ImageViewBase<BinaryPerPixelView<Image1T,Image2T,FuncT> >
  {
  private:
    Image1T m_image1;
    Image2T m_image2;
    FuncT   m_func;
  public:
    typedef typename boost::result_of<FuncT(typename Image1T::pixel_type, typename Image2T::pixel_type)>::type result_type;
    typedef typename boost::remove_cv<typename boost::remove_reference<result_type>::type>::type pixel_type;

    typedef BinaryPerPixelAccessor<typename Image1T::pixel_accessor,
                                   typename Image2T::pixel_accessor,
                                   FuncT> pixel_accessor;

    BinaryPerPixelView( Image1T const& image1, Image2T const& image2 )
      : m_image1(image1), m_image2(image2), m_func()
    {
      VW_ASSERT( m_image1.cols  ()==m_image2.cols  () &&
                 m_image1.rows  ()==m_image2.rows  () &&
                 m_image1.planes()==m_image2.planes(),
                 ArgumentErr() << "BinaryPerPixelView: Images must have same dimensions in binary image operation." );
    }

    BinaryPerPixelView( Image1T const& image1, Image2T const& image2, FuncT const& func )
      : m_image1(image1), m_image2(image2), m_func(func)
    {
      VW_ASSERT( m_image1.cols  ()==m_image2.cols  () &&
                 m_image1.rows  ()==m_image2.rows  () &&
                 m_image1.planes()==m_image2.planes(),
                 ArgumentErr() << "BinaryPerPixelView: Images must have same dimensions in binary image operation." );
    }

    inline int32 cols  () const { return m_image1.cols  (); }
    inline int32 rows  () const { return m_image1.rows  (); }
    inline int32 planes() const { return m_image1.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor(m_image1.origin(),m_image2.origin(),m_func); }
    inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { return m_func(m_image1(i,j,p),m_image2(i,j,p)); }

    /// \cond INTERNAL
    typedef BinaryPerPixelView<typename Image1T::prerasterize_type, typename Image2T::prerasterize_type, FuncT> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const { return prerasterize_type( m_image1.prerasterize(bbox), m_image2.prerasterize(bbox), m_func ); }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };

  // *******************************************************************
  // TrinaryPerPixelView
  // *******************************************************************

  // Specialized Accessor
  template <class Image1IterT, class Image2IterT, class Image3IterT, class FuncT>
  class TrinaryPerPixelAccessor {
    Image1IterT m_iter1;
    Image2IterT m_iter2;
    Image3IterT m_iter3;
    FuncT const& m_func;
  public:
    typedef typename boost::result_of<FuncT(typename Image1IterT::pixel_type,typename Image2IterT::pixel_type,typename Image3IterT::pixel_type)>::type result_type;
    typedef typename boost::remove_cv<typename boost::remove_reference<result_type>::type>::type pixel_type;
    typedef typename boost::mpl::if_<boost::mpl::and_<boost::is_same<typename Image1IterT::offset_type, typename Image2IterT::offset_type>,
                                                      boost::is_same<typename Image1IterT::offset_type, typename Image3IterT::offset_type> >,
                                     typename Image1IterT::offset_type, int32>::type offset_type;

    TrinaryPerPixelAccessor( Image1IterT const& iter1, Image2IterT const& iter2, Image2IterT const& iter3, FuncT const& func ) : m_iter1(iter1), m_iter2(iter2), m_iter3(iter3), m_func(func) {}
    inline TrinaryPerPixelAccessor& next_col  () { m_iter1.next_col();   m_iter2.next_col();   m_iter3.next_col();   return *this; }
    inline TrinaryPerPixelAccessor& prev_col  () { m_iter1.prev_col();   m_iter2.prev_col();   m_iter3.prev_col();   return *this; }
    inline TrinaryPerPixelAccessor& next_row  () { m_iter1.next_row();   m_iter2.next_row();   m_iter3.next_row();   return *this; }
    inline TrinaryPerPixelAccessor& prev_row  () { m_iter1.prev_row();   m_iter2.prev_row();   m_iter3.prev_row();   return *this; }
    inline TrinaryPerPixelAccessor& next_plane() { m_iter1.next_plane(); m_iter2.next_plane(); m_iter3.next_plane(); return *this; }
    inline TrinaryPerPixelAccessor& prev_plane() { m_iter1.prev_plane(); m_iter2.prev_plane(); m_iter3.prev_plane(); return *this; }
    inline TrinaryPerPixelAccessor& advance( offset_type di, offset_type dj, ssize_t dp=0 )
      { m_iter1.advance(di,dj,dp); m_iter2.advance(di,dj,dp); m_iter3.advance(di,dj,dp); return *this; }
    inline result_type operator*() const { return m_func(*m_iter1,*m_iter2,*m_iter3); }
  };

  // Image View Class Definition
  template <class Image1T, class Image2T, class Image3T, class FuncT>
  class TrinaryPerPixelView : public ImageViewBase<TrinaryPerPixelView<Image1T,Image2T,Image3T,FuncT> >
  {
  private:
    Image1T m_image1;
    Image2T m_image2;
    Image3T m_image3;
    FuncT m_func;
  public:
    typedef typename boost::result_of<FuncT(typename Image1T::pixel_type, typename Image2T::pixel_type, typename Image3T::pixel_type)>::type result_type;
    typedef typename boost::remove_cv<typename boost::remove_reference<result_type>::type>::type pixel_type;

    typedef TrinaryPerPixelAccessor<typename Image1T::pixel_accessor,
                                    typename Image2T::pixel_accessor,
                                    typename Image3T::pixel_accessor,
                                    FuncT> pixel_accessor;

    TrinaryPerPixelView( Image1T const& image1, Image2T const& image2, Image3T const& image3 )
      : m_image1(image1), m_image2(image2), m_image3(image3), m_func()
    {
      VW_ASSERT( m_image1.cols  ()==m_image2.cols  () &&
                 m_image1.rows  ()==m_image2.rows  () &&
                 m_image1.planes()==m_image2.planes() &&
                 m_image1.cols  ()==m_image3.cols  () &&
                 m_image1.rows  ()==m_image3.rows  () &&
                 m_image1.planes()==m_image3.planes(),
                 ArgumentErr() << "TrinaryPerPixelView: Images must have same dimensions in trinary image operation." );
    }

    TrinaryPerPixelView( Image1T const& image1, Image2T const& image2, Image3T const& image3, FuncT const& func )
      : m_image1(image1), m_image2(image2), m_image3(image3), m_func(func)
    {
      VW_ASSERT( m_image1.cols  ()==m_image2.cols  () &&
                 m_image1.rows  ()==m_image2.rows  () &&
                 m_image1.planes()==m_image2.planes() &&
                 m_image1.cols  ()==m_image3.cols  () &&
                 m_image1.rows  ()==m_image3.rows  () &&
                 m_image1.planes()==m_image3.planes(),
                 ArgumentErr() << "TrinaryPerPixelView: Images must have same dimensions in trinary image operation." );
    }

    inline int32 cols  () const { return m_image1.cols  (); }
    inline int32 rows  () const { return m_image1.rows  (); }
    inline int32 planes() const { return m_image1.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor(m_image1.origin(),m_image2.origin(),m_image3.origin(),m_func); }
    inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { return m_func(m_image1(i,j,p),m_image2(i,j,p),m_image3(i,j,p)); }

    /// \cond INTERNAL
    typedef TrinaryPerPixelView<typename Image1T::prerasterize_type, typename Image2T::prerasterize_type, typename Image3T::prerasterize_type, FuncT> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const { return prerasterize_type( m_image1.prerasterize(bbox), m_image2.prerasterize(bbox), m_image3.prerasterize(bbox), m_func ); }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };


// *******************************************************************
  // QuaternaryPerPixelView
  // *******************************************************************

  // Specialized Accessor
  template <class Image1IterT, class Image2IterT, class Image3IterT, class Image4IterT, class FuncT>
  class QuaternaryPerPixelAccessor {
    Image1IterT m_iter1;
    Image2IterT m_iter2;
    Image3IterT m_iter3;
    Image4IterT m_iter4;
    FuncT const& m_func;
  public:
    typedef typename boost::result_of<FuncT(typename Image1IterT::pixel_type,typename Image2IterT::pixel_type,
                                            typename Image3IterT::pixel_type,typename Image4IterT::pixel_type)>::type result_type;
    typedef typename boost::remove_cv<typename boost::remove_reference<result_type>::type>::type pixel_type;
    typedef typename boost::mpl::if_<boost::mpl::and_<boost::is_same<typename Image1IterT::offset_type, typename Image2IterT::offset_type>,
                                                      boost::is_same<typename Image1IterT::offset_type, typename Image3IterT::offset_type>,
                                                      boost::is_same<typename Image1IterT::offset_type, typename Image4IterT::offset_type> >,
                                     typename Image1IterT::offset_type, int32>::type offset_type;

    QuaternaryPerPixelAccessor( Image1IterT const& iter1, Image2IterT const& iter2, 
                                Image3IterT const& iter3, Image4IterT const& iter4, FuncT const& func )
                                 : m_iter1(iter1), m_iter2(iter2), m_iter3(iter3), m_iter4(iter4), m_func(func) {}
    inline QuaternaryPerPixelAccessor& next_col  () { m_iter1.next_col();   m_iter2.next_col();   m_iter3.next_col();   m_iter4.next_col();   return *this; }
    inline QuaternaryPerPixelAccessor& prev_col  () { m_iter1.prev_col();   m_iter2.prev_col();   m_iter3.prev_col();   m_iter4.prev_col();   return *this; }
    inline QuaternaryPerPixelAccessor& next_row  () { m_iter1.next_row();   m_iter2.next_row();   m_iter3.next_row();   m_iter4.next_row();   return *this; }
    inline QuaternaryPerPixelAccessor& prev_row  () { m_iter1.prev_row();   m_iter2.prev_row();   m_iter3.prev_row();   m_iter4.prev_row();   return *this; }
    inline QuaternaryPerPixelAccessor& next_plane() { m_iter1.next_plane(); m_iter2.next_plane(); m_iter3.next_plane(); m_iter4.next_plane(); return *this; }
    inline QuaternaryPerPixelAccessor& prev_plane() { m_iter1.prev_plane(); m_iter2.prev_plane(); m_iter3.prev_plane(); m_iter4.prev_plane(); return *this; }
    inline QuaternaryPerPixelAccessor& advance( offset_type di, offset_type dj, ssize_t dp=0 )
      { m_iter1.advance(di,dj,dp); m_iter2.advance(di,dj,dp); m_iter3.advance(di,dj,dp); m_iter4.advance(di,dj,dp); return *this; }
    inline result_type operator*() const { return m_func(*m_iter1,*m_iter2,*m_iter3,*m_iter4); }
  };

  // Image View Class Definition
  template <class Image1T, class Image2T, class Image3T, class Image4T, class FuncT>
  class QuaternaryPerPixelView : public ImageViewBase<QuaternaryPerPixelView<Image1T,Image2T,Image3T,Image4T,FuncT> >
  {
  private:
    Image1T m_image1;
    Image2T m_image2;
    Image3T m_image3;
    Image4T m_image4;
    FuncT m_func;
  public:
    typedef typename boost::result_of<FuncT(typename Image1T::pixel_type, typename Image2T::pixel_type, 
                                            typename Image3T::pixel_type, typename Image4T::pixel_type)>::type result_type;
    typedef typename boost::remove_cv<typename boost::remove_reference<result_type>::type>::type pixel_type;

    typedef QuaternaryPerPixelAccessor<typename Image1T::pixel_accessor,
                                       typename Image2T::pixel_accessor,
                                       typename Image3T::pixel_accessor,
                                       typename Image4T::pixel_accessor,
                                       FuncT> pixel_accessor;

    QuaternaryPerPixelView( Image1T const& image1, Image2T const& image2, Image3T const& image3, Image4T const& image4 )
      : m_image1(image1), m_image2(image2), m_image3(image3), m_image4(image4), m_func()
    {
      VW_ASSERT( m_image1.cols  ()==m_image2.cols  () &&
                 m_image1.rows  ()==m_image2.rows  () &&
                 m_image1.planes()==m_image2.planes() &&
                 m_image1.cols  ()==m_image3.cols  () &&
                 m_image1.rows  ()==m_image3.rows  () &&
                 m_image1.planes()==m_image3.planes() &&
                 m_image1.cols  ()==m_image4.cols  () &&
                 m_image1.rows  ()==m_image4.rows  () &&
                 m_image1.planes()==m_image4.planes(),
                 ArgumentErr() << "QuaternaryPerPixelView: Images must have same dimensions in quaternary image operation." );
    }

    QuaternaryPerPixelView( Image1T const& image1, Image2T const& image2, Image3T const& image3, Image4T const& image4, FuncT const& func )
      : m_image1(image1), m_image2(image2), m_image3(image3), m_image4(image4), m_func(func)
    {
      VW_ASSERT( m_image1.cols  ()==m_image2.cols  () &&
                 m_image1.rows  ()==m_image2.rows  () &&
                 m_image1.planes()==m_image2.planes() &&
                 m_image1.cols  ()==m_image3.cols  () &&
                 m_image1.rows  ()==m_image3.rows  () &&
                 m_image1.planes()==m_image3.planes() &&
                 m_image1.cols  ()==m_image4.cols  () &&
                 m_image1.rows  ()==m_image4.rows  () &&
                 m_image1.planes()==m_image4.planes(),
                 ArgumentErr() << "QuaternaryPerPixelView: Images must have same dimensions in quaternary image operation." );
    }

    inline int32 cols  () const { return m_image1.cols  (); }
    inline int32 rows  () const { return m_image1.rows  (); }
    inline int32 planes() const { return m_image1.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor(m_image1.origin(),m_image2.origin(),m_image3.origin(),m_image4.origin(),m_func); }
    inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { return m_func(m_image1(i,j,p),m_image2(i,j,p),m_image3(i,j,p),m_image4(i,j,p)); }

    /// \cond INTERNAL
    typedef QuaternaryPerPixelView<typename Image1T::prerasterize_type, typename Image2T::prerasterize_type, 
                                   typename Image3T::prerasterize_type, typename Image4T::prerasterize_type, FuncT> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const { return prerasterize_type( m_image1.prerasterize(bbox), m_image2.prerasterize(bbox),
                                                                                                  m_image3.prerasterize(bbox), m_image4.prerasterize(bbox), m_func ); }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };





} // End namespace vw

