// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file UtilityViews.h
///
/// Utility objects conforming to the view concept.
///
#ifndef __VW_IMAGE_UTILITYVIEWS_H__
#define __VW_IMAGE_UTILITYVIEWS_H__

#include <boost/smart_ptr.hpp>
#include <boost/type_traits.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>

#include <vw/Image/ImageViewBase.h>
#include <vw/Image/PixelAccessors.h>
#include <vw/Image/PerPixelViews.h>

namespace vw {

  // *******************************************************************
  // constant_view()
  // *******************************************************************

  template <class PixelT>
  struct ConstantIndexFunctor {
    typedef PixelT result_type;
    result_type m_value;
    ConstantIndexFunctor(result_type const& value) : m_value(value) {}
    result_type operator()(double i, double j, int32 p) const {
      return m_value;
    }
  };

  template <class PixelT>
  inline PerPixelIndexView<ConstantIndexFunctor<PixelT> >
  constant_view( PixelT const& value, int32 cols, int32 rows, int32 planes = 1 ) {
    typedef PerPixelIndexView<ConstantIndexFunctor<PixelT> > result_type;
    return result_type( ConstantIndexFunctor<PixelT>(value), cols, rows, planes );
  }

  template <class PixelT, class ImageT>
  inline PerPixelIndexView<ConstantIndexFunctor<PixelT> >
  constant_view( PixelT value, ImageViewBase<ImageT> const& image ) { 
    return constant_view(value,
                         image.impl().cols(),
                         image.impl().rows(),
                         image.impl().planes());
  }

  // *******************************************************************
  // pixel_index_view()
  // *******************************************************************

  struct VectorIndexFunctor {
    typedef Vector2 result_type;
    result_type operator()(double i, double j, int32 p) const {
      return Vector2(i, j);
    }  
  };

  inline PerPixelIndexView<VectorIndexFunctor>
  pixel_index_view( int32 cols, int32 rows, int32 planes = 1 ) {
    typedef PerPixelIndexView<VectorIndexFunctor> result_type;
    return result_type( VectorIndexFunctor(), cols, rows, planes );
  }

  template <class ImageT>
  inline PerPixelIndexView<VectorIndexFunctor>
  pixel_index_view( ImageViewBase<ImageT> const& image ) {
    return pixel_index_view(image.impl().cols(),
                            image.impl().rows(),
                            image.impl().planes());
  }

  // *******************************************************************
  // pixel_index3_view()
  // *******************************************************************

  struct Vector3IndexFunctor {
    typedef Vector3 result_type;
    result_type operator()(double i, double j, int32 p) const {
      return Vector3(i, j, p);
    }
  };

  inline PerPixelIndexView<Vector3IndexFunctor>
  pixel_index3_view( int32 cols, int32 rows, int32 planes = 1 ) {
    typedef PerPixelIndexView<Vector3IndexFunctor> result_type;
    return result_type( Vector3IndexFunctor(), cols, rows, planes );
  }

  template <class ImageT>
  inline PerPixelIndexView<Vector3IndexFunctor>
  pixel_index3_view( ImageViewBase<ImageT> const& image ) {
    return pixel_index3_view(image.impl().cols(),
                             image.impl().rows(),
                             image.impl().planes());
  }

  // *******************************************************************
  // uniform_noise_view(), gaussian_noise_view()
  // *******************************************************************

  template <class FuncT>
  struct IndexIndependentFunctor {
    typedef typename FuncT::result_type result_type;
    mutable FuncT m_func;
    IndexIndependentFunctor(FuncT func) : m_func(func) {}
    result_type operator()(double i, double j, int32 p) const {
      return m_func();
    }
  };

  template <class GenT>
  inline PerPixelIndexView<IndexIndependentFunctor<boost::variate_generator<GenT&, boost::uniform_01<> > > >
  uniform_noise_view( GenT& gen, int32 cols, int32 rows, int32 planes = 1) {
    typedef boost::variate_generator<GenT&, boost::uniform_01<> > vargen_type;
    typedef PerPixelIndexView<IndexIndependentFunctor<vargen_type> > return_type;
    boost::uniform_01<> dist;
    vargen_type vargen(gen, dist);
    return return_type( IndexIndependentFunctor<vargen_type>(vargen), 
                        cols, rows, planes );
  }

  template <class GenT, class ImageT>
  inline PerPixelIndexView<IndexIndependentFunctor<boost::variate_generator<GenT&, boost::uniform_01<> > > >
  uniform_noise_view( GenT& gen, ImageViewBase<ImageT> const& image ) {
    return uniform_noise_view(gen,
                              image.impl().cols(),
                              image.impl().rows(),
                              image.impl().planes());
  }

  template <class GenT>
  inline PerPixelIndexView<IndexIndependentFunctor<boost::variate_generator<GenT&, boost::normal_distribution<> > > >
  gaussian_noise_view( GenT& gen, double mean, double sigma,
                       int32 cols, int32 rows, int32 planes = 1 ) {
    typedef boost::variate_generator<GenT&, boost::normal_distribution<> > vargen_type;
    typedef PerPixelIndexView<IndexIndependentFunctor<vargen_type> > return_type;
    boost::normal_distribution<> dist(mean, sigma);
    vargen_type vargen(gen, dist);
    return return_type( IndexIndependentFunctor<vargen_type>(vargen),
                        cols, rows, planes );
  }

  template <class GenT, class ImageT>
  inline PerPixelIndexView<IndexIndependentFunctor<boost::variate_generator<GenT&, boost::normal_distribution<> > > >
  gaussian_noise_view( GenT& gen, double mean, double sigma,
                       ImageViewBase<ImageT> const& image ) {
    return gaussian_noise_view(gen, mean, sigma, 
                               image.impl().cols(), 
                               image.impl().rows(), 
                               image.impl().planes());
  }

} // namespace vw

#endif // __VW_IMAGE_IMAGEVIEW_H__
