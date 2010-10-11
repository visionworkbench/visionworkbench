// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef Transforms_H
#define Transforms_H

#include <cmath>

#include <vw/Image.h>
#include <vw/Math.h>

#include <vw/GPU/Setup.h>
#include <vw/GPU/GPUImage.h>
#include <vw/GPU/GPUProgram.h>
#include <vw/GPU/Interpolation.h>
#include <vw/GPU/EdgeExtension.h>


namespace vw { namespace GPU {

  // fixed_homography_transform

  template <class InterpT, class EdgeT>
  inline GPUImageBase fixed_homography_transform(GPUImageBase const& image, Matrix<float>& homography, InterpT, EdgeT) {
    GPUImageBase temp = image;
    temp.apply_homography(homography, InterpT(), EdgeT());
    return temp;
  }

  template <class InterpT>
  inline GPUImageBase fixed_homography_transform(GPUImageBase const& image, Matrix<float>& homography, InterpT) {
    GPUImageBase temp = image;
    temp.apply_homography(homography, InterpT(), DefaultEdgeExtension());
    return temp;
  }

  inline GPUImageBase fixed_homography_transform(GPUImageBase const& image, Matrix<float>& homography) {
    GPUImageBase temp = image;
    temp.apply_homography(homography, DefaultInterpolation(), DefaultEdgeExtension());
    return temp;
  }

  template <class PixelT, class InterpT, class EdgeT>
  inline GPUImage<PixelT> fixed_homography_transform(GPUImage<PixelT> const& image, Matrix<float>& homography, InterpT, EdgeT) {
    GPUImage<PixelT> temp = image;
    image.apply_homography(homography, InterpT(), EdgeT());
    return temp;
  }

  template <class PixelT, class InterpT>
  inline GPUImage<PixelT> fixed_homography_transform(GPUImage<PixelT> const& image, Matrix<float>& homography, InterpT) {
    GPUImage<PixelT> temp = image;
    image.apply_homography(homography, DefaultInterpolation(), DefaultEdgeExtension());
    return temp;
  }

  template <class PixelT>
  inline GPUImage<PixelT> fixed_homography_transform(GPUImage<PixelT> const& image, Matrix<float>& homography) {
    GPUImage<PixelT> temp = image;
    image.apply_homography(homography, DefaultInterpolation(), DefaultEdgeExtension());
    return temp;
  }



  // fixed_rotate

  template <class InterpT, class EdgeT>
  inline GPUImageBase fixed_rotate(GPUImageBase const& image, float theta, InterpT, EdgeT) {
    int width = image.width();
    int height = image.height();
    Matrix<float> h_translate(3, 3);
    h_translate.set_identity();
    h_translate(0,2) = width/2.0;
    h_translate(1,2) = height/2.0;
    Matrix<float> h_rotate(3, 3);
    h_rotate.set_identity();
    h_rotate(0,0) = ::cos(theta);
    h_rotate(0,1) = ::sin(theta);
    h_rotate(1,0) = -::sin(theta);
    h_rotate(1,1) = ::cos(theta);
    Matrix<float> h_translate_back(3, 3);
    h_translate_back.set_identity();
    h_translate_back(0,2) = -width/2.0;
    h_translate_back(1,2) = -height/2.0;
    Matrix<float> h = h_translate * h_rotate * h_translate_back;
    GPUImageBase temp = image;
    temp.apply_homography(h, InterpT(), EdgeT());
    return temp;
  }

  template <class InterpT>
  inline GPUImageBase fixed_rotate(GPUImageBase const& image, float theta, InterpT) {
    return fixed_rotate(image, theta, InterpT(), DefaultEdgeExtension());
  }

  inline GPUImageBase fixed_rotate(GPUImageBase const& image, float theta) {
    return fixed_rotate(image, theta, DefaultInterpolation(), DefaultEdgeExtension());
  }

 template <class PixelT, class InterpT, class EdgeT>
  inline GPUImage<PixelT> fixed_rotate(GPUImage<PixelT> const& image, float theta, InterpT, EdgeT) {
    return fixed_rotate((GPUImageBase&) image, theta, InterpT(), EdgeT());
  }

  template <class PixelT, class InterpT>
  inline GPUImage<PixelT> fixed_rotate(GPUImage<PixelT> const& image, float theta, InterpT) {
    return fixed_rotate((GPUImageBase&) image, theta, InterpT(), DefaultEdgeExtension());
  }

  template <class PixelT>
  inline GPUImage<PixelT> fixed_rotate(GPUImage<PixelT> const& image, float theta) {
    return fixed_rotate((GPUImageBase&) image, theta, DefaultInterpolation(), DefaultEdgeExtension());
  }


  // resample

 template <class InterpT, class EdgeT>
  inline GPUImageBase resample(GPUImageBase const& image, float x_scale_factor, float y_scale_factor, int new_width, int new_height, InterpT, EdgeT) {
    Matrix<float> h_resample(3, 3);
    h_resample.set_identity();
    h_resample(0,0) = 1 / x_scale_factor;
    h_resample(1,1) = 1 / y_scale_factor;
    GPUImageBase temp = image;
    if(new_width == 0) new_width = image.width();
    if(new_height == 0) new_height = image.height();
    temp.apply_homography(h_resample, InterpT(), EdgeT(), new_width, new_height);
    return temp;
  }

  template <class InterpT>
  inline GPUImageBase resample(GPUImageBase const& image, float x_scale_factor, float y_scale_factor, int new_width, int new_height, InterpT) {
    return resample(image, x_scale_factor, y_scale_factor, new_width, new_height, InterpT(), DefaultEdgeExtension());
  }

  inline GPUImageBase resample(GPUImageBase const& image, float x_scale_factor, float y_scale_factor, int new_width, int new_height) {
    return resample(image, x_scale_factor, y_scale_factor, new_width, new_height, DefaultInterpolation(), DefaultEdgeExtension());
  }

 template <class PixelT, class InterpT, class EdgeT>
  inline GPUImage<PixelT> resample(GPUImage<PixelT> const& image, float x_scale_factor, float y_scale_factor, int new_width, int new_height, InterpT, EdgeT) {
    return resample((GPUImageBase&) image, x_scale_factor, y_scale_factor, new_width, new_height, InterpT(), EdgeT());
  }

  template <class PixelT, class InterpT>
  inline GPUImage<PixelT> resample(GPUImage<PixelT> const& image, float x_scale_factor, float y_scale_factor, int new_width, int new_height, InterpT) {
    return resample((GPUImageBase&) image, x_scale_factor, y_scale_factor, new_width, new_height, InterpT(), DefaultEdgeExtension());
  }

  template <class PixelT>
  inline GPUImage<PixelT> resample(GPUImage<PixelT> const& image, float x_scale_factor, float y_scale_factor, int new_width, int new_height) {
    return resample((GPUImageBase&) image, x_scale_factor, y_scale_factor, new_width, new_height, DefaultInterpolation(), DefaultEdgeExtension());
  }

  /*
 template <class InterpT, class EdgeT>
  GPUImageBase resample(GPUImageBase const& image, float x_scale_factor, float y_scale_factor, InterpT, EdgeT) {
    int new_width = (int) floor(x_scale_factor * image.cols());
    int new_height = (int) floor(y_scale_factor * image.rows());
    return resample(image, x_scale_factor, y_scale_factor, new_width, new_height, InterpT(), DefaultEdgeExtension());
  }

  template <class InterpT>
  GPUImageBase resample(GPUImageBase const& image, float x_scale_factor, float y_scale_factor, InterpT) {
    int new_width = (int) floor(x_scale_factor * image.cols());
    int new_height = (int) floor(y_scale_factor * image.rows());
    return resample(image, x_scale_factor, y_scale_factor, new_width, new_height, InterpT(), DefaultEdgeExtension());
  }

  GPUImageBase resample(GPUImageBase const& image, float x_scale_factor, float y_scale_factor) {
    int new_width = (int) floor(x_scale_factor * image.cols());
    int new_height = (int) floor(y_scale_factor * image.rows());
    return resample(image, x_scale_factor, y_scale_factor, new_width, new_height, interpolation::Default(), DefaultEdgeExtension());
  }

 template <class PixelT, class InterpT, class EdgeT>
  inline GPUImage<PixelT> resample(GPUImage<PixelT> const& image, float x_scale_factor, float y_scale_factor, InterpT, EdgeT) {
    int new_width = (int) floor(x_scale_factor * image.cols());
    int new_height = (int) floor(x_scale_factor * image.rows());
    return resample((GPUImageBase&) image, x_scale_factor, y_scale_factor, new_width, new_height, InterpT(), EdgeT());
  }

  template <class PixelT, class InterpT>
  inline GPUImage<PixelT> resample(GPUImage<PixelT> const& image, float x_scale_factor, float y_scale_factor, InterpT) {
    int new_width = (int) floor(x_scale_factor * image.cols());
    int new_height = (int) floor(y_scale_factor * image.rows());
    return resample((GPUImageBase&) image, x_scale_factor, y_scale_factor, new_width, new_height, InterpT(), DefaultEdgeExtension());
  }

  template <class PixelT>
  inline GPUImage<PixelT> resample(GPUImage<PixelT> const& image, float x_scale_factor, float y_scale_factor) {
    int new_width = (int) floor(x_scale_factor * image.cols());
    int new_height = (int) floor(y_scale_factor * image.rows());
    return resample((GPUImageBase&) image, x_scale_factor, y_scale_factor, new_width, new_height, interpolation::Default(), DefaultEdgeExtension());
  }
  */

  // translate

 template <class InterpT, class EdgeT>
  inline GPUImageBase translate(GPUImageBase const& image, float x_translation, float y_translation, InterpT, EdgeT) {
    Matrix<float> h_translate(3, 3);
    h_translate.set_identity();
    h_translate(0,2) = x_translation;
    h_translate(1,2) = y_translation;
    GPUImageBase temp = image;
    temp.apply_homography(h_translate, InterpT(), EdgeT());
    return temp;
  }

  template <class InterpT>
  inline GPUImageBase translate(GPUImageBase const& image, float x_translation, float y_translation, InterpT) {
    return translate(image, x_translation, y_translation, InterpT(), DefaultEdgeExtension());
  }

  inline GPUImageBase translate(GPUImageBase const& image, float x_translation, float y_translation) {
    return translate(image, x_translation, y_translation, DefaultInterpolation(), DefaultEdgeExtension());
  }

 template <class PixelT, class InterpT, class EdgeT>
  inline GPUImage<PixelT> translate(GPUImage<PixelT> const& image, float x_translation, float y_translation, InterpT, EdgeT) {
    return translate((GPUImageBase&) image, x_translation, y_translation, InterpT(), EdgeT());
  }

  template <class PixelT, class InterpT>
  inline GPUImage<PixelT> translate(GPUImage<PixelT> const& image, float x_translation, float y_translation, InterpT) {
    return translate((GPUImageBase&) image, x_translation, y_translation, InterpT(), DefaultEdgeExtension());
  }

  template <class PixelT>
  inline GPUImage<PixelT> translate(GPUImage<PixelT> const& image, float x_translation, float y_translation) {
    return translate((GPUImageBase&) image, x_translation, y_translation, DefaultInterpolation(), DefaultEdgeExtension());
  }



  // free_homography_transform

  template <class InterpT, class EdgeT>
  inline GPUImageBase free_homography_transform(GPUImageBase const& image, Matrix<float>& homography, InterpT, EdgeT) {
    HomographyTransform h_functor(homography);
    Vector2 pt(0,0);
    Vector2 pt_min, pt_max, pt_new;
    pt_min = h_functor.reverse(pt);  pt_max = h_functor.reverse(pt);
    for (pt[0] = 0; pt[0] < image.width(); (pt[0])++) {
      for (pt[1] = 0; pt[1] < image.height(); (pt[1])++) {
        pt_new = h_functor.reverse(pt);
        pt_min(0) = (pt_new(0) < pt_min(0)) ? pt_new(0) : pt_min(0);
        pt_max(0) = (pt_new(0) > pt_max(0)) ? pt_new(0) : pt_max(0);
        pt_min(1) = (pt_new(1) < pt_min(1)) ? pt_new(1) : pt_min(1);
        pt_max(1) = (pt_new(1) > pt_max(1)) ? pt_new(1) : pt_max(1);
      }
    }

    Matrix<float> h_translate(3, 3);
    h_translate.set_identity();
    h_translate(0,2) = -pt_min(0);
    h_translate(1,2) = -pt_min(1);
    Matrix<float> h_out(3, 3);
    h_out = h_translate * homography;

    GPUImageBase temp = image;
    temp.apply_homography(h_out, InterpT(), EdgeT(), (int) (pt_max(0) - pt_min(0)), (int) (pt_max(1) - pt_min(1)));
    return temp;
  }

  template <class InterpT>
  inline GPUImageBase free_homography_transform(GPUImageBase const& image, Matrix<float>& homography, InterpT) {
    return free_homography_transform(image, homography, InterpT(), DefaultEdgeExtension());
  }

  inline GPUImageBase free_homography_transform(GPUImageBase const& image, Matrix<float>& homography) {
    return free_homography_transform(image, homography, DefaultInterpolation(), DefaultEdgeExtension());
  }

  template <class PixelT, class InterpT, class EdgeT>
  inline GPUImage<PixelT> free_homography_transform(GPUImage<PixelT> const& image, Matrix<float>& homography, InterpT, EdgeT) {
    return free_homography_transform((GPUImageBase&) image, homography, InterpT(), EdgeT());
  }

  template <class PixelT, class InterpT>
  inline GPUImage<PixelT> free_homography_transform(GPUImage<PixelT> const& image, Matrix<float>& homography, InterpT) {
    return free_homography_transform((GPUImageBase&) image, homography, InterpT(), DefaultEdgeExtension());
  }

  template <class PixelT>
  inline GPUImage<PixelT> free_homography_transform(GPUImage<PixelT> const& image, Matrix<float>& homography) {
    return free_homography_transform((GPUImageBase&) image, homography, DefaultInterpolation(), DefaultEdgeExtension());
  }

  // free_rotate

  template <class InterpT, class EdgeT>
  inline GPUImageBase free_rotate(GPUImageBase const& image, float theta, InterpT, EdgeT) {
    int width = image.width();
    int height = image.height();
    Matrix<float> h_translate(3, 3);
    h_translate.set_identity();
    h_translate(0,2) = width/2.0;
    h_translate(1,2) = height/2.0;
    Matrix<float> h_rotate(3, 3);
    h_rotate.set_identity();
    h_rotate(0,0) = ::cos(theta);
    h_rotate(0,1) = ::sin(theta);
    h_rotate(1,0) = -::sin(theta);
    h_rotate(1,1) = ::cos(theta);
    Matrix<float> h_translate_back(3, 3);
    h_translate_back.set_identity();
    h_translate_back(0,2) = -width/2.0;
    h_translate_back(1,2) = -height/2.0;
    Matrix<float> h = h_translate * h_rotate * h_translate_back;
    GPUImageBase temp = image;
    return free_homography_transform(image, h, InterpT(), EdgeT());
  }


  template <class InterpT>
  inline GPUImageBase free_rotate(GPUImageBase const& image, float theta, InterpT) {
    return free_rotate(image, theta, InterpT(), DefaultEdgeExtension());
  }

  inline GPUImageBase free_rotate(GPUImageBase const& image, float theta) {
    return free_rotate(image, theta, DefaultInterpolation(), DefaultEdgeExtension());
  }

 template <class PixelT, class InterpT, class EdgeT>
  inline GPUImage<PixelT> free_rotate(GPUImage<PixelT> const& image, float theta, InterpT, EdgeT) {
    return free_rotate((GPUImageBase&) image, theta, InterpT(), EdgeT());
  }

  template <class PixelT, class InterpT>
  inline GPUImage<PixelT> free_rotate(GPUImage<PixelT> const& image, float theta, InterpT) {
    return free_rotate((GPUImageBase&) image, theta, InterpT(), DefaultEdgeExtension());
  }

  template <class PixelT>
  inline GPUImage<PixelT> free_rotate(GPUImage<PixelT> const& image, float theta) {
    return free_rotate((GPUImageBase&) image, theta, DefaultInterpolation(), DefaultEdgeExtension());
  }


} } // namespaces GPU, vw


#endif
