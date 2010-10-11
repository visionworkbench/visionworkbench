// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_STEREO_GAUSSIAN_MIXTURE_COMPONENT__
#define __VW_STEREO_GAUSSIAN_MIXTURE_COMPONENT__

#include <vw/Stereo/MixtureComponent.h>
#include <vw/Stereo/AffineMixtureComponent.h>

#include <vw/Image.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Math.h>

namespace vw {
namespace stereo {

  template<class PixelT, class PrecisionT>
  class GaussianMixtureComponent : public MixtureComponentBase<GaussianMixtureComponent<PixelT, PrecisionT> > {
  public:
    template<class ImageT>
    GaussianMixtureComponent(ImageViewBase<ImageT> const& image,
                             Vector2i window_size,
                             PrecisionT mu_0, PrecisionT sigma_0, PrecisionT sigma_min) : image_r(image.impl()), image_c(window_size(0), window_size(1)),
      m(mu_0), m0(mu_0), s(sigma_0), s0(sigma_0), s_min(sigma_min), w(window_size(0), window_size(1)), w_local(image.impl().cols(), image.impl().rows()),
      p(window_size(0), window_size(1)) {
      sqrt_2_pi = sqrt(2*M_PI);
      fill(w, 1.);
      fill(p, .5);
      l_likelihood = INFINITY;
      image_c = crop(transform(edge_extend(image_r, ZeroEdgeExtension()), T), window);
    }

    void reset(BBox2i left_window,
               ImageView<PrecisionT> const& g) {
      fill(w, 1.);
      fill(p, .5);
      s = s0;
      m = m0;
      window = left_window;
      image_c = crop(transform(edge_extend(image_r, ZeroEdgeExtension()), T), window);
      //w_gaussian = g;
    }

    void print_status(std::string const& prefix) {
      std::cout << prefix << "mu = " << m << std::endl;
      std::cout << prefix << "sigma = " << s << std::endl;
    }

    void fit_parameters() {
      // weights only used here, so transform the weights into the
      // correct space only here

      //BBox2i image_size(0, 0, image_r.cols(), image_r.rows());
      //BBox2i bbox = compute_transformed_bbox_fast(translate(w, window.min().x(), window.min().y()), inverse(T));
      //bbox += inverse(T).forward(window.min());
      //std::cout << "1: " << bbox << std::endl;
      //bbox.crop(image_size);

      /*
        ImageView<PrecisionT> w_adj = crop(transform(crop(translate(edge_extend(w, ZeroEdgeExtension()), window.min().x(), window.min().y()),
        image_size),
        inverse(T)),
        bbox);
      */
      PrecisionT sum_weights_adj = sum_of_pixel_values(w);
      //PrecisionT sum_weights = sum_of_pixel_values(crop(transform(w, T), window));

      m = sum_of_pixel_values(w*image_c)/sum_weights_adj;
      if(m != m) { //TODO: handle this better?
          /*
          ImageWindow win = create_window("weights");

          show_image(win, w_adj);
          usleep(1*1000*1000);
          destroy_window(win);



          std::cout << "bbox: " << bbox << std::endl;
          std::cout << "image_r: " << sum_of_pixel_values(crop(image_r, bbox)) << std::endl;
          std::cout << "sum_w: " << sum_of_pixel_values(w) << std::endl;
          std::cout << "sum_w_adj: " << sum_of_pixel_values(w_adj) << std::endl;
          std::cout << m << " -> " << m0 << std::endl;
          */
        m = m0;
      }

      err = image_c - m;
      s = sqrt(sum_of_pixel_values(w*pow(err, 2))/sum_weights_adj);
      if(s != s) //TODO: handle this better?
        s = s0;

      s = std::max<PrecisionT>(s, s_min); // threshold the variance so it never gets too small

      l_likelihood = sum_of_pixel_values(.5*w*pow(err/s,2)) + sum_weights_adj*(log(s) + log(sqrt_2_pi));

    }

    void update_posterior() {
      p = exp(-.5*pow((image_c - m)/s, 2))/(s*sqrt_2_pi); // gaussian posterior density
    }

    void set_affine_transform(AffineTransformOrigin const& T_new) {
      T = T_new;
      image_c = crop(transform(edge_extend(image_r, ZeroEdgeExtension()), T), window);
    }

    PrecisionT log_likelihood() { return l_likelihood; }
    PrecisionT mu() { return m; }
    PrecisionT sigma() { return s; }

    ImageView<PrecisionT> const& errors() const { return err; }
    ImageView<PrecisionT> & weights() { return w; }
    ImageView<PrecisionT> const& weights() const { return w; }

    // TODO: make this more efficient for reading
    ImageView<PrecisionT> prob() const { return p; }


    private:
      ImageViewRef<PixelT> image_r;
      ImageView<PixelT> image_c; // clipped image
      BBox2i window;

      PrecisionT m, m0;
      PrecisionT s, s0;
      PrecisionT s_min;
      PrecisionT sqrt_2_pi;
      PrecisionT l_likelihood;

      ImageView<PrecisionT> w;  // weights
      ImageView<PrecisionT> w_gaussian;  // gaussian window
      ImageView<PrecisionT> err;  // errors
      ImageView<PrecisionT> w_local; // weights in local coordinate frame
      ImageView<PrecisionT> p; // posterior probabilities of data in window

      AffineTransformOrigin T;

    };

}}

#endif//__VW_STEREO_GAUSSIAN_MIXTURE_COMPONENT__
