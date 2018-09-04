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


#ifndef __GAMMAMIXTURECOMPONENT_H__
#define __GAMMAMIXTURECOMPONENT_H__

#include <vw/Stereo/MixtureComponent.h>
#include <vw/Stereo/AffineMixtureComponent.h>

#include <vw/Image/ImageViewRef.h>
#include <limits>

namespace vw {
namespace stereo {

  template<class PixelT, class PrecisionT>
  class GammaMixtureComponent : public MixtureComponentBase<GammaMixtureComponent<PixelT, PrecisionT> > {
  private:
    double digamma(double k) {
      if(k >= 8)
        return log(k) - (1 + (1 - (1/10. - 1/(21*k*k))/(k*k))/(6*k))/(2*k);
      else
        return digamma(k+1) - 1/k;
    }

    double trigamma(double k) {
      if(k >= 8)
        return (1 + (1 + (1 - (1/5. - 1/(7*k*k))/(k*k))/(3*k))/(2*k))/k;
      else
        return trigamma(k+1) + 1/(k*k);
    }

  public:
    template<class ImageT>
    GammaMixtureComponent(ImageViewBase<ImageT> const& image,
                          Vector2i window_size,
                          PrecisionT k_0, PrecisionT theta_0,
                          PrecisionT pk_min, PrecisionT ptheta_min) : image_r(image.impl()), image_c(window_size(0), window_size(1)),
      _k(k_0), k0(k_0),
      _theta(theta_0), theta0(theta_0), k_min(pk_min), theta_min(ptheta_min),
      w(window_size(0), window_size(1)), w_gaussian(window_size(0), window_size(1)), w_local(image.impl().cols(), image.impl().rows()),
      p(window_size(0), window_size(1)) {

      fill(w, 1.);
      fill(p, .5);
      l_likelihood = INFINITY;

      k_max = 200;
      //k_min = 1.; //TODO: set these limits in a less adhoc way
      var_min = 1e-3; // 1e-2 worked best

      image_r = clamp(image_r, 1e-14, INFINITY);
      image_c = clamp(crop(transform(edge_extend(image_r, ZeroEdgeExtension()), T), window), 1e-14, INFINITY);
      l_set = false;
    }

    void reset(BBox2i left_window, ImageView<PrecisionT> const& g) {
      fill(w, 1.);
      fill(p, .5);
      _k = k0;
      _theta = theta0;
      window = left_window;
      image_c = clamp(crop(transform(edge_extend(image_r, ZeroEdgeExtension()), T), window), 1e-14, INFINITY);
      w_gaussian = g;
    }

    void fit_parameters() { // weights only used here, so transform the weights into the correct space only here
      BBox2i image_size(0, 0, image_r.cols(), image_r.rows());
      BBox2i bbox = compute_transformed_bbox_fast(translate(edge_extend(w, ZeroEdgeExtension()), window.min().x(), window.min().y()).get_size(), inverse(T));
      bbox.max().x() += 1; //for some reason compute transformed bbox outputs an inclusive bottom-right boundry point....
      bbox.max().y() += 1;
      bbox += inverse(T).forward(window.min());
      bbox.crop(image_size);

      /*
        ImageView<PrecisionT> w_adj = crop(transform(crop(translate(edge_extend(w, ZeroEdgeExtension()), window.min().x(), window.min().y()),
                                                          image_size),
                                                     inverse(T)),
                                           bbox);
      */
      ImageViewRef<PrecisionT> w_adj = w;

      PrecisionT sum_weights_adj = sum_of_pixel_values(w_adj);

      if(sum_weights_adj < 1e-10)
        return;

        PrecisionT k_last = _k;
        PrecisionT theta_last = _theta;
        double l_likelihood_old = -sum_of_pixel_values(w*((_k-1)*log(image_c) - image_c/_theta)) + sum_weights_adj*(_k*log(_theta) + lgamma(_k));

        PrecisionT mu = 0, mu_log = 0;
        mu = sum_of_pixel_values(w*image_c)/sum_weights_adj;
        mu_log = sum_of_pixel_values(w*log(image_c))/sum_weights_adj;

        PrecisionT s = log(mu) - mu_log;
        if(s <= 0) { // this is to keep things well conditioned; the equation log(k) - digamma(k) = s cannot be solved for s <= 0
          s = 1e-6;
          //std::cout << "s negative! K became too low!" << std::endl;
        }
        _k = (3 - s + sqrt(pow(s-3, 2) + 24*s))/(12*s);
        //_k = std::max<PrecisionT>(_k, k_min); // threshold
        //_k = std::min<PrecisionT>(_k, k_max); // threshold
        if(_k != _k) {
          _k = k_last;
        }

        // now use newton-raphson to get true MLE estimates
        double delta = std::numeric_limits<double>::max();
        int iter = 0;
        while(delta > 1e-7 && iter < 5) {
          double k_old = _k;
          _k = _k - (log(_k) - digamma(_k) - s)/(1/_k - trigamma(_k));
          delta = fabs(k_old - _k);
          iter++;
        }
        _k = std::max<PrecisionT>(_k, k_min); // threshold
        _k = std::min<PrecisionT>(_k, k_max); // threshold
        _theta = mu/_k;
        _theta = std::max<PrecisionT>(_theta, var_min/sqrt(_k)); // threshold

        /*
        if(_theta != _theta) {
          _theta = theta_last;
        }
        */
        //_theta = std::max<PrecisionT>(_theta, theta_min); // threshold



        l_likelihood = -sum_of_pixel_values(w*((_k-1)*log(image_c) -
                                               image_c/_theta))
          + sum_weights_adj*(_k*log(_theta) + lgamma(_k));

        l_set = true;

        if(l_likelihood_old /*+ 1e-10*/ < l_likelihood) {
          //std::cout << "LIKELIHOOD WENT UP IN Gamma by " << l_likelihood - l_likelihood_old<< std::endl;
          //std::cout << "reverting k: " << _k << " -> " << k_last << std::endl;
          //std::cout << "reverting theta: " << _theta << " -> " << theta_last << std::endl;
          //std::cout << "s = " << s << std::endl;
          _k = k_last;
          _theta = theta_last;
          l_likelihood = -sum_of_pixel_values(w*((_k-1)*log(image_c) -
                                                 image_c/_theta))
            + sum_weights_adj*(_k*log(_theta) + lgamma(_k));
        }

        /*
        if(l_likelihood > std::numeric_limits<double>::max() || l_likelihood != l_likelihood) {
          std::cout << "l is bad: " << l_likelihood << std::endl;
          std::cout << "k: " << _k << std::endl;
          std::cout << "theta: " << _theta << std::endl;
          std::cout << "s: " << s << std::endl;
          std::cout << "mu: " << mu << std::endl;
          std::cout << "mu_log: " << mu_log << std::endl;
          std::cout << "sum_weights_adj: " << sum_weights_adj << std::endl;
        }
        */
      }

      void print_status(std::string const& prefix) {
        std::cout << prefix << "k = " << _k << std::endl;
        std::cout << prefix << "theta = " << _theta << std::endl;
      }

      void update_posterior() {
        p = exp((_k-1)*log(image_c) - image_c/_theta - lgamma(_k) - _k*log(_theta));
      }

      void set_affine_transform(AffineTransformOrigin const& T_new) {
        T = T_new;
        image_c = clamp(crop(transform(edge_extend(image_r, ZeroEdgeExtension()), T), window), 1e-14, INFINITY);
      }

      PrecisionT log_likelihood() {
        if(!l_set) {
          PrecisionT sum_weights = sum_of_pixel_values(w); // this is the normalization constant for the weights
          l_likelihood = -sum_of_pixel_values(w*((_k-1)*log(image_c) -
                                                 image_c/_theta))
            + sum_weights*(_k*log(_theta) + lgamma(_k));
          l_set = true;
        }
        return l_likelihood;
      }

      PrecisionT k() { return _k; }
      PrecisionT theta() { return _theta; }

      ImageView<PrecisionT> & weights() { return w; }
      ImageView<PrecisionT> const& weights() const { return w; }
      ImageView<PrecisionT> const& prob() const { return p; }

    private:
      ImageViewRef<PixelT> image_r;
      ImageView<PixelT> image_c; // clipped image

      BBox2i window;

      PrecisionT _k, k0;
      PrecisionT _theta, theta0;
      PrecisionT k_min, theta_min;
      PrecisionT var_min;
      PrecisionT k_max;
      PrecisionT l_likelihood;

      ImageView<PrecisionT> w;  // weights
      ImageView<PrecisionT> w_gaussian;  // gaussian window
      ImageView<PrecisionT> err;  // errors
      ImageView<PrecisionT> w_local; // weights in local coordinate frame
      ImageView<PrecisionT> p; // posterior probabilities of data in window

      AffineTransformOrigin T;
      bool l_set;
    };

}}

#endif
