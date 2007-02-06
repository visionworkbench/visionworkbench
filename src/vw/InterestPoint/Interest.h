#ifndef _INTEREST_POINT_INTEREST_H_
#define _INTEREST_POINT_INTEREST_H_

#include <vector>
#include <stdio.h>

#include <vw/Image/Manipulation.h>
#include <vw/Image/Algorithms.h>
#include <vw/InterestPoint/InterestData.h>
#include <vw/FileIO.h>

namespace vw { namespace ip {

  /// Returns "cornerness" image, where the local maxima correspond to corners.
  /// By default uses Noble measure of corner strength (requires no tuning).
  /// Also supports Harris measure if positive k is specified (typical values:
  /// 0.04 <= k <= 0.15).
  template <class T>
  void harris_interest(ImageInterestData<T> const& data, T k = -1.0, T scale = 1.0) {
    typedef ImageView<T> Image;

    // Calculate elements of Harris matrix
    std::vector<T> kernel;
    generate_gaussian_kernel(kernel, scale, 0);
    printf("Rows=%i, Cols=%i\n", data.grad_x.rows(), data.grad_x.cols());
    printf("scale=%f, kernel size=%i... ", scale, kernel.size());
    Image Ix2 = separable_convolution_filter(data.grad_x * data.grad_x,
					     kernel, kernel);
    Image Iy2 = separable_convolution_filter(data.grad_y * data.grad_y,
					     kernel, kernel);
    Image Ixy = separable_convolution_filter(data.grad_x * data.grad_y,
					     kernel, kernel);

    // Estimate "cornerness"
    printf("Estimating cornerness... ");
    Image trace = Ix2 + Iy2;
    Image det = Ix2 * Iy2 - Ixy * Ixy;
    if (k < 0) {
      // Noble measure (preferred)
      data.interest = det / (trace + 0.000001);
    } else {
      // Standard Harris corner measure
      data.interest = det - k * trace * trace;
    }
    printf("Leaving\n");
  }

  template <class T>
  ImageView<T> harris_interest(ImageView<T> const& image, T k = -1.0, T scale = 1.0) {
    ImageInterestData<T> data(image);
    harris_interest(data, k, scale);
    return data.interest;
  }

  // Scale-normalized Laplacian of Gaussian interest measure
  template <class T>
  void log_interest(ImageInterestData<T> const& data, T scale = 1.0) {
    data.interest = scale * laplacian_filter(data.src);
  }

  template <class T>
  ImageView<T> log_interest(ImageView<T> const& image, T scale = 1.0) {
    ImageInterestData<T> data(image);
    log_interest(data, scale);
    return data.interest;
  }

  template <class T>
  class InterestBase {
  protected:
    PeakType type;

  public:
    virtual int compute_interest(ImageInterestData<T> const& data, T scale = 1.0) {
      return 0;
    }

    PeakType peak_type() { return type; }
  };

  template <class T>
  class HarrisInterest : public InterestBase<T> {
  protected:
    T k;
    T v2; // Relative integration scale parameter (squared)

  public:
    HarrisInterest(T k_in = -1.0, T v2_in = 2.0) : k(k_in), v2(v2_in) { this->type = IP_MAX; }

    virtual int compute_interest(ImageInterestData<T> const& data, T scale = -1.0) {
      if (scale < 0)
	harris_interest(data, k);
      else
	harris_interest(data, k, scale / v2);
      return 0;
    }
  };

  template <class T>
  class LoGInterest : public InterestBase<T> {
  public:
    LoGInterest() { this->type = IP_MINMAX; }

    virtual int compute_interest(ImageInterestData<T> const& data, T scale = 1.0) {
      log_interest(data, scale);
      return 0;
    }
  };


/* // Code graveyard - CRTP experiment
  template <class ImplT>
  struct InterestBase {
    typedef typename ImplT::value_type valT;

    /// Returns the derived implementation type.
    ImplT& impl() { return *static_cast<ImplT*>(this); }

    /// Returns the derived implementation type.
    ImplT const& impl() const { return *static_cast<ImplT const*>(this); }

    int compute_interest(const ImageInterestData<valT>& data,
                         valT scale = 1.0) {
        return impl().compute_interest(data, scale);
    }

    inline PeakType get_peak_type() {
        return ImplT::peak_type;
    }
  };

  template <class T>
  class HarrisInterest : public InterestBase<HarrisInterest<T> > {
  private:
    T k;
    T v2; // Relative integration scale parameter (squared)

  public:
    typedef T value_type;
    static const PeakType peak_type = IP_MAX;

    HarrisInterest(T k_in = -1.0, T v2_in = 2.0) : k(k_in), v2(v2_in) {}

    int compute_interest(ImageInterestData<T> const& data, T scale = -1.0) {
      if (scale < 0)
	harris_interest(data, k);
      else
	harris_interest(data, k, scale / v2);
      return 0;
    }
  };
*/
} } //namespace vw::ip

#endif
