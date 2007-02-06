#ifndef __INTERESTPOINT_SIFTDESCRIPTOR_H__
#define __INTERESTPOINT_SIFTDESCRIPTOR_H__

#include <vw/InterestPoint/Descriptor.h>
#include <vw/InterestPoint/Detector.h>
#include <vw/InterestPoint/WeightedHistogram.h>

namespace vw { 
namespace ip {

  template <class T>
  class SIFT_Descriptor : public DescriptorBase<SIFT_Descriptor<T> >
  {
    ImageInterestData<T> data;

  public:
    typedef T real_type;

    void set_source(const ImageView<T>& src) {
      data.set_source(src);
    }

    int compute_descriptor_from_support(InterestPoint& pt,
					const ImageView<T>& ori_region,
					const ImageView<T>& mag_region) {
      pt.descriptor.set_size(128);
      T histograms[4][4][8];
      T bin_size = M_PI / (T)4;

      // Apply gaussian weighting to magnitudes
      ImageView<T> weight_region, kernel;
      make_gaussian_kernel_2d(kernel, 9.0, 16);
      weight_region = mag_region * kernel;

      for (int i = 0; i < 16; i++) {
	for (int j = 0; j < 16; j++) {
	  T ori = ori_region(i, j);
	  T weight = weight_region(i, j);

	}
      }

      pt.descriptor = normalize(pt.descriptor);

      return 0;
    }

    int compute_descriptors( std::vector<InterestPoint>& points,
			     const ImageView<T>& source ) {
      ImageView<T> ori_region;
      ImageView<T> mag_region;

      set_source(source);
      for (unsigned i = 0; i <points.size(); i++){
	get_support(ori_region, points[i], data.ori);
	get_support(mag_region, points[i], data.mag);
	
	// This function is specialized for different descriptors.
	compute_descriptor_from_support( points[i], ori_region, mag_region);
      }
      return 0;
    }

    int get_support(ImageView<T>& support,
		    const InterestPoint& pt,
		    const ImageView<T>& source) {
      vw::ip::get_support(support, pt, source, 16); //16x16 sample
    }
  };

}} // namespace vw::ip

#endif //__INTERESTPOINT_SIFTDESCRIPTOR_H__
