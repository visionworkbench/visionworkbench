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


/// \file ImageComposite.h
///
/// A view class that represents a composite mosaic of images.
///
#ifndef __VW_MOSAIC_DISKIMAGEPYRAMID_H__
#define __VW_MOSAIC_DISKIMAGEPYRAMID_H__

#include <iostream>
#include <list>
#include <set>
#include <vector>
#include <cstdlib>

#include <vw/Core/ProgressCallback.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Cartography/GeoReferenceUtils.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/DiskImageResourceGDAL.h>
#include <vw/FileIO/FileUtils.h>
#include <vw/FileIO/DiskImageUtils.h>
#include <vw/Image/AntiAliasing.h>
#include <vw/Image/MaskViews.h>
#include <vw/Math/Statistics.h>

namespace vw { namespace mosaic {

  namespace fs = boost::filesystem;


  // TODO: Move to file utils!
  /// Get a filename by simply replacing current extension with given suffix
  inline std::string filename_from_suffix1(std::string const& input_file,
                                           std::string const& suffix){
    std::string prefix = vw::prefix_from_filename(input_file);
    std::string output_file = prefix + suffix;
    return output_file;
  }

  /// Get a filename by  replacing current extension with given suffix,
  /// and making the file be in the current directory.
  inline std::string filename_from_suffix2(std::string const& input_file,
                                           std::string const& suffix){
    std::string prefix = vw::prefix_from_filename(input_file);
    boost::filesystem::path p(input_file);
    prefix = p.stem().string();
    std::string output_file = prefix + suffix;
    return output_file;
  }

  // Compile-intensive: read_channels results in 12 template instantiations.
  // Single-channel images are read into Image<double>, while
  // multi-channel in Image<Vector>, and we skip reading extra channels.
  // Read a file, or one of its GDAL overview levels (overview_level >= 0), as a
  // single-channel image. A negative overview_level reads the full resolution.
  template<class PixelT>
  typename boost::enable_if<boost::mpl::or_<boost::is_same<PixelT,double>,
                            boost::is_same<PixelT,float>>,
                            ImageViewRef<PixelT>>::type
  custom_read(std::string const& file, int overview_level = -1){
    int num_channels = get_num_channels(file);
    if (overview_level < 0) {
      // Bugfix: Sometimes the image may actually have multiple channels
      // but we choose to read one channel only in the caller of this
      // function. In that case use read_channels() rather than
      // DiskImageView(). The latter would cause a crash.
      if (num_channels == 1)
        return DiskImageView<PixelT>(file);
      else
        return vw::select_channel(vw::read_channels<1, PixelT>(file, 0), 0);
    }
    // Read the given overview, opening it as if it were the full image.
    boost::shared_ptr<vw::DiskImageResourceGDAL>
      rsrc(new vw::DiskImageResourceGDAL(file, overview_level));
    if (num_channels == 1)
      return DiskImageView<PixelT>(rsrc);
    if (num_channels == 2)
      return vw::select_channel(DiskImageView<vw::Vector<PixelT,2>>(rsrc), 0);
    if (num_channels == 3)
      return vw::select_channel(DiskImageView<vw::Vector<PixelT,3>>(rsrc), 0);
    if (num_channels == 4)
      return vw::select_channel(DiskImageView<vw::Vector<PixelT,4>>(rsrc), 0);
    vw_throw(vw::NoImplErr()
             << "Reading an overview of an image with more than 4 "
             << "non-uint8 bands is not supported.\n");
    return DiskImageView<PixelT>(rsrc); // not reached
  }
  // The multi-channel (vector pixel) version.
  template<class PixelT>
  typename boost::disable_if<boost::mpl::or_<boost::is_same<PixelT,double>,
                             boost::is_same<PixelT,float>>,
                             ImageViewRef<PixelT>>::type
  custom_read(std::string const& file, int overview_level = -1){
    if (overview_level < 0)
      return vw::read_channels<vw::math::VectorSize<PixelT>::value, typename PixelT::value_type>(file, 0);
    // Read the given overview. It preserves the band count of the full image.
    boost::shared_ptr<vw::DiskImageResourceGDAL>
      rsrc(new vw::DiskImageResourceGDAL(file, overview_level));
    return DiskImageView<PixelT>(rsrc);
  }

  // TODO: Clean up!
  // Gets called for scalar pixel types (double, float, uint8)
  template<class PixelT>
  typename boost::enable_if< boost::mpl::or_< boost::is_same<PixelT,double>,
                             boost::mpl::or_< boost::is_same<PixelT,float>,
                             boost::is_same<PixelT,vw::uint8> > >,
                             ImageViewRef< PixelMask<PixelT> > >::type
  create_custom_mask(ImageViewRef<PixelT> & img, double nodata_val){
    return create_mask(img, nodata_val);
  }
  // Gets called for compound pixel types (Vector<u8, N>)
  template<class PixelT>
  typename boost::disable_if< boost::mpl::or_< boost::is_same<PixelT,double>,
                              boost::mpl::or_< boost::is_same<PixelT,float>,
                              boost::is_same<PixelT,vw::uint8> > >,
                              ImageViewRef< PixelMask<PixelT> > >::type
  create_custom_mask(ImageViewRef<PixelT> & img, double nodata_val){
    // No nodata was provided (nodata_val is NaN): mask nothing. Casting NaN to
    // an integer channel type is undefined (in practice yields 0), which would
    // wrongly mask pure-black pixels.
    if (std::isnan(nodata_val))
      return pixel_cast<PixelMask<PixelT>>(img);
    PixelT mask_pixel;
    mask_pixel.set_all(nodata_val);
    return create_mask(img, mask_pixel);
  }

  // For scalar pixel types (double, float, uint8): use range masking if valid
  // range is provided, otherwise use nodata masking.
  template<class PixelT>
  typename boost::enable_if<boost::mpl::or_<boost::is_same<PixelT, double>,
                            boost::mpl::or_<boost::is_same<PixelT, float>,
                            boost::is_same<PixelT, vw::uint8>>>,
                            ImageViewRef<PixelMask<PixelT>>>::type
  maskForPyramid(ImageViewRef<PixelT> & img, double nodata_val,
                 float valid_min, float valid_max) {
    if (!std::isnan(valid_min) && !std::isnan(valid_max))
      return create_mask(img, PixelT(valid_min), PixelT(valid_max));
    return create_mask(img, nodata_val);
  }

  // For compound pixel types (Vector<u8, N>): always use nodata masking.
  // Range masking is not supported for multi-channel pixels.
  template<class PixelT>
  typename boost::disable_if<boost::mpl::or_<boost::is_same<PixelT, double>,
                             boost::mpl::or_<boost::is_same<PixelT, float>,
                             boost::is_same<PixelT, vw::uint8>>>,
                             ImageViewRef<PixelMask<PixelT>>>::type
  maskForPyramid(ImageViewRef<PixelT> & img, double nodata_val,
                 float /*valid_min*/, float /*valid_max*/) {
    // No nodata was provided (nodata_val is NaN): mask nothing. Casting NaN to
    // an integer channel type is undefined (in practice yields 0), which would
    // wrongly mask pure-black pixels out of the anti-aliased pyramid average.
    if (std::isnan(nodata_val))
      return pixel_cast<PixelMask<PixelT>>(img);
    PixelT mask_pixel;
    mask_pixel.set_all(nodata_val);
    return create_mask(img, mask_pixel);
  }

  // TODO: Cleanup!
  /// If output_file exists, is not older than input_file,
  /// and has given numbers of rows and columns, don't overwrite it.
  inline bool overwrite_if_no_good(std::string const& input_file,
                                   std::string const& output_file,
                                   int cols=-1, int rows=-1){

    fs::path input_path(input_file);
    std::time_t input_time = fs::last_write_time(input_path);

    bool overwrite = true;
    if (fs::exists(output_file)) {
      try{
        DiskImageView<double> curr(output_file);
        fs::path output_path(output_file);
        std::time_t output_time = fs::last_write_time(output_path);
        if (output_time < input_time){
          overwrite = true; // too old
          return overwrite;
        }

        if (cols >= 0 && rows >= 0) {
          if (curr.cols() == cols && curr.rows() == rows){
            overwrite = false;
          }else{
            overwrite = true;
          }
        }else{
          overwrite = false;
        }

      }catch(...){
        overwrite = true;
      }
    }

    return overwrite;
  }

  // Find approximate min/max for valid pixels, ignoring outliers.
  // Uses the lowest pyramid level. Only used for single-channel images
  // that need [min, max] to [0, 255] rescaling for display. Multi-channel
  // uint8 images map directly to RGB without rescaling.
  // The file is read as double regardless of on-disk type, with no rescaling,
  // so nodata comparisons (which use double) are exact.
  inline vw::Vector2
  approx_bounds_nocache(std::string const& file, bool has_nodata,
                        double nodata_val,
                        float valid_min = std::numeric_limits<float>::quiet_NaN(),
                        float valid_max = std::numeric_limits<float>::quiet_NaN(),
                        int overview_level = -1) {

    double big = std::numeric_limits<double>::max();
    bool has_valid_range = !std::isnan(valid_min) && !std::isnan(valid_max);

    // Sample the coarsest overview when one is given, so this stays fast.
    boost::shared_ptr<vw::DiskImageResource> rsrc;
    if (overview_level < 0)
      rsrc = vw::DiskImageResourcePtr(file);
    else
      rsrc.reset(new vw::DiskImageResourceGDAL(file, overview_level));
    rsrc->set_rescale(false);
    DiskImageView<double> img(rsrc);

    double num_samples = 250.0; // too many samples makes this slow
    int delta_col = (int)std::max(ceil(img.cols() / num_samples), 1.0);
    int delta_row = (int)std::max(ceil(img.rows() / num_samples), 1.0);
    std::vector<double> vals;
    vals.reserve(num_samples * num_samples);
    for (int col = 0; col < img.cols(); col += delta_col) {
      for (int row = 0; row < img.rows(); row += delta_row) {
        double v = img(col, row);
        if (std::isnan(v)) continue;
        if (has_nodata && v == nodata_val) continue;
        if (has_valid_range && (v < valid_min || v > valid_max)) continue;
        vals.push_back(v);
      }
    }

    vw::Vector2 bounds(-big, big);
    if (vals.empty()) return bounds;

    // Find the bounds using percentiles
    double pct = 0.01; // 1% outliers at either end
    double outlier_factor = 4;
    double b = 0.0, e = 0.0;
    vw::math::find_outlier_brackets(vals, pct, outlier_factor, b, e);

    // Tighten the bounds
    std::sort(vals.begin(), vals.end());
    b = std::max(b, vals[0]);
    e = std::min(e, vals[vals.size()-1]);

    // If the bounds are the same, expand them
    if (b == e) {
      b -= 1;
      e += 1;
    }

    return vw::Vector2(b, e);
  }
  
  /// A class to manage very large images and their subsampled
  /// versions in a pyramid. The most recently accessed tiles are
  /// cached in memory. Caching is handled by use of the
  /// DiskImageView class. Constructing this class creates a temporary
  /// file on disk for each level of the pyramid.
  template <class PixelT>
  class DiskImagePyramid {

  public:
    typedef typename DiskImageView<PixelT>::pixel_type pixel_type;

    // Constructor. Note that we use NaN as nodata if not available,
    // that has the effect of not accidentally setting some pixels to nodata.
    // If valid_min and valid_max are not NaN, pixels outside that range
    // are treated as nodata (needed for ISIS .cub special pixels).
    DiskImagePyramid(std::string const& base_file = "",
		     vw::GdalWriteOptions const& opt = vw::GdalWriteOptions(),
		     int lowest_resolution_subimage_num_pixels = 1000*1000,
		     int subsample = 2,
		     float valid_min = std::numeric_limits<float>::quiet_NaN(),
		     float valid_max = std::numeric_limits<float>::quiet_NaN());

    // Given a region (at full resolution) and a scale factor, compute
    // the portion of the image in the region, subsampled by a factor no
    // more than the input scale factor.  Also return the precise subsample
    // factor used and the region at that scale level.
    void get_image_clip(double scale_in, BBox2i region_in,
			ImageView<PixelT> & clip, double & scale_out, BBox2i & region_out) const;

    // Find the right pyramid level to use for the given sub scale
    int pyramidLevel(double sub_scale) const;
    
    vw::Vector2 approx_bounds() const {
      return m_approx_bounds;
    }
    
    ~DiskImagePyramid() {}

    // These all describe the highest resolution pyramid layer
    int32 cols  () const { return m_pyramid[0].cols(); }
    int32 rows  () const { return m_pyramid[0].rows(); }
    int32 planes() const { return m_pyramid[0].planes(); }

    std::vector<int> const& scales() const { return m_scales; }

    std::vector<ImageViewRef<PixelT>> const& pyramid() const { return m_pyramid; }
    
    double get_nodata_val() const { return(m_nodata_val); }

    /// Return the highest resolution pyramid layer
    ImageViewRef<PixelT>        bottom()       { return m_pyramid[0]; }
    ImageViewRef<PixelT> const& bottom() const { return m_pyramid[0]; }

    std::set<std::string> const& get_temporary_files() const {return m_temporary_files;}

  private:

    vw::GdalWriteOptions m_opt;

    // The subsample factor to go to the next level of the pyramid (must be >= 2).
    int m_subsample;

    // The maxiumum number of pixels in the coarsest level of the pyramid
    // (keep on downsampling until getting to this number or under it).
    int m_lowest_resolution_subimage_num_pixels;

    //  The pyramid. Largest images come earlier.
    std::vector<ImageViewRef<PixelT>> m_pyramid;

    // The files (stored on disk) containing the images in the pyramid.
    std::vector<std::string> m_pyramid_files;

    // We may wipe these at the end
    std::vector<std::string> m_cached_files;

    double m_nodata_val;
    float m_valid_min;
    float m_valid_max;
    std::vector<int> m_scales;
    
    /// Contains all the temporary image files we create
    std::set<std::string> m_temporary_files;

    vw::Vector2 m_approx_bounds;
  };

  //#################################################################################
  // Function definitions

  template <class PixelT>
  DiskImagePyramid<PixelT>::DiskImagePyramid(std::string const& base_file,
                                             vw::GdalWriteOptions const& opt,
                                             int lowest_resolution_subimage_num_pixels,
                                             int subsample,
                                             float valid_min,
                                             float valid_max):
    m_opt(opt), m_subsample(subsample),
    m_lowest_resolution_subimage_num_pixels(lowest_resolution_subimage_num_pixels),
    m_nodata_val(std::numeric_limits<double>::quiet_NaN()),
    m_valid_min(valid_min), m_valid_max(valid_max) {

    if (base_file.empty())
      return;

    if (subsample < 2)
      vw_throw( ArgumentErr()
                << "Must subsample by a factor of at least 2.\n");

    if (lowest_resolution_subimage_num_pixels < 4)
      vw_throw( ArgumentErr()
                << "The lowest resolution subimage must be at least 2x2 in size.\n");

    m_pyramid.push_back(custom_read<PixelT>(base_file));

    m_pyramid_files.push_back(base_file);
    m_scales.push_back(1);

    // Get the nodata value, if present
    bool has_nodata = read_nodata_val(base_file, m_nodata_val);

    // Decimation factors (2, 4, 8, ...), until the coarsest level is small enough.
    std::vector<int> factors;
    int scale = 1;
    while (double(m_pyramid[0].cols()/scale) * double(m_pyramid[0].rows()/scale)
           > m_lowest_resolution_subimage_num_pixels) {
      scale *= subsample;
      factors.push_back(scale);
    }

    if (!factors.empty()) {
      vw_out() << "Construct an image pyramid for: " << base_file << "\n";

      // Build overviews with GDAL (gdaladdo), an external .ovr sidecar. GDAL
      // preserves the type, nodata, and mask (so ISIS specials are handled).
      std::string proxy_dir
        = (fs::temp_directory_path() / "asp_pyramids").string();
      try {
        fs::create_directories(proxy_dir);
      } catch (...) {}

      bool ok = cartography::build_overviews(base_file, factors, "AVERAGE",
                                             proxy_dir);
      if (ok) {
        // Each coarser level reads the matching overview of the same file.
        for (size_t k = 0; k < factors.size(); k++) {
          m_pyramid.push_back(custom_read<PixelT>(base_file, (int)k));
          m_pyramid_files.push_back(base_file);
          m_scales.push_back(factors[k]);
        }
        // The overview sidecar is a temporary file, removable on exit.
        m_temporary_files.insert(base_file + ".ovr");
      } else {
        // No overviews: display reads the full image at all zoom levels (slower).
        factors.clear();
      }
    }

    // This is expensive, so cache it going forward. Sample the coarsest level.
    int coarsest_ov = (int)factors.size() - 1; // -1 when there are no overviews
    m_approx_bounds = approx_bounds_nocache(base_file,
                                            has_nodata, m_nodata_val,
                                            m_valid_min, m_valid_max,
                                            coarsest_ov);
  }

  // Find the right pyramid level to use for the given sub scale
  template <class PixelT>
  int DiskImagePyramid<PixelT>::pyramidLevel(double sub_scale) const {
    int level = 0;
    while (1) {
      if (level+1 >= (int)m_scales.size()) break; // last level
      if (m_scales[level+1] > sub_scale)  break; // too coarse
      level++;
    }

    return level;
  }
  
  template <class PixelT>
  void DiskImagePyramid<PixelT>::get_image_clip(double scale_in, BBox2i region_in,
                                                ImageView<PixelT> & clip,
                                                double & scale_out, BBox2i & region_out) const {

    if (m_pyramid.empty())
      vw_throw( ArgumentErr() << "Uninitialized image pyramid.\n");

    // Find the right pyramid level to use
    int level = pyramidLevel(scale_in);

    // Debug aid (gated): print which pyramid level serves each display scale.
    static bool debug_pyramid = (getenv("ASP_PYRAMID_DEBUG") != NULL);
    if (debug_pyramid)
      std::cout << "ASP_PYRAMID_DEBUG scale_in=" << scale_in
                << " level=" << level
                << " scale=" << m_scales[level]
                << " file=" << m_pyramid_files[level] << "\n";

    region_in.crop(bounding_box(m_pyramid[0]));
    scale_out  = m_scales[level];
    region_out = region_in/scale_out;
    region_out.crop(bounding_box(m_pyramid[level]));

    clip = crop(m_pyramid[level], region_out);
  }

}} // End namespace vw::mosaic

#endif // __VW_MOSAIC_DISKIMAGEPYRAMID_H__
