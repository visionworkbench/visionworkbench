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

#include <vw/Core/ProgressCallback.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Cartography/GeoReferenceUtils.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/FileUtils.h>
#include <vw/FileIO/DiskImageUtils.h>
#include <vw/Image/AntiAliasing.h>
#include <vw/Image/MaskViews.h>
#include <vw/Math/Statistics.h>

#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/type_traits/is_integral.hpp>

#include <cmath>
#include <limits>

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
  // Scalar (arithmetic) pixel types are read into a single-channel image,
  // while compound (Vector) pixel types are read into Image<Vector>, and we
  // skip reading extra channels. Using boost::is_arithmetic (rather than an
  // explicit double/float list) lets integer scalar types such as uint8 and
  // uint16 also build native pyramids.
  template<class PixelT>
  typename boost::enable_if<boost::is_arithmetic<PixelT>,
                            ImageViewRef<PixelT>>::type
  custom_read(std::string const& file){
    // Bugfix: Sometimes the image may actually have multiple channels
    // but we choose to read one channel only in the caller of this
    // function. In that case use read_channels() rather than
    // DiskImageView(). The latter would cause a crash.
    int num_channels = get_num_channels(file);
    if (num_channels == 1)
      return DiskImageView<PixelT>(file);
    else
      return vw::select_channel(vw::read_channels<1, PixelT>(file, 0), 0);
  }
  template<class PixelT>
  typename boost::disable_if<boost::is_arithmetic<PixelT>,
                             ImageViewRef<PixelT>>::type
  custom_read(std::string const& file){
    return vw::read_channels<vw::math::VectorSize<PixelT>::value, typename PixelT::value_type>(file, 0);
  }

  // True when nodata_val cannot be represented exactly in an integral
  // PixelT: NaN (the "no nodata" sentinel), a value outside the type's
  // range (e.g. a uint16 file whose metadata declares nodata -32768), or a
  // fractional value (e.g. nodata 0.5): converting the first two to PixelT
  // is undefined behavior, and truncating a fraction would change which
  // pixels match, masking legitimate pixels either way. Callers treat such
  // a nodata as "no nodata". Floating-point PixelT represents NaN and any
  // nodata directly.
  template<class PixelT>
  typename boost::enable_if<boost::is_integral<PixelT>, bool>::type
  unrepresentableNodata(double nodata_val) {
    // Compare against the exclusive upper bound 2^digits, which is always
    // exactly representable in double; double(max()) rounds UP for 64-bit
    // types, which would let max()+1 through. lowest() (0 or -2^digits) is
    // always exact.
    double upper = std::ldexp(1.0, std::numeric_limits<PixelT>::digits);
    return std::isnan(nodata_val) ||
           nodata_val < double(std::numeric_limits<PixelT>::lowest()) ||
           nodata_val >= upper ||
           nodata_val != std::floor(nodata_val);
  }
  template<class PixelT>
  typename boost::disable_if<boost::is_integral<PixelT>, bool>::type
  unrepresentableNodata(double) { return false; }

  // Build a valid-range mask, converting the floating-point bounds to
  // PixelT safely. For an integral PixelT the bounds are snapped to the
  // integers the comparison actually distinguishes (ceil for the lower,
  // floor for the upper) and clamped into the type's domain, since casting
  // an out-of-range or infinite bound is undefined behavior and truncating
  // a fractional one changes which pixels count as valid. A range wholly
  // outside the domain marks every pixel invalid.
  template<class PixelT>
  typename boost::enable_if<boost::is_integral<PixelT>,
                            ImageViewRef<PixelMask<PixelT>>>::type
  rangeMaskForPyramid(ImageViewRef<PixelT> & img,
                      float valid_min, float valid_max) {
    double lo = std::ceil (double(valid_min));
    double hi = std::floor(double(valid_max));
    double lowest = double(std::numeric_limits<PixelT>::lowest());
    // Exclusive upper bound 2^digits: always exact in double, unlike
    // double(max()) which rounds up for 64-bit types.
    double upper = std::ldexp(1.0, std::numeric_limits<PixelT>::digits);
    if (lo < lowest) lo = lowest;
    if (hi >= upper) hi = std::nextafter(upper, 0.0);
    if (lo > hi) // range does not intersect the type's domain
      return create_mask(img, PixelT(1), PixelT(0));
    return create_mask(img, PixelT(lo), PixelT(hi));
  }
  template<class PixelT>
  typename boost::disable_if<boost::is_integral<PixelT>,
                             ImageViewRef<PixelMask<PixelT>>>::type
  rangeMaskForPyramid(ImageViewRef<PixelT> & img,
                      float valid_min, float valid_max) {
    return create_mask(img, PixelT(valid_min), PixelT(valid_max));
  }

  // TODO: Clean up!
  // Gets called for scalar (arithmetic) pixel types (double, float, uint8,
  // uint16, ...)
  template<class PixelT>
  typename boost::enable_if< boost::is_arithmetic<PixelT>,
                             ImageViewRef< PixelMask<PixelT> > >::type
  create_custom_mask(ImageViewRef<PixelT> & img, double nodata_val){
    // See unrepresentableNodata(): mark all pixels valid rather than convert
    // a NaN or out-of-range nodata to an integer PixelT.
    if (unrepresentableNodata<PixelT>(nodata_val))
      return pixel_cast<PixelMask<PixelT>>(img);
    return create_mask(img, nodata_val);
  }
  // Gets called for compound pixel types (Vector<u8, N>)
  template<class PixelT>
  typename boost::disable_if< boost::is_arithmetic<PixelT>,
                              ImageViewRef< PixelMask<PixelT> > >::type
  create_custom_mask(ImageViewRef<PixelT> & img, double nodata_val){
    PixelT mask_pixel;
    mask_pixel.set_all(nodata_val);
    return create_mask(img, mask_pixel);
  }

  // For scalar (arithmetic) pixel types (double, float, uint8, uint16, ...):
  // use range masking if valid range is provided, otherwise use nodata masking.
  template<class PixelT>
  typename boost::enable_if<boost::is_arithmetic<PixelT>,
                            ImageViewRef<PixelMask<PixelT>>>::type
  maskForPyramid(ImageViewRef<PixelT> & img, double nodata_val,
                 float valid_min, float valid_max) {
    if (!std::isnan(valid_min) && !std::isnan(valid_max))
      return rangeMaskForPyramid(img, valid_min, valid_max);
    // See unrepresentableNodata(): mark all pixels valid rather than convert
    // a NaN or out-of-range nodata to an integer PixelT, which is undefined
    // behavior and masked legitimate pixels (e.g. zeros) out of the average.
    if (unrepresentableNodata<PixelT>(nodata_val))
      return pixel_cast<PixelMask<PixelT>>(img);
    return create_mask(img, nodata_val);
  }

  // For compound pixel types (Vector<u8, N>): always use nodata masking.
  // Range masking is not supported for multi-channel pixels.
  template<class PixelT>
  typename boost::disable_if<boost::is_arithmetic<PixelT>,
                             ImageViewRef<PixelMask<PixelT>>>::type
  maskForPyramid(ImageViewRef<PixelT> & img, double nodata_val,
                 float /*valid_min*/, float /*valid_max*/) {
    PixelT mask_pixel;
    mask_pixel.set_all(nodata_val);
    return create_mask(img, mask_pixel);
  }

  // TODO: Cleanup!
  /// If output_file exists, is not older than input_file,
  /// and has given numbers of rows and columns, don't overwrite it.
  inline bool overwrite_if_no_good(std::string const& input_file,
                                   std::string const& output_file,
                                   int cols=-1, int rows=-1,
                                   ChannelTypeEnum expected_channel_type
                                     = VW_CHANNEL_UNKNOWN,
                                   int expected_num_channels = -1){

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

        // A cached file of the wrong pixel type must be regenerated, or a
        // pyramid whose type changed (e.g. toggling an option that keeps
        // the input's integer type) would silently reuse it. The channel
        // type alone does not identify the pixel type, so also compare the
        // total band count (channels x planes covers either file layout).
        if (expected_channel_type != VW_CHANNEL_UNKNOWN ||
            expected_num_channels > 0) {
          boost::shared_ptr<DiskImageResource>
            rsrc(DiskImageResource::open(output_file));
          if (expected_channel_type != VW_CHANNEL_UNKNOWN &&
              rsrc->channel_type() != expected_channel_type)
            return true;
          if (expected_num_channels > 0 &&
              rsrc->channels() * rsrc->planes() != expected_num_channels)
            return true;
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
                        float valid_max = std::numeric_limits<float>::quiet_NaN()) {

    double big = std::numeric_limits<double>::max();
    bool has_valid_range = !std::isnan(valid_min) && !std::isnan(valid_max);

    boost::shared_ptr<vw::DiskImageResource> rsrc(vw::DiskImageResourcePtr(file));
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

    cartography::GeoReference georef;
    bool has_georef = cartography::read_georeference(georef, base_file);

    // Keep making more pyramid levels until they are small enough
    int level = 0;
    int scale = 1;
    while (double(m_pyramid[level].cols()) * double(m_pyramid[level].rows())
           > m_lowest_resolution_subimage_num_pixels) {

      // The name of the file at the current scale
      std::ostringstream os;
      scale *= subsample;
      os <<  "_sub" << scale << ".tif";
      std::string suffix = os.str();

      // If a valid pixel range was provided and the pixel type is scalar,
      // use range masking. Otherwise fall back to exact nodata comparison.
      ImageViewRef<PixelMask<PixelT>> masked
        = maskForPyramid(m_pyramid[level], m_nodata_val,
                         m_valid_min, m_valid_max);
      double sub_scale   = 1.0/subsample;
      int    tile_size   = 256;
      int    sub_threads = 1;

      // Resample the image at the current pyramid level.
      // TODO: resample_aa is a hacky thingy. Need to understand
      // what is a good way of resampling.
      // Note that below we cast the channels to double for resampling,
      // then cast back to current pixel type for saving.
      // Fill value for pixels invalidated by the mask. When an integer
      // scalar PixelT cannot represent m_nodata_val (NaN or out of range),
      // converting it is undefined behavior; use a value-initialized pixel
      // instead (it is never consumed, as maskForPyramid marks all pixels
      // valid then). Floating-point and compound PixelT keep the existing
      // conversion.
      PixelT nodata_pixel = PixelT();
      if (!unrepresentableNodata<PixelT>(m_nodata_val))
        set_all(nodata_pixel, m_nodata_val);
      ImageViewRef<PixelT> unmasked
        = block_rasterize(cache_tile_aware_render
                          (pixel_cast<PixelT>
                           (apply_mask
                            (resample_aa
                             (channel_cast<double>(masked), sub_scale),
                             nodata_pixel
                             )),
                           Vector2i(tile_size,tile_size) * sub_scale
                           ), Vector2i(tile_size,tile_size), sub_threads);
      
      // Write the current image.
      if (has_georef)
        georef = resample(georef, sub_scale);

      // If the file exists, and has the right size, and is not too old,
      // don't write it again
      std::string curr_file = filename_from_suffix1(base_file, suffix);
      ChannelTypeEnum pixel_channel_type
        = ChannelTypeID<typename CompoundChannelType<PixelT>::type>::value;
      int pixel_num_channels = CompoundNumChannels<PixelT>::value;
      bool will_write = overwrite_if_no_good(base_file, curr_file,
                                             unmasked.cols(), unmasked.rows(),
                                             pixel_channel_type,
                                             pixel_num_channels);

      if (will_write) {

        if (level == 0)
          vw_out() << "Construct an image pyramid for: " << base_file << "\n";

        TerminalProgressCallback tpc("vw", ": ");
        vw_out() << "Writing: " << curr_file << "\n";
        try{
          cartography::block_write_gdal_image(curr_file, unmasked, has_georef, georef,
                      has_nodata, m_nodata_val, opt, tpc);
        } catch(...) {
          vw_out() << "Failed to write: " << curr_file << "\n";
          curr_file = filename_from_suffix2(base_file, suffix);
          will_write = overwrite_if_no_good(base_file, curr_file,
                    unmasked.cols(), unmasked.rows(), pixel_channel_type,
                    pixel_num_channels);
          if (will_write) {
            vw_out() << "Writing: " << curr_file << "\n";
            cartography::block_write_gdal_image(curr_file, unmasked, has_georef, georef,
                  has_nodata, m_nodata_val, opt, tpc);
          }
        }
      }

      if (!will_write)
        vw_out() << "Using existing subsampled image: " << curr_file << "\n";

      // Note that m_pyramid contains a handle to DiskImageView.
      // DiskImageView's implementation will make it possible to
      // cache in memory the most recently used tiles of all
      // the images in the pyramid.
      m_pyramid_files.push_back(curr_file);
    
      m_temporary_files.insert(curr_file);
      m_pyramid.push_back(DiskImageView<PixelT>(curr_file));
      m_scales.push_back(scale);

      level++;
    } // End level creation loop

    // This is expensive, so cache it going forward
    m_approx_bounds = approx_bounds_nocache(m_pyramid_files.back(),
                                            has_nodata, m_nodata_val,
                                            m_valid_min, m_valid_max);
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

    region_in.crop(bounding_box(m_pyramid[0]));
    scale_out  = m_scales[level];
    region_out = region_in/scale_out;
    region_out.crop(bounding_box(m_pyramid[level]));

    clip = crop(m_pyramid[level], region_out);
  }

}} // End namespace vw::mosaic

#endif // __VW_MOSAIC_DISKIMAGEPYRAMID_H__
