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
#include <vector>
#include <list>

#include <vw/Core/Cache.h>
#include <vw/Core/ProgressCallback.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Cartography/GeoReferenceUtils.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/FileUtils.h>
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

  // TODO: Clean up!
  /// Single-channel images are read into Image<double>, while
  /// multi-channel in Image<Vector>, and we skip reading extra channels.
  template<class PixelT>
  typename boost::enable_if<boost::is_same<PixelT,double>, ImageViewRef<PixelT> >::type
  custom_read(std::string const& file){
    return DiskImageView<PixelT>(file);
  }
  template<class PixelT>
  typename boost::disable_if<boost::is_same<PixelT,double>, ImageViewRef<PixelT> >::type
  custom_read(std::string const& file){
    return vw::read_channels<vw::math::VectorSize<PixelT>::value, typename PixelT::value_type>(file, 0);
  }


  // TODO: Clean up!
  // Gets called for PixelT == double
  template<class PixelT>
  typename boost::enable_if< boost::mpl::or_< boost::is_same<PixelT,double>, 
                             boost::is_same<PixelT,vw::uint8> >, ImageViewRef< PixelMask<PixelT> > >::type
  create_custom_mask(ImageViewRef<PixelT> & img, double nodata_val){
    return create_mask_less_or_equal(img, nodata_val);
  }
  // Gets called for PixelT == Vector<u8, 1> and Vector<u8, 3>
  template<class PixelT>
  typename boost::disable_if< boost::mpl::or_< boost::is_same<PixelT,double>, 
                              boost::is_same<PixelT,vw::uint8> >, ImageViewRef< PixelMask<PixelT> > >::type
  create_custom_mask(ImageViewRef<PixelT> & img, double nodata_val){
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

  /// Logic to find some approximate values for the valid pixels, ignoring the worst
  /// outliers. Use the lowest pyramid level.
  template <class PixelT>
  typename boost::enable_if<boost::is_same<PixelT,double>, vw::Vector2 >::type
  get_approx_bounds_noclass(std::string const& file){
    
    double big = std::numeric_limits<double>::max();
    
    // Get the nodata value, if present
    double nodata_val = -big;
    bool has_nodata = read_nodata_val(file, nodata_val);
    
    // To make life easy
    if (std::isnan(nodata_val)) nodata_val = -big;
    
    std::vector<double> vals;
    DiskImageView<PixelT> img(file);
    for (int col = 0; col < img.cols(); col++) {
      for (int row = 0; row < img.rows(); row++) {
	
	if (std::isnan(img(col, row))) continue;
	if (has_nodata && img(col, row) == nodata_val) continue;
	vals.push_back(img(col, row));
      }
    }
    
    vw::Vector2 bounds(-big, big);
    if (vals.empty()) return bounds;
    
    //  Find sane bounds, even if they are somewhat exaggerated 
    double pct = 0.25;
    double outlier_factor = 4;
    double b, e;
    vw::math::find_outlier_brackets(vals, pct, outlier_factor, b, e);

    return vw::Vector2(b, e);
  }
  
  template <class PixelT>
  typename boost::disable_if<boost::is_same<PixelT,double>, vw::Vector2 >::type
  get_approx_bounds_noclass(std::string const& file){
    return vw::Vector2(); // multi-channel image
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
    DiskImagePyramid(std::string const& base_file = "",
		     cartography::GdalWriteOptions const& opt = cartography::GdalWriteOptions(),
		     int top_image_max_pix = 1000*1000,
		     int subsample = 2);

    // Given a region (at full resolution) and a scale factor, compute
    // the portion of the image in the region, subsampled by a factor no
    // more than the input scale factor.  Also return the precise subsample
    // factor used and the region at that scale level.
    void get_image_clip(double scale_in, BBox2i region_in,
			ImageView<PixelT> & clip, double & scale_out, BBox2i & region_out) const;

    vw::Vector2 get_approx_bounds() const {
      return get_approx_bounds_noclass<PixelT>(m_pyramid_files.back());
    }
    
    ~DiskImagePyramid() {}

    // These all describe the highest resolution pyramid layer
    int32 cols  () const { return m_pyramid[0].cols(); }
    int32 rows  () const { return m_pyramid[0].rows(); }
    int32 planes() const { return m_pyramid[0].planes(); }

    double get_nodata_val() const {return(m_nodata_val);}

    /// Return the highest resolution pyramid layer
    ImageViewRef<PixelT>        bottom()       { return m_pyramid[0]; }
    ImageViewRef<PixelT> const& bottom() const { return m_pyramid[0]; }

    std::set<std::string> const& get_temporary_files() const {return m_temporary_files;}

  private:

    cartography::GdalWriteOptions m_opt;

    // The subsample factor to go to the next level of the pyramid (must be >= 2).
    int m_subsample;

    // The maxiumum number of pixels in the coarsest level of the pyramid
    // (keep on downsampling until getting to this number or under it).
    int m_top_image_max_pix;

    //  The pyramid. Largest images come earlier.
    std::vector< ImageViewRef<PixelT> > m_pyramid;

    // The files (stored on disk) containing the images in the pyramid.
    std::vector<std::string> m_pyramid_files;

    // We may wipe these at the end
    std::vector<std::string> m_cached_files;

    double m_nodata_val;
    std::vector<int> m_scales;
    
    /// Contains all the temporary image files we create
    std::set<std::string> m_temporary_files; 
  };


  //#################################################################################
  // Function definitions

  template <class PixelT>
  DiskImagePyramid<PixelT>::DiskImagePyramid(std::string const& base_file,
                                             cartography::GdalWriteOptions const& opt,
                                             int top_image_max_pix,
                                             int subsample):
    m_opt(opt), m_subsample(subsample),
    m_top_image_max_pix(top_image_max_pix),
    m_nodata_val(std::numeric_limits<double>::quiet_NaN()) {
    
    if (base_file.empty())
      return;

    if (subsample < 2) {
      vw_throw( ArgumentErr() << "Must subsample by a factor of at least 2.\n");
    }

    if (top_image_max_pix < 4) {
      vw_throw( ArgumentErr() << "The image at the top of the pyramid must be at least 2x2 in size.\n");
    }

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
    while (double(m_pyramid[level].cols())*double(m_pyramid[level].rows())
           > m_top_image_max_pix ){

      // The name of the file at the current scale
      std::ostringstream os;
      scale *= subsample;
      os <<  "_sub" << scale << ".tif";
      std::string suffix = os.str();

      if (level == 0) {
        vw_out() << "Detected large image: " << base_file  << "." << std::endl;
        vw_out() << "Will construct an image pyramid on disk."  << std::endl;
      }

      ImageViewRef< PixelMask<PixelT> > masked
        = create_custom_mask(m_pyramid[level], m_nodata_val);
      double sub_scale   = 1.0/subsample;
      int    tile_size   = 256;
      int    sub_threads = 1;

      // Resample the image at the current pyramid level.
      // TODO: resample_aa is a hacky thingy. Need to understand
      // what is a good way of resampling.
      // Note that below we cast the channels to double for resampling,
      // then cast back to current pixel type for saving.
      PixelT nodata_pixel;
      set_all(nodata_pixel, m_nodata_val);
      ImageViewRef<PixelT> unmasked
        = block_rasterize
        (cache_tile_aware_render
          (pixel_cast<PixelT>
          (apply_mask
            (resample_aa
            (channel_cast<double>(masked), sub_scale),
            nodata_pixel
            )),
          Vector2i(tile_size,tile_size) * sub_scale
          ), Vector2i(tile_size,tile_size), sub_threads
        );
      
      // Write the current image.
      if (has_georef)
        georef = resample(georef, sub_scale);

      // If the file exists, and has the right size, and is not too old,
      // don't write it again
      std::string curr_file = filename_from_suffix1(base_file, suffix);
      bool will_write = overwrite_if_no_good(base_file, curr_file,
                                             unmasked.cols(), unmasked.rows());

      if (will_write) {
        TerminalProgressCallback tpc("vw", ": ");
        vw_out() << "Writing: " << curr_file << std::endl;
        try{
          cartography::block_write_gdal_image(curr_file, unmasked, has_georef, georef,
                      has_nodata, m_nodata_val, opt, tpc);
        }catch(...){
          vw_out() << "Failed to write: " << curr_file << "\n";
          curr_file = filename_from_suffix2(base_file, suffix);
          will_write = overwrite_if_no_good(base_file, curr_file,
                    unmasked.cols(), unmasked.rows());
          if (will_write) {
            vw_out() << "Writing: " << curr_file << std::endl;
            cartography::block_write_gdal_image(curr_file, unmasked, has_georef, georef,
                  has_nodata, m_nodata_val, opt, tpc);
          }
        }
      }


      if (!will_write)
        vw_out() << "Using existing subsampled image: " << curr_file << std::endl;

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

  }

  template <class PixelT>
  void DiskImagePyramid<PixelT>::get_image_clip(double scale_in, BBox2i region_in,
                                                ImageView<PixelT> & clip,
                                                double & scale_out, BBox2i & region_out) const {

    if (m_pyramid.empty())
      vw_throw( ArgumentErr() << "Uninitialized image pyramid.\n");

    // Find the right pyramid level to use
    int level = 0;
    while (1) {
      if (level+1 >= (int)m_scales.size()) break; // last level
      if (m_scales[level+1] > scale_in)  break; // too coarse
      level++;
    }

    //vw_out() << "Reading: " << m_pyramid_files[level] << std::endl;

    region_in.crop(bounding_box(m_pyramid[0]));
    scale_out  = m_scales[level];
    region_out = region_in/scale_out;
    region_out.crop(bounding_box(m_pyramid[level]));

    clip = crop(m_pyramid[level], region_out);
  }

}} // End namespace vw::mosaic

#endif // __VW_MOSAIC_DISKIMAGEPYRAMID_H__
