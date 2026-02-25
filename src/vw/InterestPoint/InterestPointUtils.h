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


/// \file InterestPointUtils.h
///
/// Keep here some low-level utilities to avoid cluttering InterestData.h
/// which is included in many places.
///

// TODO(oalexan1): Once many ip finding functions are no longer templates,
// do the same to these ones.

#ifndef __INTEREST_POINT_UTILS_H__
#define __INTEREST_POINT_UTILS_H__

#include <vw/InterestPoint/InterestPoint.h>

#include <vw/Image/Algorithms.h>
#include <vw/Image/MaskViews.h>
#include <vw/Image/ImageMath.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Math/RandomSet.h>

namespace vw {
namespace ip {

/// Remove all pixels within a radius of the border or nodata.
// Note: A nodata pixel is one for which pixel <= nodata.
template <class ImageT>
void remove_ip_near_nodata(vw::ImageViewBase<ImageT> const& image,   double nodata,
                           vw::ip::InterestPointList      & ip_list, int    radius ){

  size_t prior_ip = ip_list.size();

  typedef ImageView<typename ImageT::pixel_type> CropImageT;
  const int width = 2*radius+1;
  CropImageT subsection(width,width);   

  // Get shrunk bounding box
  BBox2i bound = bounding_box( image.impl() );
  bound.contract(radius);

  // Loop through all the points
  for ( ip::InterestPointList::iterator ip = ip_list.begin(); ip != ip_list.end(); ++ip ) {

    // Remove the point if it was too close to the image borders
    if ( !bound.contains( Vector2i(ip->ix,ip->iy) ) ) {
      ip = ip_list.erase(ip);
      ip--;
      continue;
    }

    // Remove the point if any invalid pixels are too close to it
    subsection = crop( image.impl(), ip->ix-radius, ip->iy-radius, width, width );
    for ( typename CropImageT::iterator pixel = subsection.begin();
          pixel != subsection.end(); pixel++ ) {
      if (*pixel == nodata) {
        ip = ip_list.erase(ip);
        ip--;
        break;
      }
    }
  }
  VW_OUT( DebugMessage, "asp" ) << "Removed " << prior_ip - ip_list.size()
                                << " interest points due to their proximity to nodata values."
                                << std::endl << "Nodata value used " << nodata << std::endl;
} // End function remove_ip_near_nodata

/// Draw a helpful image showing where IPs are detected in an image.
/// - If reduce_scale is set the IP circle markers will be smaller.
template <typename ImageT>
void write_ip_debug_image(std::string const& out_file_name,
                          vw::ImageViewBase<ImageT> const& image,
                          InterestPointList const& ip,
                          bool has_nodata, double nodata,
                          bool force_full_res=false,
                          bool reduce_scale  =false) {

  // Scale the images to keep the size down below 1500x1500 so they draw quickly.
  float sub_scale  = sqrt(1500.0 * 1500.0 / float(image.impl().cols() * image.impl().rows()));
	sub_scale /= 2;
  if ((sub_scale > 1) || force_full_res)
    sub_scale = 1;

  vw_out() << "Writing debug image: " << out_file_name << " with downsample: "<< sub_scale << "\n";
  
  // Get the min and max point intensity
  vw_out(InfoMessage,"interest_point") << "\t > Gathering statistics:\n";
  float min = 1e30, max = -1e30;
  for ( InterestPointList::const_iterator point = ip.begin(); point != ip.end(); ++point ) {
    if ( point->interest > max )
      max = point->interest;
    if ( point->interest < min )
      min = point->interest;
  }
  float diff = max - min;

  // Create a muted version of the input image
  ImageView<PixelRGB<uint8> > oimage;
  const double IMAGE_FADE_PERCENTAGE = 0.7;
  if ( has_nodata ) {
    ImageView<PixelMask<float> > gray = resample(create_mask(image, nodata), sub_scale);
    oimage = pixel_cast_rescale<PixelRGB<uint8> >(apply_mask(normalize(gray) * IMAGE_FADE_PERCENTAGE));
  } else {
    oimage = pixel_cast_rescale<PixelRGB<uint8> >(normalize(resample(image, sub_scale)) * IMAGE_FADE_PERCENTAGE);
  }

  // Draw each of the interest points on to the muted input image.
  for ( InterestPointList::const_iterator point = ip.begin(); point != ip.end(); ++point ) {
    float norm_i = (point->interest - min)/diff;
    PixelRGB<uint8> color(0,0,0);
    if ( norm_i < .5 ) {
      // Form of red
      color.r() = 255;
      color.g() = (uint8)(2*norm_i*255);
    } else {
      // Form of green
      color.g() = 255;
      color.r() = 255 - (uint8)(2*(norm_i-.5)*255);
    }

    // Marking point w/ Dot
    int ix = point->ix*sub_scale;
    int iy = point->iy*sub_scale;
    int x  = point->x*sub_scale;
    int y  = point->y*sub_scale;
    oimage(ix,iy) = color;
    
    double scale = point->scale;
    if (reduce_scale)
      scale = 0.12*scale;

    // Circling point based on the scale
    for (float a = 0; a < 6; a+=.392 ) {
      float a_d = a + .392; // ?
      Vector2i start( int(2*scale*cos(a  )+x),
                      int(2*scale*sin(a  )+y) );
      Vector2i end(   int(2*scale*cos(a_d)+x),
                      int(2*scale*sin(a_d)+y) );
      draw_line( oimage, color, start, end );
    }
  } // End loop through interest points

  boost::scoped_ptr<DiskImageResource> rsrc( DiskImageResource::create(out_file_name,
                                                                       oimage.format() ) );
  vw_out(InfoMessage,"interest_point") << "\t > Writing out ip debug image: "
                                       << out_file_name << "\n";
  write_image( *rsrc, oimage);
}

template <class T>
void pick_pair_subset(std::vector<T> & ip1,
                      std::vector<T> & ip2,
                      int max_pairwise_matches) {

  if (max_pairwise_matches >= 0 && (int)ip1.size() > max_pairwise_matches) {

    std::vector<int> subset;
    vw::math::pick_random_indices_in_range(ip1.size(), max_pairwise_matches, subset);
    std::sort(subset.begin(), subset.end()); // sort the indices; not strictly necessary
    
    std::vector<T> ip1_full, ip2_full;
    ip1_full.swap(ip1);
    ip2_full.swap(ip2);
    
    ip1.resize(max_pairwise_matches);
    ip2.resize(max_pairwise_matches);
    for (size_t it = 0; it < subset.size(); it++) {
      ip1[it] = ip1_full[subset[it]];
      ip2[it] = ip2_full[subset[it]];
    }
  }
  
  return;
}

}} // namespace vw::ip

#endif //__INTEREST_POINT_UTILS_H__
