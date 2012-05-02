// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Misc.h

#ifndef __VW_PHOTOMETRY_MISC_H__
#define __VW_PHOTOMETRY_MISC_H__

#include <sys/stat.h>
#include <math.h>
#include <time.h>
#include <vector>

#include <vw/Image/ImageViewBase.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Photometry/Reconstruct.h>

namespace vw {
  namespace photometry {

    void upsample_uint8_image(std::string output_file, std::string input_file, int upsampleFactor);
    
  //upsamples a DEM by an upsampleFactor
  void upsample_image(std::string output_file, std::string input_file, int upsampleFactor);

  //subsamples a geo referenced tiff image by two
  void subsample_image(std::string output_file, std::string input_file);

  // Given two images and two georeferences, this function picks a set
  // of matching pixel samples between the two images.  It rejects
  // pixels that are not valid, and it should probably also reject
  // pixels that are near saturation (though it does not yet!).
  template <class ViewT>
  std::vector<Vector4> sample_images(ImageViewBase<ViewT> const& image1,
                                     ImageViewBase<ViewT> const& image2,
                                     cartography::GeoReference const& geo1,
                                     cartography::GeoReference const& geo2,
                                     int num_samples,
                                     std::string const& DEM_file,
                                     std::vector<Vector3> *normalArray,
                                     std::vector<Vector3> *xyzArray );


  /// Erases a file suffix if one exists and returns the base string
  static std::string prefix_from_filename(std::string const& filename) {
    std::string result = filename;
    int index = result.rfind(".");
    if (index != -1)
      result.erase(index, result.size());
    return result;
  }

  /// Erases a file suffix if one exists and returns the base string less3 characters
  static std::string prefix_less3_from_filename(std::string const& filename) {
    std::string result = filename;
    int index = result.rfind(".");
    if (index != -1)
      result.erase(index-3, result.size()+3);
    return result;
  }

  /// Erases a file suffix if one exists and returns the base string less3 characters
  static std::string sufix_from_filename(std::string const& filename) {
    std::string result = filename;
    int index = result.rfind("/");
    if (index != -1)
      result.erase(0, index);
    return result;
  }

  //reads the tiff DEM into a 3D coordinate
  //pos is a Vector2 of pixel coordinates, GR is georeference
  template <class ViewT>
  Vector3 pixel_to_cart (Vector2 pos, ImageViewBase<ViewT> const& img,
                         cartography::GeoReference GR);


  std::vector<std::string> parse_command_arguments(int argc, char *argv[] );

  void getTileCornersWithoutPadding(// Inputs
                                    int numCols, int numRows,
                                    cartography::GeoReference const& geoRef,
                                    double tileSize, int pixelPadding,
                                    // Outputs
                                    double & min_x, double & max_x,
                                    double & min_y, double & max_y
                                    );
  
    void applyPaddingToTileCorners(// Inputs
                                   cartography::GeoReference const& geoRef,
                                   int pixelPadding,
                                   double min_x, double max_x,
                                   double min_y, double max_y,
                                   // Outputs
                                   double & min_x_padded, double & max_x_padded,
                                   double & min_y_padded, double & max_y_padded);
    
    void readDEMTilesIntersectingBox(// Inputs
                                     double noDEMDataValue,
                                     Vector2 boxNW, Vector2 boxSE,
                                     std::vector<std::string> const& DEMTiles,
                                     // Outputs
                                     ImageView<PixelGray<float> > & combinedDEM,
                                     cartography::GeoReference    & combinedDEM_geo);

    void listTifsInDir(const std::string & dirName,
                       std::vector<std::string> & tifsInDir
                       );

    void writeSunAndSpacecraftPosition(std::string prefix, std::string sunFile, std::string spacecraftFile,
                                       Vector3 sunPosition, Vector3 spacecraftPosition);
    
    std::string getFirstElevenCharsFromFileName(std::string fileName);
    
    void indexFilesByKey(std::string dirName, std::map<std::string, std::string> & index);

    void enforceUint8Img(std::string imgName);

    bool readNoDEMDataVal(std::string DEMFile, float & noDEMDataValue);
    
    void maskPixels(std::string imgFile, std::string maskFile, double shadowThresh, std::string outDir);

  template <class pixelInType, class pixelOutType>
  bool getSubImageWithMargin(// Inputs
                             Vector2 begLonLat, Vector2 endLonLat,
                             std::string imageFile,
                             // Outputs
                             ImageView<pixelOutType>   & subImage,
                             cartography::GeoReference & sub_geo){
    
    // Read from disk only the portion of image whose pixel lon lat
    // values are in between begLonLat and endLonLat.  Adjust the
    // georeference accordingly. This is necessary to save on memory
    // for large images.
    
    // IMPORTANT: We assume that the input image is of pixelInType, and
    // we will convert the portion we read into pixelOutType.  Also we
    // read a few more pixels on each side, to help later with
    // interpolation.
  
    int extra = 2;
    DiskImageView<pixelInType> input_img(imageFile);
    //std::cout << "Reading: " << imageFile << std::endl;
    
    cartography::GeoReference input_geo;
    read_georeference(input_geo, imageFile);

    Vector2 begPix = input_geo.lonlat_to_pixel(begLonLat);
    Vector2 endPix = input_geo.lonlat_to_pixel(endLonLat);

    int beg_col = std::max(0,                (int)floor(begPix(0)) - extra);
    int end_col = std::min(input_img.cols(), (int)ceil(endPix(0))  + extra);
    int beg_row = std::max(0,                (int)floor(begPix(1)) - extra);
    int end_row = std::min(input_img.rows(), (int)ceil(endPix(1))  + extra);
  
    if (beg_col >= end_col || beg_row >= end_row) return false;
  
    subImage.set_size(end_col - beg_col, end_row - beg_row);
    for (int row = beg_row; row < end_row; row++){
      for (int col = beg_col; col < end_col; col++){
        // Note the pixelInType to pixelOutType conversion below
        subImage(col - beg_col, row - beg_row) = input_img(col, row); 
      }
    }
  
    // In sub_geo, the (0, 0) pixel will where (beg_col, beg_row) is in input_geo.
    sub_geo = vw::cartography::crop(input_geo, beg_col, beg_row);

    //system("echo getSubImageWithMargin top is $(top -u $(whoami) -b -n 1|grep lt-reconstruct)");
  
    return true;
  }

  template<class T>
  void dumpImageToFile(T & img, std::string fileName){
    
    printf("dumping %s\n", fileName.c_str());
    
    std::ofstream fs(fileName.c_str());
    fs.precision(20);
    
    for (int k = 0 ; k < img.rows(); ++k) {
      for (int l = 0; l < img.cols(); ++l) {
        fs << (double)img(l, k) << " ";
      }
      fs << std::endl;
    }
    fs.close();
  }

}} // end vw::photometry

namespace vw {
namespace photometry {

  template <class ImageT>
  bool isGoodPixel(ImageT const& maskImg, cartography::GeoReference maskGeo, 
                   Vector2 lonlat, double t){
  
    Vector2 maskPix = maskGeo.lonlat_to_pixel(lonlat);
    float x = maskPix(0), y = maskPix(1);
  
    int x0 = (int)floor(x), x1 = (int)ceil(x);
    int y0 = (int)floor(y), y1 = (int)ceil(y);
    bool isGood = ((x >= 0) && (x <= maskImg.cols()-1) &&
                   (y >= 0) && (y <= maskImg.rows()-1) && 
                   is_valid(maskImg(x0, y0)) && (float)maskImg(x0, y0) >= t &&
                   is_valid(maskImg(x0, y1)) && (float)maskImg(x0, y1) >= t &&
                   is_valid(maskImg(x1, y0)) && (float)maskImg(x1, y0) >= t &&
                   is_valid(maskImg(x1, y1)) && (float)maskImg(x1, y1) >= t
                   );
  
    return isGood;
  }

  // Mask the pixels in a given image using pixels from a mask image.
  // All the process is lazy, the full images are never stored in memory.

  template <class ImageT>
  class MaskImage : public ImageViewBase<MaskImage<ImageT> > {

    typedef typename boost::mpl::if_<IsFloatingPointIndexable<ImageT>, double, int32>::type offset_type;

    ImageT m_inputImg;
    const ImageViewRef<PixelMask<PixelGray<uint8> > > & m_maskImg;
    double m_shadowThresh;
    cartography::GeoReference m_imgGeo;
    cartography::GeoReference m_maskGeo;
    
  public:
    typedef typename ImageT::pixel_type pixel_type;
    typedef const pixel_type result_type;
    typedef ProceduralPixelAccessor<MaskImage> pixel_accessor;

    MaskImage(ImageViewRef<PixelMask<PixelGray<uint8> > > const& inputImg,
              ImageT const& inputImg_prerasterized,
              ImageViewRef<PixelMask<PixelGray<uint8> > > const& maskImg,
              double shadowThresh, 
              cartography::GeoReference const& imgGeo,
              cartography::GeoReference const& maskGeo):
      m_inputImg(inputImg_prerasterized),
      m_maskImg(maskImg),
      m_shadowThresh(shadowThresh),
      m_imgGeo(imgGeo),
      m_maskGeo(maskGeo){}
    
    inline int32 cols() const { return m_inputImg.cols(); }
    inline int32 rows() const { return m_inputImg.rows(); }
    inline int32 planes() const { return m_inputImg.planes(); }

    inline pixel_accessor origin() const { return pixel_accessor(*this); }

    inline result_type operator()( offset_type i, offset_type j, int32 p=0 ) const {

      if  ((float)m_inputImg(i, j) == 0) result_type(0.0);
      
      // Note: Below we take into account that the lonlat of the image
      // may be shifted by 360 degrees from the lonlat of the mask.
      // This is a bug fix.
      Vector2 lonlat = m_imgGeo.pixel_to_lonlat(Vector2(i, j));
      double t = m_shadowThresh;
      if (!isGoodPixel(m_maskImg, m_maskGeo, lonlat, t)                   &&
          !isGoodPixel(m_maskImg, m_maskGeo, lonlat + Vector2(360, 0), t) && 
          !isGoodPixel(m_maskImg, m_maskGeo, lonlat - Vector2(360, 0), t)
          ){
        return result_type(0.0);
      }else{
        return m_inputImg(i, j);
      }
      
    }
    
    /// \cond INTERNAL
    typedef MaskImage<typename ImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
      return prerasterize_type(m_inputImg, m_inputImg.prerasterize(bbox), m_maskImg, m_shadowThresh, 
                               m_imgGeo, m_maskGeo);
    }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };

  // --------------------------------------------------------------------------
  // Functional API. See the description of the MaskImage for more info.
  // --------------------------------------------------------------------------
  template <class ImageT>
  MaskImage<ImageT>
  mask_image(ImageViewBase<ImageT> const& inputImg,
             ImageViewBase<ImageT> const& maskImg,
             double shadowThresh,
             cartography::GeoReference const& imgGeo,
             cartography::GeoReference const& maskGeo) {
    return MaskImage<ImageT>(inputImg.impl(), inputImg.impl(), maskImg.impl(), shadowThresh, imgGeo, maskGeo);
  }

} // namespace photometry

  /// \cond INTERNAL
  // Type traits
  template <class ImageT>
  struct IsFloatingPointIndexable< photometry::MaskImage<ImageT> > : public IsFloatingPointIndexable<ImageT> {};
  /// \endcond
}

#endif//__VW_PHOTOMETRY_MISC_H__
