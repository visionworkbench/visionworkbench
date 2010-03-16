// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <string>
#include <fstream>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>

#include <boost/operators.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/fstream.hpp>
namespace fs = boost::filesystem;

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
using namespace vw;
using namespace vw::math;
using namespace vw::cartography;

#include <math.h>
#include <vw/Photometry/Reconstruct.h>


void ComputeSaveShadowMap(std::string origfile, std::string shadowMapFile, GlobalParams globalParams)
{
  DiskImageView<PixelMask<PixelGray<uint8> > > originalImage(origfile);
  ImageView<PixelMask<PixelGray<uint8> > > shadowImage(originalImage.cols(), originalImage.rows());


  GeoReference originalGeo;
  read_georeference(originalGeo, origfile);

  for (unsigned k=0; k < originalImage.rows(); ++k) {
      for (unsigned l=0; l < originalImage.cols(); ++l) {


         Vector2 sample_pix(l,k);

         if ( is_valid(originalImage(l,k)) ) {
           shadowImage(l, k) = 0;
           if (originalImage(l, k) < globalParams.shadowThresh){
               shadowImage(l, k) = 255;
           }
         }

      }
    }

    write_georeferenced_image(shadowMapFile,
                              channel_cast<uint8>(shadowImage),
                              originalGeo, TerminalProgressCallback("{Core}","Processing:"));
}
