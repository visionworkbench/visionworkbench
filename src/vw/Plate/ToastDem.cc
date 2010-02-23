// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/ToastDem.h>
#include <vw/Image/Interpolation.h>
using namespace vw;

// ---------------------------------------------------------------------------------
//                                WWT DEM Indices
// ---------------------------------------------------------------------------------

namespace {

  // DEM triangle indices for the upper right (0-90 degrees) and lower
  // left (180-270 degrees) quadrants of the TOAST space.
  static int const ur_ll_quadrant_u_toast_indices[] = { 0,16,32,0,16,32,0,16,32,8,8,16,8,0,8,12,12,16,12,16,12,8,4,4,12,12,8,4,0,4,8,12,12,4,4,0,4,8,4,14,14,16,14,16,14,12,10,10,14,14,12,14,16,14,12,10,10,14,14,16,14,12,14,4,2,2,6,6,8,6,8,6,4,6,6,10,10,8,10,8,10,12,14,14,10,10,12,2,0,2,4,6,6,2,2,0,2,4,2,12,14,14,10,10,8,10,8,10,12,10,10,2,2,0,2,0,2,4,6,6,2,2,4,6,8,6,4,2,2,6,6,8,6,4,6,24,24,32,24,16,24,28,28,32,28,32,28,24,20,20,28,28,24,20,16,20,24,28,28,20,20,16,20,24,20,30,30,32,30,32,30,28,26,26,30,30,28,30,32,30,28,26,26,30,30,32,30,28,30,20,18,18,22,22,24,22,24,22,20,22,22,26,26,24,26,24,26,28,30,30,26,26,28,18,16,18,20,22,22,18,18,16,18,20,18,28,30,30,26,26,24,26,24,26,28,26,26,18,18,16,18,16,18,20,22,22,18,18,20,22,24,22,20,18,18,22,22,24,22,20,22,8,8,16,8,0,8,12,12,16,12,16,12,8,4,4,12,12,8,4,0,4,8,12,12,4,4,0,4,8,4,14,14,16,14,16,14,12,10,10,14,14,12,14,16,14,12,10,10,14,14,16,14,12,14,4,2,2,6,6,8,6,8,6,4,6,6,10,10,8,10,8,10,12,14,14,10,10,12,2,0,2,4,6,6,2,2,0,2,4,2,12,14,14,10,10,8,10,8,10,12,10,10,2,2,0,2,0,2,4,6,6,2,2,4,6,8,6,4,2,2,6,6,8,6,4,6,24,24,32,24,16,24,28,28,32,28,32,28,24,20,20,28,28,24,20,16,20,24,28,28,20,20,16,20,24,20,30,30,32,30,32,30,28,26,26,30,30,28,30,32,30,28,26,26,30,30,32,30,28,30,20,18,18,22,22,24,22,24,22,20,22,22,26,26,24,26,24,26,28,30,30,26,26,28,18,16,18,20,22,22,18,18,16,18,20,18,28,30,30,26,26,24,26,24,26,28,26,26,18,18,16,18,16,18,20,22,22,18,18,20,22,24,22,20,18,18,22,22,24,22,20,22};

  static int const ur_ll_quadrant_v_toast_indices[] = {0,0,0,16,16,16,32,32,32,0,8,8,8,8,16,8,12,12,4,4,0,4,0,4,8,4,4,12,12,16,12,16,12,8,4,4,12,12,8,12,14,14,10,10,8,10,8,10,12,10,10,2,2,0,2,0,2,4,6,6,2,2,4,2,0,2,4,6,6,2,2,0,2,4,2,4,2,2,6,6,8,6,8,6,4,6,6,14,14,16,14,16,14,12,10,10,14,14,12,14,16,14,12,10,10,14,14,16,14,12,14,4,2,2,6,6,8,6,8,6,4,6,6,10,10,8,10,8,10,12,14,14,10,10,12,0,8,8,8,8,16,8,12,12,4,4,0,4,0,4,8,4,4,12,12,16,12,16,12,8,4,4,12,12,8,12,14,14,10,10,8,10,8,10,12,10,10,2,2,0,2,0,2,4,6,6,2,2,4,2,0,2,4,6,6,2,2,0,2,4,2,4,2,2,6,6,8,6,8,6,4,6,6,14,14,16,14,16,14,12,10,10,14,14,12,14,16,14,12,10,10,14,14,16,14,12,14,4,2,2,6,6,8,6,8,6,4,6,6,10,10,8,10,8,10,12,14,14,10,10,12,16,24,24,24,24,32,24,28,28,20,20,16,20,16,20,24,20,20,28,28,32,28,32,28,24,20,20,28,28,24,28,30,30,26,26,24,26,24,26,28,26,26,18,18,16,18,16,18,20,22,22,18,18,20,18,16,18,20,22,22,18,18,16,18,20,18,20,18,18,22,22,24,22,24,22,20,22,22,30,30,32,30,32,30,28,26,26,30,30,28,30,32,30,28,26,26,30,30,32,30,28,30,20,18,18,22,22,24,22,24,22,20,22,22,26,26,24,26,24,26,28,30,30,26,26,28,16,24,24,24,24,32,24,28,28,20,20,16,20,16,20,24,20,20,28,28,32,28,32,28,24,20,20,28,28,24,28,30,30,26,26,24,26,24,26,28,26,26,18,18,16,18,16,18,20,22,22,18,18,20,18,16,18,20,22,22,18,18,16,18,20,18,20,18,18,22,22,24,22,24,22,20,22,22,30,30,32,30,32,30,28,26,26,30,30,28,30,32,30,28,26,26,30,30,32,30,28,30,20,18,18,22,22,24,22,24,22,20,22,22,26,26,24,26,24,26,28,30,30,26,26,28};

  // DEM triangle indices for the upper left (90-180 degrees) and lower
  // left (270-360 degrees) quadrants of the TOAST space.
  static int const ul_lr_quadrant_u_toast_indices[] = {0,16,32,0,16,32,0,16,32,8,0,8,8,8,16,4,0,4,8,12,12,4,4,0,4,8,4,12,12,16,12,16,12,8,4,4,12,12,8,2,0,2,4,6,6,2,2,0,2,4,2,12,14,14,10,10,8,10,8,10,12,10,10,2,2,0,2,0,2,4,6,6,2,2,4,6,8,6,4,2,2,6,6,8,6,4,6,14,14,16,14,16,14,12,10,10,14,14,12,14,16,14,12,10,10,14,14,16,14,12,14,4,2,2,6,6,8,6,8,6,4,6,6,10,10,8,10,8,10,12,14,14,10,10,12,24,16,24,24,24,32,20,16,20,24,28,28,20,20,16,20,24,20,28,28,32,28,32,28,24,20,20,28,28,24,18,16,18,20,22,22,18,18,16,18,20,18,28,30,30,26,26,24,26,24,26,28,26,26,18,18,16,18,16,18,20,22,22,18,18,20,22,24,22,20,18,18,22,22,24,22,20,22,30,30,32,30,32,30,28,26,26,30,30,28,30,32,30,28,26,26,30,30,32,30,28,30,20,18,18,22,22,24,22,24,22,20,22,22,26,26,24,26,24,26,28,30,30,26,26,28,8,0,8,8,8,16,4,0,4,8,12,12,4,4,0,4,8,4,12,12,16,12,16,12,8,4,4,12,12,8,2,0,2,4,6,6,2,2,0,2,4,2,12,14,14,10,10,8,10,8,10,12,10,10,2,2,0,2,0,2,4,6,6,2,2,4,6,8,6,4,2,2,6,6,8,6,4,6,14,14,16,14,16,14,12,10,10,14,14,12,14,16,14,12,10,10,14,14,16,14,12,14,4,2,2,6,6,8,6,8,6,4,6,6,10,10,8,10,8,10,12,14,14,10,10,12,24,16,24,24,24,32,20,16,20,24,28,28,20,20,16,20,24,20,28,28,32,28,32,28,24,20,20,28,28,24,18,16,18,20,22,22,18,18,16,18,20,18,28,30,30,26,26,24,26,24,26,28,26,26,18,18,16,18,16,18,20,22,22,18,18,20,22,24,22,20,18,18,22,22,24,22,20,22,30,30,32,30,32,30,28,26,26,30,30,28,30,32,30,28,26,26,30,30,32,30,28,30,20,18,18,22,22,24,22,24,22,20,22,22,26,26,24,26,24,26,28,30,30,26,26,28};

  static int const ul_lr_quadrant_v_toast_indices[] = {0,0,0,16,16,16,32,32,32,0,8,8,8,16,8,8,12,12,4,4,0,4,0,4,8,4,4,12,16,12,8,4,4,12,12,16,12,8,12,12,14,14,10,10,8,10,8,10,12,10,10,2,2,0,2,0,2,4,6,6,2,2,4,2,0,2,4,6,6,2,2,0,2,4,2,4,2,2,6,6,8,6,8,6,4,6,6,14,16,14,12,10,10,14,14,16,14,12,14,4,2,2,6,6,8,6,8,6,4,6,6,14,14,16,14,16,14,12,10,10,14,14,12,10,8,10,12,14,14,10,10,8,10,12,10,0,8,8,8,16,8,8,12,12,4,4,0,4,0,4,8,4,4,12,16,12,8,4,4,12,12,16,12,8,12,12,14,14,10,10,8,10,8,10,12,10,10,2,2,0,2,0,2,4,6,6,2,2,4,2,0,2,4,6,6,2,2,0,2,4,2,4,2,2,6,6,8,6,8,6,4,6,6,14,16,14,12,10,10,14,14,16,14,12,14,4,2,2,6,6,8,6,8,6,4,6,6,14,14,16,14,16,14,12,10,10,14,14,12,10,8,10,12,14,14,10,10,8,10,12,10,16,24,24,24,32,24,24,28,28,20,20,16,20,16,20,24,20,20,28,32,28,24,20,20,28,28,32,28,24,28,28,30,30,26,26,24,26,24,26,28,26,26,18,18,16,18,16,18,20,22,22,18,18,20,18,16,18,20,22,22,18,18,16,18,20,18,20,18,18,22,22,24,22,24,22,20,22,22,30,32,30,28,26,26,30,30,32,30,28,30,20,18,18,22,22,24,22,24,22,20,22,22,30,30,32,30,32,30,28,26,26,30,30,28,26,24,26,28,30,30,26,26,24,26,28,26,16,24,24,24,32,24,24,28,28,20,20,16,20,16,20,24,20,20,28,32,28,24,20,20,28,28,32,28,24,28,28,30,30,26,26,24,26,24,26,28,26,26,18,18,16,18,16,18,20,22,22,18,18,20,18,16,18,20,22,22,18,18,16,18,20,18,20,18,18,22,22,24,22,24,22,20,22,22,30,32,30,28,26,26,30,30,32,30,28,30,20,18,18,22,22,24,22,24,22,20,22,22,30,30,32,30,32,30,28,26,26,30,30,28,26,24,26,28,30,30,26,26,24,26,28,26};

  static const int num_toast_indices = 513;
}

// returns false if the type didn't exist and was skipped (not necessarily an error!)
bool vw::platefile::make_toast_dem_tile(
    const ToastDemWriter& writer,
    const PlateFile& platefile, int32 col, int32 row, int32 level, int32 transaction_id) {

  typedef PixelGrayA<int16> Pixel;

  // First, we need to determine which set of triangle indices to use
  // for this tile.  For upper right and lower left quadrants, we use
  // the first set.  For upper left and lower right quadrants, we use
  // the second set.  This ensures that the long leg of the triangles
  // that cut across the square tile are oriented such that the match
  // up with the equator.
  int const* u_toast_indices;
  int const* v_toast_indices;
  int midpoint = pow(2,level)/2;
  if ( (col >= midpoint && row < midpoint) ||   // upper right
       (col < midpoint && row >= midpoint) ) {  // lower left
    u_toast_indices = ur_ll_quadrant_u_toast_indices;
    v_toast_indices = ur_ll_quadrant_v_toast_indices;
  } else { // upper left or lower right
    u_toast_indices = ul_lr_quadrant_u_toast_indices;
    v_toast_indices = ul_lr_quadrant_v_toast_indices;
  }

  // Next, we need to determine how many levels of difference there is
  // between the image tile size and the DEM tile size.  The latter is
  // 32x32.
  VW_DEBUG_ASSERT(platefile->default_tile_size() > 32,
                  LogicErr() << "Platefile tile size must be larger than 32x32.");

  // How far do we move up the tree to find the already-subsampled images?
  //int level_difference = log(platefile.default_tile_size()/32) / log(2);
  // TODO: use this to skip the interpolation
  int level_difference = 0;

  ImageView<Pixel> tile;
  try {
    // Read the tile & prepare an interpolation view for sampling it.
    platefile.read(tile, col, row, level, transaction_id);
  } catch (TileNotFoundErr &e) {
    // Do nothing if the tile does not exist
    return false;
  }

  const uint64 data_size = num_toast_indices*2;

  // Reduce heap pressure by allocating this up here.
  boost::shared_array<uint8> data(new uint8[data_size]);

  InterpolationView<ImageView<Pixel>, BilinearInterpolation> interp_tile(tile,
                                                                         BilinearInterpolation());
  uint32 dem_level = level + level_difference;
  uint32 tile_diff = pow(2,level_difference);

  union I16 {
    int16 i16;
    uint8 u8[2];
  };
  // Just make sure no funny alignment games are going on.
  BOOST_STATIC_ASSERT(sizeof(I16) == 2);

  for (unsigned v = 0; v < tile_diff; ++v) {
    for (unsigned u = 0; u < tile_diff; ++u) {
      int32 dem_col = col * tile_diff + u;
      int32 dem_row = row * tile_diff + v;

      // Iterate over the triangle vertex arrays above, writing the DEM
      // values in INTEL byte order to disk.
      for (int32 i = 0; i < num_toast_indices; ++i) {
        float u_sample_index = u * ((interp_tile.cols() - 1) / float(tile_diff)) + 
          float(u_toast_indices[i]) / 32.0 * (interp_tile.cols() - 1) / float(tile_diff);
        float v_sample_index = v * ((interp_tile.rows() - 1) / float(tile_diff)) +
          float(v_toast_indices[i]) / 32.0 * (interp_tile.rows() - 1) / float(tile_diff);

        I16 value;

        // TODO: Think about what to do if the pixel has alpha!
        value.i16 = interp_tile( u_sample_index, v_sample_index ).v();
#if VW_BYTE_ORDER == VW_BIG_ENDIAN
        // spec for toast dem files says "intel" byte order (little-endian)
        // so if we're on big endian, swap them.
        std::swap(value.u8[0], value.u8[1]);
#endif
        // write lsb
        data[i*2]   = value.u8[0];
        // write msb
        data[i*2+1] = value.u8[1];
      }

      // We've batched a dem tile worth of data. Call the writer.
      writer(data, data_size, dem_col, dem_row, dem_level, transaction_id);
    }
  }
  return true;
}

namespace {
  struct DemFilesystem : public vw::platefile::ToastDemWriter {

    const std::string& base_output_name;

    DemFilesystem(const std::string& base_output_name) : base_output_name(base_output_name) {}

    void operator()(const boost::shared_array<uint8> data, uint64 data_size, int32 dem_col, int32 dem_row, int32 dem_level, int32 transaction_id) const {
      // Create the level directory (if it doesn't exist)
      std::ostringstream ostr;
      ostr << base_output_name
        << "/" << dem_level
        << "/" << dem_col
        << "/" << dem_row;

      fs::path output_file(ostr.str());
      fs::create_directories(output_file.branch_path());

      // Open the file for writing
      std::ofstream of(output_file.file_string().c_str(), std::ios::binary);
      if (!of.is_open())
        vw_throw(IOErr() << "Failed to open toast dem tile for writing");

      of.write(reinterpret_cast<const char*>(data.get()), data_size);
    }
  };
}

void vw::platefile::save_toast_dem_tile(std::string base_output_name, 
                                        boost::shared_ptr<PlateFile> platefile,
                                        int32 col, int32 row, int32 level,
                                        int32 transaction_id) {

  DemFilesystem writer(base_output_name);
  make_toast_dem_tile(writer, *platefile, col, row, level, transaction_id);
}
