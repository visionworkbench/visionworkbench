// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/ToastDem.h>
#include <vw/Plate/TileManipulation.h>
#include <vw/Plate/PlateFile.h>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>

using namespace vw;
namespace fs = boost::filesystem;

// ---------------------------------------------------------------------------------
//                                WWT DEM Indices
// ---------------------------------------------------------------------------------

namespace {

  // DEM triangle indices for the upper right (0-90 degrees) and lower
  // left (180-270 degrees) quadrants of the TOAST space.
  static uint32 const ur_ll_quadrant_u_toast_indices[] = { 0,16,32,0,16,32,0,16,32,8,8,16,8,0,8,12,12,16,12,16,12,8,4,4,12,12,8,4,0,4,8,12,12,4,4,0,4,8,4,14,14,16,14,16,14,12,10,10,14,14,12,14,16,14,12,10,10,14,14,16,14,12,14,4,2,2,6,6,8,6,8,6,4,6,6,10,10,8,10,8,10,12,14,14,10,10,12,2,0,2,4,6,6,2,2,0,2,4,2,12,14,14,10,10,8,10,8,10,12,10,10,2,2,0,2,0,2,4,6,6,2,2,4,6,8,6,4,2,2,6,6,8,6,4,6,24,24,32,24,16,24,28,28,32,28,32,28,24,20,20,28,28,24,20,16,20,24,28,28,20,20,16,20,24,20,30,30,32,30,32,30,28,26,26,30,30,28,30,32,30,28,26,26,30,30,32,30,28,30,20,18,18,22,22,24,22,24,22,20,22,22,26,26,24,26,24,26,28,30,30,26,26,28,18,16,18,20,22,22,18,18,16,18,20,18,28,30,30,26,26,24,26,24,26,28,26,26,18,18,16,18,16,18,20,22,22,18,18,20,22,24,22,20,18,18,22,22,24,22,20,22,8,8,16,8,0,8,12,12,16,12,16,12,8,4,4,12,12,8,4,0,4,8,12,12,4,4,0,4,8,4,14,14,16,14,16,14,12,10,10,14,14,12,14,16,14,12,10,10,14,14,16,14,12,14,4,2,2,6,6,8,6,8,6,4,6,6,10,10,8,10,8,10,12,14,14,10,10,12,2,0,2,4,6,6,2,2,0,2,4,2,12,14,14,10,10,8,10,8,10,12,10,10,2,2,0,2,0,2,4,6,6,2,2,4,6,8,6,4,2,2,6,6,8,6,4,6,24,24,32,24,16,24,28,28,32,28,32,28,24,20,20,28,28,24,20,16,20,24,28,28,20,20,16,20,24,20,30,30,32,30,32,30,28,26,26,30,30,28,30,32,30,28,26,26,30,30,32,30,28,30,20,18,18,22,22,24,22,24,22,20,22,22,26,26,24,26,24,26,28,30,30,26,26,28,18,16,18,20,22,22,18,18,16,18,20,18,28,30,30,26,26,24,26,24,26,28,26,26,18,18,16,18,16,18,20,22,22,18,18,20,22,24,22,20,18,18,22,22,24,22,20,22};

  static uint32 const ur_ll_quadrant_v_toast_indices[] = {0,0,0,16,16,16,32,32,32,0,8,8,8,8,16,8,12,12,4,4,0,4,0,4,8,4,4,12,12,16,12,16,12,8,4,4,12,12,8,12,14,14,10,10,8,10,8,10,12,10,10,2,2,0,2,0,2,4,6,6,2,2,4,2,0,2,4,6,6,2,2,0,2,4,2,4,2,2,6,6,8,6,8,6,4,6,6,14,14,16,14,16,14,12,10,10,14,14,12,14,16,14,12,10,10,14,14,16,14,12,14,4,2,2,6,6,8,6,8,6,4,6,6,10,10,8,10,8,10,12,14,14,10,10,12,0,8,8,8,8,16,8,12,12,4,4,0,4,0,4,8,4,4,12,12,16,12,16,12,8,4,4,12,12,8,12,14,14,10,10,8,10,8,10,12,10,10,2,2,0,2,0,2,4,6,6,2,2,4,2,0,2,4,6,6,2,2,0,2,4,2,4,2,2,6,6,8,6,8,6,4,6,6,14,14,16,14,16,14,12,10,10,14,14,12,14,16,14,12,10,10,14,14,16,14,12,14,4,2,2,6,6,8,6,8,6,4,6,6,10,10,8,10,8,10,12,14,14,10,10,12,16,24,24,24,24,32,24,28,28,20,20,16,20,16,20,24,20,20,28,28,32,28,32,28,24,20,20,28,28,24,28,30,30,26,26,24,26,24,26,28,26,26,18,18,16,18,16,18,20,22,22,18,18,20,18,16,18,20,22,22,18,18,16,18,20,18,20,18,18,22,22,24,22,24,22,20,22,22,30,30,32,30,32,30,28,26,26,30,30,28,30,32,30,28,26,26,30,30,32,30,28,30,20,18,18,22,22,24,22,24,22,20,22,22,26,26,24,26,24,26,28,30,30,26,26,28,16,24,24,24,24,32,24,28,28,20,20,16,20,16,20,24,20,20,28,28,32,28,32,28,24,20,20,28,28,24,28,30,30,26,26,24,26,24,26,28,26,26,18,18,16,18,16,18,20,22,22,18,18,20,18,16,18,20,22,22,18,18,16,18,20,18,20,18,18,22,22,24,22,24,22,20,22,22,30,30,32,30,32,30,28,26,26,30,30,28,30,32,30,28,26,26,30,30,32,30,28,30,20,18,18,22,22,24,22,24,22,20,22,22,26,26,24,26,24,26,28,30,30,26,26,28};

  // DEM triangle indices for the upper left (90-180 degrees) and lower
  // left (270-360 degrees) quadrants of the TOAST space.
  static uint32 const ul_lr_quadrant_u_toast_indices[] = {0,16,32,0,16,32,0,16,32,8,0,8,8,8,16,4,0,4,8,12,12,4,4,0,4,8,4,12,12,16,12,16,12,8,4,4,12,12,8,2,0,2,4,6,6,2,2,0,2,4,2,12,14,14,10,10,8,10,8,10,12,10,10,2,2,0,2,0,2,4,6,6,2,2,4,6,8,6,4,2,2,6,6,8,6,4,6,14,14,16,14,16,14,12,10,10,14,14,12,14,16,14,12,10,10,14,14,16,14,12,14,4,2,2,6,6,8,6,8,6,4,6,6,10,10,8,10,8,10,12,14,14,10,10,12,24,16,24,24,24,32,20,16,20,24,28,28,20,20,16,20,24,20,28,28,32,28,32,28,24,20,20,28,28,24,18,16,18,20,22,22,18,18,16,18,20,18,28,30,30,26,26,24,26,24,26,28,26,26,18,18,16,18,16,18,20,22,22,18,18,20,22,24,22,20,18,18,22,22,24,22,20,22,30,30,32,30,32,30,28,26,26,30,30,28,30,32,30,28,26,26,30,30,32,30,28,30,20,18,18,22,22,24,22,24,22,20,22,22,26,26,24,26,24,26,28,30,30,26,26,28,8,0,8,8,8,16,4,0,4,8,12,12,4,4,0,4,8,4,12,12,16,12,16,12,8,4,4,12,12,8,2,0,2,4,6,6,2,2,0,2,4,2,12,14,14,10,10,8,10,8,10,12,10,10,2,2,0,2,0,2,4,6,6,2,2,4,6,8,6,4,2,2,6,6,8,6,4,6,14,14,16,14,16,14,12,10,10,14,14,12,14,16,14,12,10,10,14,14,16,14,12,14,4,2,2,6,6,8,6,8,6,4,6,6,10,10,8,10,8,10,12,14,14,10,10,12,24,16,24,24,24,32,20,16,20,24,28,28,20,20,16,20,24,20,28,28,32,28,32,28,24,20,20,28,28,24,18,16,18,20,22,22,18,18,16,18,20,18,28,30,30,26,26,24,26,24,26,28,26,26,18,18,16,18,16,18,20,22,22,18,18,20,22,24,22,20,18,18,22,22,24,22,20,22,30,30,32,30,32,30,28,26,26,30,30,28,30,32,30,28,26,26,30,30,32,30,28,30,20,18,18,22,22,24,22,24,22,20,22,22,26,26,24,26,24,26,28,30,30,26,26,28};

  static uint32 const ul_lr_quadrant_v_toast_indices[] = {0,0,0,16,16,16,32,32,32,0,8,8,8,16,8,8,12,12,4,4,0,4,0,4,8,4,4,12,16,12,8,4,4,12,12,16,12,8,12,12,14,14,10,10,8,10,8,10,12,10,10,2,2,0,2,0,2,4,6,6,2,2,4,2,0,2,4,6,6,2,2,0,2,4,2,4,2,2,6,6,8,6,8,6,4,6,6,14,16,14,12,10,10,14,14,16,14,12,14,4,2,2,6,6,8,6,8,6,4,6,6,14,14,16,14,16,14,12,10,10,14,14,12,10,8,10,12,14,14,10,10,8,10,12,10,0,8,8,8,16,8,8,12,12,4,4,0,4,0,4,8,4,4,12,16,12,8,4,4,12,12,16,12,8,12,12,14,14,10,10,8,10,8,10,12,10,10,2,2,0,2,0,2,4,6,6,2,2,4,2,0,2,4,6,6,2,2,0,2,4,2,4,2,2,6,6,8,6,8,6,4,6,6,14,16,14,12,10,10,14,14,16,14,12,14,4,2,2,6,6,8,6,8,6,4,6,6,14,14,16,14,16,14,12,10,10,14,14,12,10,8,10,12,14,14,10,10,8,10,12,10,16,24,24,24,32,24,24,28,28,20,20,16,20,16,20,24,20,20,28,32,28,24,20,20,28,28,32,28,24,28,28,30,30,26,26,24,26,24,26,28,26,26,18,18,16,18,16,18,20,22,22,18,18,20,18,16,18,20,22,22,18,18,16,18,20,18,20,18,18,22,22,24,22,24,22,20,22,22,30,32,30,28,26,26,30,30,32,30,28,30,20,18,18,22,22,24,22,24,22,20,22,22,30,30,32,30,32,30,28,26,26,30,30,28,26,24,26,28,30,30,26,26,24,26,28,26,16,24,24,24,32,24,24,28,28,20,20,16,20,16,20,24,20,20,28,32,28,24,20,20,28,28,32,28,24,28,28,30,30,26,26,24,26,24,26,28,26,26,18,18,16,18,16,18,20,22,22,18,18,20,18,16,18,20,22,22,18,18,16,18,20,18,20,18,18,22,22,24,22,24,22,20,22,22,30,32,30,28,26,26,30,30,32,30,28,30,20,18,18,22,22,24,22,24,22,20,22,22,30,30,32,30,32,30,28,26,26,30,30,28,26,24,26,28,30,30,26,26,24,26,28,26};

  static const size_t num_toast_indices = 513;

  // Access individual bytes without violating aliasing by using a pointer
  union I16 {
    int16 i16;
    uint8 u8[2];
  };
  // Just make sure no funny alignment games are going on.
  BOOST_STATIC_ASSERT(sizeof(I16) == 2);

  union F32 {
    float32 f32;
    uint8 u8[4];
  };
  BOOST_STATIC_ASSERT(sizeof(F32) == 4);

}

namespace vw {
namespace platefile {

  template <class PixelT>
  bool toast_dem_tile_helper(const ToastDemWriter& writer,
                             const PlateFile& platefile,
                             int32 col, int32 row, int32 level,
                             int32 level_difference,
                             TransactionOrNeg input_transaction_id,
                             Transaction output_transaction_id) {

    // First, we need to determine which set of triangle indices to use
    // for this tile.  For upper right and lower left quadrants, we use
    // the first set.  For upper left and lower right quadrants, we use
    // the second set.  This ensures that the long leg of the triangles
    // that cut across the square tile are oriented such that the match
    // up with the equator.
    uint32 const* u_toast_indices;
    uint32 const* v_toast_indices;
    int32 midpoint = 1 << (level-1);
    if ( (col >= midpoint && row < midpoint) ||   // upper right
         (col < midpoint && row >= midpoint) ) {  // lower left
      u_toast_indices = ur_ll_quadrant_u_toast_indices;
      v_toast_indices = ur_ll_quadrant_v_toast_indices;
    } else { // upper left or lower right
      u_toast_indices = ul_lr_quadrant_u_toast_indices;
      v_toast_indices = ul_lr_quadrant_v_toast_indices;
    }

    ImageView<PixelT> src_tile;
    try {
      // Read the tile & prepare the interpolation view for sampling it.
      platefile.read(src_tile, col, row, level, input_transaction_id);
    } catch (TileNotFoundErr &e) {
      // Do nothing if the tile does not exist
      return false;
    }

    VW_ASSERT( src_tile.cols() == 257 && src_tile.rows() == 257,
               ArgumentErr() << "Could not process " << src_tile.cols() << "x" << src_tile.rows()
               << "tile.  Source data for toast DEM tiles must be 257x257 pixels.");

    const uint32 src_level = level;
    const uint32 dst_level = src_level + 3;
    const uint32 dst_region_offset = 8;

    // Reduce heap pressure by allocating this up here and reusing
    static const size_t tile_bytes = num_toast_indices*sizeof(typename PixelChannelType<PixelT>::type);
    boost::shared_array<uint8> data(new uint8[tile_bytes]);
    for (size_t i = 0; i < tile_bytes; ++i) {
      data[i] = 0;
    }

    // Iterate over the 64 subtiles.
    for (uint32 subtile_r = 0; subtile_r < 8; ++subtile_r) {
      for (uint32 subtile_c = 0; subtile_c < 8; ++subtile_c) {
        const uint32 dst_row = row * dst_region_offset + subtile_r;
        const uint32 dst_col = col * dst_region_offset + subtile_c;

        // For the top level tiles (where level_difference == 0), we
        // override the u_base and v_base and start from the
        // beginning of the tile.
        const uint32 u_base = level_difference == 0 ? 0 : (subtile_c * 32);
        const uint32 v_base = level_difference == 0 ? 0 : (subtile_r * 32);

        for (size_t i = 0; i < num_toast_indices; ++i) {
          // Again for the top three levels of tiles, we do some
          // special casing here.  We aren't doing any subtiles, so we
          // instead create a toast dem tile from a complete image
          // tile by decimating it.
          const uint32 u = level_difference == 0 ? (u_base + u_toast_indices[i]*8)
                                                 : u_base + u_toast_indices[i];
          const uint32 v = level_difference == 0 ? (v_base + v_toast_indices[i]*8)
                                                 : v_base + v_toast_indices[i];

          // Handle 16-bit integer data types
          if (sizeof(typename PixelChannelType<PixelT>::type) == 2) {
            // TODO: Think about what to do if the pixel has alpha!
            I16 value;
            value.i16 = static_cast<int16>(src_tile(u,v).v());

            // spec for toast dem files says "intel" byte order (little-endian)
            // so if we're on big endian, swap them.
#if VW_BYTE_ORDER == VW_BIG_ENDIAN
            std::swap(value.u8[0], value.u8[1]);
#endif
            // write lsb
            data[i*2]   = value.u8[0];
            // write msb
            data[i*2+1] = value.u8[1];

          // Handle 32-bit floating point data types
          } else if (sizeof(typename PixelChannelType<PixelT>::type) == 4) {
            F32 value;
            value.f32 = src_tile(u,v).v();

            // Write the float to the data stream
            data[i*4]   = value.u8[0];
            data[i*4+1] = value.u8[1];
            data[i*4+2] = value.u8[2];
            data[i*4+3] = value.u8[3];
          }
        }

        // We've batched a dem tile worth of data. Call the writer.
        if (level_difference == 0) {
          writer(data, tile_bytes, col, row, level, output_transaction_id);
        } else {
          writer(data, tile_bytes, dst_col, dst_row, dst_level, output_transaction_id);
        }

      }
    }

    return true;

    // -------------------------- OLD CODE ------------------------
    // NOTE: This code worked a little bit more "automatically", but
    // the interpolation that it did resulted in some unsightly
    // hairline seams in the WWT terrain.  The above code is much more
    // custom-built to WWT terrain generation and avoids the
    // interpolation problem, but at the expense of generality.  (the
    // generality is hardly really necessary here, though...)

    //     int src_level = level;
    //     int dst_level = src_level + level_difference;
    //     int dst_region_offset = 1 << level_difference;

    //     // Reduce heap pressure by allocating this up here and reusing
    //     static const uint64 tile_bytes = num_toast_indices*sizeof(typename PixelChannelType<PixelT>::type);
    //     boost::shared_array<uint8> data(new uint8[tile_bytes]);

    //     // Iterate over the 64 subtiles.
    //     for (int v = 0; v < dst_region_offset; ++v) {
    //       for (int u = 0; u < dst_region_offset; ++u) {
    //         int dst_col = col * dst_region_offset + u;
    //         int dst_row = row * dst_region_offset + v;

    //         float slice_cols = (float(src_tile.cols())-1.0) / dst_region_offset;
    //         float slice_rows = (float(src_tile.rows())-1.0) / dst_region_offset;

    //         float subtile_base_col = slice_cols * u;
    //         float subtile_base_row = slice_rows * v;

    //         // Create an interpolation view to sample from.
    //         InterpolationView<EdgeExtensionView<ImageView<PixelT>, ConstantEdgeExtension>,
    //           BilinearInterpolation> interp_img = interpolate(src_tile,
    //                                                           BilinearInterpolation(),
    //                                                           ConstantEdgeExtension());

    //         // Iterate over the triangle vertex arrays above, writing the DEM
    //         // values in INTEL byte order to disk.
    //         // std::cout << "--> " << u << " " << v << " -- "
    //         //           << subtile_base_col << " " << subtile_base_row << "\n";
    //         for (int32 i = 0; i < num_toast_indices; ++i) {
    //           float u_sample_index = float(u_toast_indices[i]) * slice_cols / 32.0;
    //           float v_sample_index = float(v_toast_indices[i]) * slice_rows / 32.0;

    //           // std::cout << "[" << (subtile_base_col + u_sample_index) << " "
    //           //           << (subtile_base_row + v_sample_index) << "]   ";

    //           // Handle 16-bit integer data types
    //           if (sizeof(typename PixelChannelType<PixelT>::type) == 2) {
    //             // TODO: Think about what to do if the pixel has alpha!
    //             I16 value;
    //             value.i16 = interp_img( subtile_base_col + u_sample_index,
    //                                     subtile_base_row + v_sample_index ).v();
    // #if VW_BYTE_ORDER == VW_BIG_ENDIAN
    //             // spec for toast dem files says "intel" byte order (little-endian)
    //             // so if we're on big endian, swap them.
    //             std::swap(value.u8[0], value.u8[1]);
    // #endif
    //             // write lsb
    //             data[i*2]   = value.u8[0];
    //             // write msb
    //             data[i*2+1] = value.u8[1];

    //           // Handle 32-bit floating point data types
    //           } else if (sizeof(typename PixelChannelType<PixelT>::type) == 4) {
    //             F32 value;
    //             value.f32 = interp_img( subtile_base_col + u_sample_index,
    //                                     subtile_base_row + v_sample_index ).v();
    //             // Write the float to the data stream
    //             data[i*4]   = value.u8[0];
    //             data[i*4+1] = value.u8[1];
    //             data[i*4+2] = value.u8[2];
    //             data[i*4+3] = value.u8[3];
    //           }
    //         }
    //         // std::cout << "\n";

    //         // We've batched a dem tile worth of data. Call the writer.
    //         writer(data, tile_bytes, dst_col, dst_row, dst_level, output_transaction_id);
    //     }
    //   }
    //    return true;
}
}} // namespace vw::platefile

// returns false if the type didn't exist and was skipped (not necessarily an error!)
bool vw::platefile::make_toast_dem_tile(const ToastDemWriter& writer,
                                        const PlateFile& platefile,
                                        int32 col, int32 row, int32 level,
                                        int32 level_difference,
                                        TransactionOrNeg input_transaction_id,
                                        Transaction output_transaction_id) {

  if (platefile.channel_type() == VW_CHANNEL_INT16) {
    return toast_dem_tile_helper<PixelGrayA<int16> >(writer, platefile, col, row, level,
                                                     level_difference,
                                                     input_transaction_id, output_transaction_id);
  } else if (platefile.channel_type() == VW_CHANNEL_FLOAT32) {
    return toast_dem_tile_helper<PixelGrayA<float32> >(writer, platefile, col, row, level,
                                                       level_difference,
                                                       input_transaction_id, output_transaction_id);
  } else {
    std::cout << "Could not convert to toast_dem.  "
              << "Unsupported channel type in platefile: " << platefile.channel_type() << "\n";
    exit(0);
  }
  return false;

}

namespace {
  struct DemFilesystem : public vw::platefile::ToastDemWriter {

    const std::string& base_output_name;

    DemFilesystem(const std::string& base_output_name) : base_output_name(base_output_name) {}

    void operator()(const boost::shared_array<uint8> data, uint64 data_size, int32 dem_col, int32 dem_row, int32 dem_level, vw::platefile::Transaction /*transaction_id*/) const {
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
                                        int32 level_difference,
                                        TransactionOrNeg transaction_id) {

  DemFilesystem writer(base_output_name);
  make_toast_dem_tile(writer, *platefile, col, row, level, level_difference,
                      transaction_id, 0); // output_transaction_id doesn't matter here.
}

boost::shared_array<uint8> vw::platefile::toast_dem_null_tile(uint64& output_tile_size) {

    static const uint64 tile_bytes = num_toast_indices*2;
    boost::shared_array<uint8> null_tile(new uint8[tile_bytes]);
    memset(null_tile.get(), 0, tile_bytes);
    output_tile_size = tile_bytes;
    return null_tile;
}
