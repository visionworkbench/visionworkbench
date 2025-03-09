/// \file FileIO/FileTypes.h
///
/// Functions to handle file types and extensions.
///
#ifndef __VW_FILEIO_FILETYPES_H__
#define __VW_FILEIO_FILETYPES_H__

#include <string>
#include <vector>

namespace vw {

/// Return lower-case extension, such as .cub (so the dot is also present).
std::string get_extension(std::string const& input, bool make_lower = true);

/// Return true if the image_file/camera file combination represents a SPOT5 camera file.
/// - Returns false if the camera_file input is empty.
bool has_spot5_extension(std::string const& image_file, std::string const& camera_file="");

/// Returns true if the file has an extension which can contain a pinhole camera model
bool has_pinhole_extension(std::string const& input);

/// Returns true if the file has an extension which is tif or ntf
bool has_tif_or_ntf_extension(std::string const& input);

/// Returns true for a shapefile
bool has_shp_extension(std::string const& input);

/// If it ends with _rpc.txt or _RPC.TXT
bool has_rpc_txt_extension(std::string const& input);

/// Returns true if the filename ends in .isd or .json.
bool has_isd_extension(std::string const& path);

/// Returns true if the file has an extension which can contain an image
bool has_image_extension(std::string const& input);

/// Returns true if the file has an extension which can contain a camera model
bool has_cam_extension(std::string const& input);

/// Makes a vector containing all files in the input vector with an extension.
/// - If prune_input_list is set, matching files are removed from the input list.
std::vector<std::string>
get_files_with_ext(std::vector<std::string>& files, std::string const& ext, 
                   bool prune_input_list);

} // namespace vw

#endif // __VW_FILETYPES_H__
