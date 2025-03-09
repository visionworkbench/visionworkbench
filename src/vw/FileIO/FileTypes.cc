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

/// \file FileTypes.cc
///
/// Functions to handle file types and extensions.
///

#include <vw/FileIO/FileTypes.h>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

namespace fs = boost::filesystem;

namespace vw {

std::string get_extension(std::string const& input, bool make_lower) {
  boost::filesystem::path ipath( input );
  std::string ext = ipath.extension().string();
  if (make_lower)
    boost::algorithm::to_lower(ext);
  return ext;
}

bool has_spot5_extension(std::string const& image_file, std::string const& camera_file) {
  // First check the image file
  std::string image_ext = vw::get_extension(image_file);
  boost::algorithm::to_lower(image_ext);
  if ((image_ext == ".bip") || (image_ext == ".bil") || (image_ext == ".bsq"))
    return true;
  // If no camera file was provided it cannot be a Spot5 file
  if ((camera_file.empty()) || (camera_file == image_file))
    return false;
  // The Spot5 file is the last thing we check
  const std::string camera_ext = vw::get_extension(camera_file);
  return ((camera_ext == ".DIM") || (camera_ext == ".dim"));
}

bool has_pinhole_extension(std::string const& input) {
  std::string ext = get_extension(input);
  if (ext == ".cahvor"  || ext == ".cahv"    ||
       ext == ".pin"     || ext == ".pinhole" ||
       ext == ".tsai"    || ext == ".cmod"    ||
       ext == ".cahvore")
    return true;
  return false;
}

bool has_tif_or_ntf_extension(std::string const& input) {
  std::string ext = vw::get_extension(input);
  if (ext == ".tif"  || ext == ".tiff" || 
      ext == ".ntf" || ext == ".nitf")
    return true;
  return false;
}

bool has_shp_extension(std::string const& input) {
  std::string ext = vw::get_extension(input);
  if (ext == ".shp")
    return true;
  return false;
}

// If it ends with _rpc.txt or _RPC.TXT
bool has_rpc_txt_extension(std::string const& input) {
  if (boost::iends_with(input, "_rpc.txt")) 
    return true;
  return false;
}

bool has_isd_extension(std::string const& path) {
  std::string ext = vw::get_extension(path);
  return ((ext == ".json") || (ext == ".isd"));
}

bool has_image_extension(std::string const& input) {
  std::string ext = vw::get_extension(input);
  if (ext == ".tif"  || ext == ".tiff" || 
      ext == ".ntf" || ext == ".nitf"  ||
      ext == ".png"  || ext == ".jpeg" ||
      ext == ".jpg"  || ext == ".jp2"  ||
      ext == ".img"  || ext == ".cub"  ||
      ext == ".bip"  || ext == ".bil"  || 
      ext == ".bsq")
    return true;
  return false;
}

bool has_cam_extension(std::string const& input) {
  std::string ext = vw::get_extension(input);
  if (vw::has_pinhole_extension(input) ||
      vw::has_isd_extension(input)     ||
      ext == ".cub" || ext == ".xml" || ext == ".dim" ||
      ext == ".rpb" || vw::has_rpc_txt_extension(input))
    return true;
  return false;
}

std::vector<std::string>
get_files_with_ext(std::vector<std::string>& files, std::string const& ext,
                   bool prune_input_list) {
  std::vector<std::string> files_with_ext;
  std::vector<std::string>::iterator it = files.begin();
  while (it != files.end()) {
    if (boost::iends_with(boost::to_lower_copy(*it), ext)) { // Match
      files_with_ext.push_back(*it);
      if (prune_input_list) // Clear match from the input list
        it = files.erase(it);
      else
        it++;
    } else {// No match
      it++;
    }
  } // End loop through input list

  return files_with_ext;
}

} // End namespace vw
