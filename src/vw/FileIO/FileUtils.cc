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


/// \file FileUtils.cc
///
/// An abstract base class referring to an image on disk.
///

#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <vw/FileIO/FileUtils.h>
#include <vw/Core/Log.h>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

namespace fs = boost::filesystem;

// Make the specified file to be relative to the specified directory.
fs::path vw::make_file_relative_to_dir(fs::path const file, fs::path const dir) {
  if (file.has_root_path()){
    if (file.root_path() != dir.root_path()) {
      return file;
    } else {
      return make_file_relative_to_dir(file.relative_path(), dir.relative_path());
    }
  } else {
    if (dir.has_root_path()) {
      fs::path file2 = fs::complete(file);
      return make_file_relative_to_dir(file2.relative_path(), dir.relative_path());
    } else {
      typedef fs::path::const_iterator path_iterator;
      path_iterator file_it = file.begin();
      path_iterator dir_it = dir.begin();
      while ( file_it != file.end() && dir_it != dir.end() ) {
        if (*file_it != *dir_it) break;
        ++file_it; ++dir_it;
      }
      fs::path result;
      for (; dir_it != dir.end(); ++dir_it) {
        result /= "..";
      }
      for (; file_it != file.end(); ++file_it) {
        result /= *file_it;
      }
      return result;
    }
  }
}

// Remove file name extension
std::string vw::prefix_from_filename(std::string const& filename) {
  std::string result = filename;
  int index = result.rfind(".");
  if (index != -1)
    result.erase(index, result.size());
  return result;
}

std::string vw::get_extension( std::string const& input ) {
  boost::filesystem::path ipath( input );
  std::string ext = ipath.extension().string();
  boost::algorithm::to_lower(ext);
  return ext;
}

// If prefix is "dir/out", create directory "dir"
void vw::create_out_dir(std::string out_prefix){

  fs::path out_prefix_path(out_prefix);
  if (out_prefix_path.has_parent_path()) {
    if (!fs::is_directory(out_prefix_path.parent_path())) {
      vw_out() << "\nCreating output directory: "
               << out_prefix_path.parent_path() << std::endl;
      fs::create_directories(out_prefix_path.parent_path());
    }
  }

  return;
}
