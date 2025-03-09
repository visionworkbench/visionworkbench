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

namespace vw {

// Make the specified file to be relative to the specified directory.
fs::path make_file_relative_to_dir(fs::path const file, fs::path const dir) {
  if (file.has_root_path()){
    if (file.root_path() != dir.root_path()) {
      return file;
    } else {
      return make_file_relative_to_dir(file.relative_path(), dir.relative_path());
    }
  } else {
    if (dir.has_root_path()) {
      fs::path file2 = fs::absolute(file);
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

std::string prefix_from_filename(std::string const& filename) {
  std::string result = filename;
  int index = result.rfind(".");
  if (index != -1)
    result.erase(index, result.size());
  return result;
}

extern "C" void openblas_set_num_threads(int num_threads);

void init_blas (){
  // Must use just one thread in BLAS in each thread
  // of our software.
  openblas_set_num_threads(1);
}
  
// If prefix is "dir/out", create directory "dir".
// If the output prefix does not contain a directory, we will do nothing,
// as the output will be written to the current directory.
void create_out_dir(std::string out_prefix){

  fs::path out_prefix_path(out_prefix);
  if (out_prefix_path.has_parent_path()) {
    if (!fs::is_directory(out_prefix_path.parent_path())) {
      vw_out() << "\nCreating output directory: "
               << out_prefix_path.parent_path() << ".\n";
      fs::create_directories(out_prefix_path.parent_path());
    }
  } 

  // It is convenient to put this code here as this function is called
  // from many places. 
  init_blas();
  
  return;
}

std::string strip_directory( std::string const& input){
 boost::filesystem::path p(input);
 return p.filename().string();
}

std::string strip_directory_and_extension( std::string const& input){
 boost::filesystem::path p(input);
 return p.stem().string();
}


std::string get_folder(std::string const& input) {
  boost::filesystem::path p(input);
  if (boost::filesystem::is_directory(p)) // Handle inputs like "/usr/local"
    return input;
  return p.parent_path().string();
}

size_t get_files_in_folder(std::string              const& folder,
                           std::vector<std::string>      & output,
                           std::string              const& ext)
{
  output.clear();
  
  // Handle invalid inputs
  if(!boost::filesystem::exists(folder) || !boost::filesystem::is_directory(folder)) 
    return 0;

  boost::filesystem::directory_iterator it(folder);
  boost::filesystem::directory_iterator endit;

  if (ext != ""){ // Check the extension
    while(it != endit) {
        if(boost::filesystem::is_regular_file(*it) && it->path().extension() == ext) 
          output.push_back(it->path().filename().string());
        ++it;
    }
  }
  else{ // No extension check
    while(it != endit) {
        if(boost::filesystem::is_regular_file(*it)) 
          output.push_back(it->path().filename().string());
        ++it;
    }
  }
  return output.size();
}


size_t get_files_with_prefix(std::string              const& prefix,
                               std::vector<std::string>    & output) {
  output.clear();
  
  // Get a list of all the files in the same folder
  std::string folder = get_folder(prefix);
  std::vector<std::string> all_files;
  get_files_in_folder(folder, all_files);

  std::string name = strip_directory(prefix);

  // Find all the files that match the prefix
  for (size_t i=0; i<all_files.size(); ++i) {
    if (all_files[i].find(name) == 0) {
      fs::path p(folder);
      p /= all_files[i]; // Return the full path to the file
      output.push_back(p.string());
    }
  }
  return output.size();
}

} // End namespace vw
