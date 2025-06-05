/// \file FileIO/Fileutils.h
///
/// Functions to handle files.
///
#ifndef __VW_FILEIO_FILEUTILS_H__
#define __VW_FILEIO_FILEUTILS_H__

#include <string>
#include <vector>
#include <boost/filesystem/path.hpp>

namespace vw{

  boost::filesystem::path make_file_relative_to_dir(boost::filesystem::path const file,
                                                    boost::filesystem::path const dir);

  /// Remove file name extension
  std::string prefix_from_filename(std::string const& filename);

  /// If prefix is "dir/out", create directory "dir"
  void create_out_dir(std::string out_prefix);

  /// Strip the directory out of a file path
  std::string strip_directory( std::string const& input);

  /// Strip the directory and extension out of a file path
  std::string strip_directory_and_extension( std::string const& input);
  
  /// Returns the folder of the provided path.
  std::string get_folder(std::string const& input);

  /// Populate a list of all files in a directory.
  /// - Returns the number of files found.
  /// - If an extension is passed in, files must match the extension.
  size_t get_files_in_folder(std::string              const& folder,
                            std::vector<std::string>      & output,
                            std::string              const& ext="");

  /// Find all of the files that start with the provided prefix.
  /// - Returns the number of files found.
  size_t get_files_with_prefix(std::string              const& prefix,
                               std::vector<std::string>      & output);
  
} // namespace vw

#endif // __VW_FILEIO_FILEUTILS_H__
