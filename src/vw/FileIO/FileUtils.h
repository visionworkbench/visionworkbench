/// \file FileIO/Fileutils.h
///
/// An abstract base class referring to an image on disk.
///
#ifndef __VW_FILEIO_FILEUTILS_H__
#define __VW_FILEIO_FILEUTILS_H__

#include <string>
#include <boost/filesystem/path.hpp>

namespace vw{

  boost::filesystem::path make_file_relative_to_dir(boost::filesystem::path const file,
                                                    boost::filesystem::path const dir);

  /// Remove file name extension
  std::string prefix_from_filename(std::string const& filename);

  /// Return lower-case extension, such as .cub (so the dot is also present).
  std::string get_extension(std::string const& input, bool make_lower = true);

  /// If prefix is "dir/out", create directory "dir"
  void create_out_dir(std::string out_prefix);

  /// Return true if the image_file/camera file combination represents a SPOT5 camera file.
  /// - Returns false if the camera_file input is empty.
  bool has_spot5_extension(std::string const& image_file, std::string const& camera_file="");
  
  
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
