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

/// \file MatcherIO.cc
///
/// Functions for reading and writing interest point match files.
///
#include <vw/InterestPoint/MatcherIO.h>
#include <vw/InterestPoint/Matcher.h>
#include <vw/InterestPoint/InterestData.h>

#include <vw/FileIO/FileUtils.h>

#include <boost/filesystem/operations.hpp>

namespace fs = boost::filesystem;

namespace vw { namespace ip {

/// Return the basename with no extension, and shorten it if need be.
/// If filename starts with out_prefix followed by dash, strip those.
std::string strip_path(std::string out_prefix, std::string filename) {

  std::string ss = out_prefix + "-";
  size_t found = filename.find(ss);

  // Check for out_prefix being non-empty, otherwise it gives a wrong answer
  if (out_prefix != "" && found != std::string::npos)
    filename.erase(found, ss.length());

  // Find the basename using boost
  filename = fs::path(filename).filename().string();

  // Find the last dot in the file name. If not found, set it to the length of
  // the string.
  size_t dot = filename.rfind(".");
  if (dot == std::string::npos)
    dot = filename.size();
    
  // Find the substring until the dot
  filename = filename.substr(0, dot);

  return filename;
}

std::string match_filename(std::string const& out_prefix,
                           std::string const& input_file1,
                           std::string const& input_file2) {

  std::string name1 = strip_path(out_prefix, input_file1);
  std::string name2 = strip_path(out_prefix, input_file2);

  // Filenames longer than this must be chopped, as too long names cause
  // problems later with boost. But shorter names may lead to inconsistencies.
  // So give a warning.
  int max_len = 60;
  if (name1.size() >= max_len) {
    vw_out(WarningMessage) << "Warning: Shortening the long file part: " << input_file1 
      << ". In case of problems, use shorter input file names.\n";
    name1 = name1.substr(0, max_len);
  }
  if (name2.size() >= max_len) {
    vw_out(WarningMessage) << "Warning: Shortening the long file part: " 
      << input_file2 << ". In case of problems, use shorter input file names.\n";
    name2 = name2.substr(0, max_len);
  }

  std::string suffix = name1 + "__" + name2 + ".match";
  
  if (out_prefix == "")
    return suffix;
  
  return out_prefix + "-" + suffix;
}

/// Convert match file name to clean match file name.
std::string clean_match_filename(std::string const& match_file) {
  std::string clean_match_file = fs::path(match_file).replace_extension("").string();
  clean_match_file += "-clean.match";
  return clean_match_file;
}

/// The name of the clean match file.
std::string clean_match_filename(std::string const& out_prefix,
                                 std::string const& input_file1,
                                 std::string const& input_file2) {

  return clean_match_filename(match_filename(out_prefix, input_file1, input_file2));
}

std::string ip_filename(std::string const& out_prefix,
                        std::string const& input_file) {
  return out_prefix + "-" + strip_path(out_prefix, input_file) + ".vwip";
}

void ip_filenames(std::string const& out_prefix,
                  std::string const& input_file1,
                  std::string const& input_file2,
                  std::string & output_ip1,
                  std::string & output_ip2){
  output_ip1 = out_prefix + "-" + strip_path(out_prefix, input_file1) + ".vwip";
  output_ip2 = out_prefix + "-" + strip_path(out_prefix, input_file2) + ".vwip";
}

  void write_lowe_ascii_ip_file(std::string ip_file, InterestPointList ip) {
    vw::create_out_dir(ip_file);

    size_t num_pts = ip.size();
    if (num_pts == 0)
      vw_throw(IOErr() << "Attempted to write Lowe SIFT format interest point file with an empty list of interest points.");

    size_t size = ip.front().descriptor.size();

    // Write out detected interest points to file.
    FILE *out = fopen(ip_file.c_str(), "w");
    fprintf(out, "%u %u\n", uint32(num_pts), uint32(size));
    for (InterestPointList::iterator i = ip.begin(); i != ip.end(); ++i) {
      float orientation = i->orientation;
      while (orientation > M_PI) orientation -= 2 * M_PI;
      while (orientation < -M_PI) orientation += 2 * M_PI;
      fprintf(out, "%.2f %.2f %.2f %.3f", i->y, i->x, i->scale, orientation);
      for (size_t element = 0; element < size; ++element) {
        if (element % 20 == 0) fprintf(out, "\n");
        fprintf(out, " %u", (uint8)(i->descriptor[element] * 255.0));
      }
      fprintf(out, "\n");
    }
    fclose(out);
  }

  inline void write_ip_record(std::ofstream &f, InterestPoint const& p) {
    f.write((char*)&(p.x), sizeof(p.x));
    f.write((char*)&(p.y), sizeof(p.y));
    f.write((char*)&(p.ix), sizeof(p.ix));
    f.write((char*)&(p.iy), sizeof(p.iy));
    f.write((char*)&(p.orientation), sizeof(p.orientation));
    f.write((char*)&(p.scale), sizeof(p.scale));
    f.write((char*)&(p.interest), sizeof(p.interest));
    f.write((char*)&(p.polarity), sizeof(p.polarity));
    f.write((char*)&(p.octave), sizeof(p.octave));
    f.write((char*)&(p.scale_lvl), sizeof(p.scale_lvl));
    uint64 size = p.size();
    f.write((char*)(&size), sizeof(uint64));
    for (size_t i = 0; i < p.descriptor.size(); ++i)
      f.write((char*)&(p.descriptor[i]), sizeof(p.descriptor[i]));
  }

  inline InterestPoint read_ip_record(std::ifstream &f) {
    if (!f)
      vw::vw_throw(vw::IOErr() << "Failed to read interest point from file.");

    InterestPoint ip;
    f.read((char*)&(ip.x), sizeof(ip.x));
    f.read((char*)&(ip.y), sizeof(ip.y));
    f.read((char*)&(ip.ix), sizeof(ip.ix));
    f.read((char*)&(ip.iy), sizeof(ip.iy));
    f.read((char*)&(ip.orientation), sizeof(ip.orientation));
    f.read((char*)&(ip.scale), sizeof(ip.scale));
    f.read((char*)&(ip.interest), sizeof(ip.interest));
    f.read((char*)&(ip.polarity), sizeof(ip.polarity));
    f.read((char*)&(ip.octave), sizeof(ip.octave));
    f.read((char*)&(ip.scale_lvl), sizeof(ip.scale_lvl));

    uint64 size = 0; // Must initialize to avoid undefined behavior if reading failed
    f.read((char*)&(size), sizeof(uint64));
    if (!f)
      return ip; // Nothing to read

    ip.descriptor = Vector<double>(size);
    for (size_t i = 0; i < size; ++i)
      f.read((char*)&(ip.descriptor[i]), sizeof(ip.descriptor[i]));
    return ip;
  }

  void write_binary_ip_file(std::string ip_file, InterestPointList ip) {
    vw::create_out_dir(ip_file);

    std::ofstream f;
    f.open(ip_file.c_str(), std::ios::binary | std::ios::out);
    InterestPointList::iterator iter = ip.begin();
    uint64 size = ip.size();
    f.write((char*)&size, sizeof(uint64));
    for (; iter != ip.end(); ++iter)
      write_ip_record(f, *iter);
    f.close();
  }

  std::vector<InterestPoint> read_binary_ip_file(std::string ip_file) {
    std::vector<InterestPoint> result;

    std::ifstream f;
    f.open(ip_file.c_str(), std::ios::binary | std::ios::in);
    if (!f.is_open())
      vw_throw(IOErr() << "Failed to open \"" << ip_file << "\" as VWIP file.");

    uint64 size = 0;
    f.read((char*)&size, sizeof(uint64));
    if (!f)
      return result;

    for (size_t i = 0; i < size; ++i)
      result.push_back(read_ip_record(f));
    f.close();
    return result;
  }

  InterestPointList read_binary_ip_file_list(std::string ip_file) {
    InterestPointList result;

    std::ifstream f;
    f.open(ip_file.c_str(), std::ios::binary | std::ios::in);
    if (!f.is_open())
      vw_throw(IOErr() << "Failed to open \"" << ip_file << "\" as VWIP file.");

    uint64 size = 0;
    f.read((char*)&size, sizeof(uint64));
    if (!f)
      return result;

    for (size_t i = 0; i < size; ++i)
      result.push_back(read_ip_record(f));
    f.close();
    return result;
  }

  // Routines for reading & writing interest point match files
  void write_binary_match_file(std::string match_file,
                               std::vector<InterestPoint> const& ip1, std::vector<InterestPoint> const& ip2) {

    vw::create_out_dir(match_file);

    std::ofstream f;
    f.open(match_file.c_str(), std::ios::binary | std::ios::out);
    std::vector<InterestPoint>::const_iterator iter1 = ip1.begin();
    std::vector<InterestPoint>::const_iterator iter2 = ip2.begin();
    uint64 size1 = ip1.size();
    uint64 size2 = ip2.size();
    if (size1 != size2)
      vw_throw(IOErr()
               << "The vectors of matching interest points must have the same size.\n");

    f.write((char*)&size1, sizeof(uint64));
    f.write((char*)&size2, sizeof(uint64));
    for (; iter1 != ip1.end(); ++iter1)
      write_ip_record(f, *iter1);
    for (; iter2 != ip2.end(); ++iter2)
      write_ip_record(f, *iter2);
    f.close();
  }

  void read_binary_match_file(std::string match_file,
                              std::vector<InterestPoint> &ip1,
                              std::vector<InterestPoint> &ip2) {
    ip1.clear();
    ip2.clear();

    std::ifstream f;
    f.open(match_file.c_str(), std::ios::binary | std::ios::in);

    // Allow match files to not exist. That because we do not write empty match files.
    if (!f.is_open())
       return;

    // But if the file exists, we must be able to read from it
    uint64 size1 = 0, size2 = 0;
    if (f)
     f.read((char*)&size1, sizeof(uint64));
    else
      vw::vw_throw(vw::IOErr() << "Failed to read match file: " << match_file);
    if (f)
      f.read((char*)&size2, sizeof(uint64));
    else
      vw::vw_throw(vw::IOErr() << "Failed to read match file: " << match_file);

    if (size1 != size2)
      vw_throw(IOErr()
               << "The vectors of matching interest points must have the same size.\n");

    if (!f) {
      // This is a bugfix for a crash. Apparently junk was being read.
      return;
    }

    for (size_t i = 0; i < size1; ++i)
      ip1.push_back(read_ip_record(f));
    for (size_t i = 0; i < size2; ++i)
      ip2.push_back(read_ip_record(f));
    f.close();
  }

}} // namespace vw::ip
