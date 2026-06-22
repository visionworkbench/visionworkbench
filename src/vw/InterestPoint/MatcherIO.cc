// __BEGIN_LICENSE__
//  Copyright (c) 2006-2026, United States Government as represented by the
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
#include <vw/InterestPoint/InterestPoint.h>

#include <vw/FileIO/FileUtils.h>

#include <boost/filesystem/operations.hpp>
#include <cstdint>
#include <iomanip>
#include <sstream>

namespace fs = boost::filesystem;

namespace vw { namespace ip {

// The file system limit on the length of a single path component is 255 bytes
// on most file systems. Budget the generated names against this, with a margin.
const size_t g_name_budget = 240;

// A fixed, portable 64-bit FNV-1a hash, as 16 hex digits. This keeps shortened
// file names distinct across platforms and runs, unlike std::hash, which
// differs between standard library implementations.
static std::string fnv1a64Hex(std::string const& s) {
  uint64_t h = 14695981039346656037ULL;
  for (unsigned char c: s) {
    h ^= static_cast<uint64_t>(c);
    h *= 1099511628211ULL;
  }
  std::ostringstream os;
  os << std::hex << std::setw(16) << std::setfill('0') << h;
  return os.str();
}

// Shorten a name to at most max_len characters. If it already fits, it is
// returned unchanged. Otherwise it is reduced to a leading portion plus a
// 64-bit hash of the full name, so distinct names stay distinct. See the
// documentation on match file naming.
static std::string shortenName(std::string const& name, size_t max_len) {
  if (name.size() <= max_len)
    return name;
  std::string hash = fnv1a64Hex(name); // 16 chars, keeps names distinct
  // The callers cap the output prefix so that a full hash always fits (see
  // match_filename and ip_filename). If there is no room for a leading part and
  // an underscore, return just the hash, so it is not itself truncated.
  if (max_len <= hash.size() + 1)
    return hash.substr(0, (max_len < hash.size() ? max_len : hash.size()));
  // Reserve room for an underscore and the hash.
  size_t lead = max_len - hash.size() - 1;
  return name.substr(0, lead) + "_" + hash;
}

// The .vwip file name for a standalone image, with no output prefix. The input
// is a path with the extension already removed. The directory is kept, and the
// file name is shortened if needed to fit the file system limit. This is used
// by ipfind and ipmatch, and shares the shortening logic with the prefixed
// file names. See the documentation on match file naming.
std::string shorten_vwip_name(std::string const& path_no_ext) {
  fs::path p(path_no_ext);
  size_t cap = (g_name_budget > 5) ? (g_name_budget - 5) : 1; // ".vwip"
  std::string name = shortenName(p.filename().string(), cap) + ".vwip";
  if (p.has_parent_path())
    return (p.parent_path() / name).string();
  return name;
}

/// Return the basename with no extension. If filename starts with out_prefix
/// followed by dash, strip those.
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

// If the basename of out_prefix is so long that it would not leave room for a
// full hash on each shortened name, truncate that basename, keeping the
// directory. This moves the shortening from the per-pair names to the shared
// prefix, so the names keep full, distinct hashes. Only triggers for an
// extremely long output prefix. See the documentation on match file naming.
static std::string capPrefixBase(std::string const& out_prefix, size_t reserve) {
  if (out_prefix == "")
    return out_prefix;
  fs::path p(out_prefix);
  std::string base = p.filename().string();
  size_t max_base = (g_name_budget > reserve) ? (g_name_budget - reserve) : 0;
  if (base.size() <= max_base)
    return out_prefix;
  base = base.substr(0, max_base);
  fs::path dir = p.parent_path();
  return dir.empty() ? base : (dir / base).string();
}

std::string match_filename(std::string const& out_prefix,
                           std::string const& input_file1,
                           std::string const& input_file2,
                           bool plain_text) {

  std::string name1 = strip_path(out_prefix, input_file1);
  std::string name2 = strip_path(out_prefix, input_file2);

  std::string ext = plain_text ? ".txt" : ".match";
  const size_t hash_len = 16;

  // The match file name must fit within the file system limit on the length of
  // a single path component. Only the part of the prefix after the last
  // directory separator counts toward the component length. If that prefix is
  // extremely long, cap it so a full hash still fits on each name, reserving
  // room for two hashes and the separators. See the documentation on match
  // file naming.
  std::string prefix = capPrefixBase(out_prefix, 1 + 2 + ext.size() + 2 * hash_len);

  // Budget the two image names against the limit. Shorten them with a hash if
  // they would not fit, in a way that keeps them distinct.
  std::string prefix_base = fs::path(prefix).filename().string();
  size_t overhead = prefix_base.size() + 2 + ext.size(); // "__" and extension
  if (prefix != "")
    overhead += 1; // the dash after the prefix
  size_t cap = (overhead + 2 < g_name_budget) ? (g_name_budget - overhead) / 2 : hash_len;

  if (name1.size() > cap || name2.size() > cap)
    vw_out(WarningMessage) << "Shortening long file names to fit the file "
      << "system limit, for: " << input_file1 << " and " << input_file2 << ".\n";

  name1 = shortenName(name1, cap);
  name2 = shortenName(name2, cap);

  std::string suffix = name1 + "__" + name2 + ext;

  if (prefix == "")
    return suffix;

  return prefix + "-" + suffix;
}

/// Convert match file name to clean match file name.
std::string clean_match_filename(std::string const& match_file,
                                   bool plain_text) {
  std::string clean_match_file = fs::path(match_file).replace_extension("").string();
  if (plain_text)
    clean_match_file += "-clean.txt";
  else
    clean_match_file += "-clean.match";
  return clean_match_file;
}

/// The name of the clean match file.
std::string clean_match_filename(std::string const& out_prefix,
                                 std::string const& input_file1,
                                 std::string const& input_file2,
                                 bool plain_text) {

  return clean_match_filename(match_filename(out_prefix, input_file1, input_file2, plain_text),
                              plain_text);
}

std::string ip_filename(std::string const& out_prefix,
                        std::string const& input_file) {
  std::string name = strip_path(out_prefix, input_file);

  // Shorten the name if, combined with the prefix, it would exceed the file
  // system limit on the length of a single path component. This single name
  // (one image) is budgeted on its own, independently of the match file (two
  // images). If the prefix is extremely long, cap it so a full hash still fits.
  // See the documentation on match file naming.
  const size_t hash_len = 16;
  std::string prefix = capPrefixBase(out_prefix, 1 + 5 + hash_len); // dash, ".vwip", hash
  std::string prefix_base = fs::path(prefix).filename().string();
  size_t overhead = prefix_base.size() + 1 + 5; // dash and ".vwip"
  size_t cap = (overhead < g_name_budget) ? (g_name_budget - overhead) : hash_len;
  name = shortenName(name, cap);

  return prefix + "-" + name + ".vwip";
}

void ip_filenames(std::string const& out_prefix,
                  std::string const& input_file1,
                  std::string const& input_file2,
                  std::string & output_ip1,
                  std::string & output_ip2) {
  output_ip1 = ip_filename(out_prefix, input_file1);
  output_ip2 = ip_filename(out_prefix, input_file2);
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
  for (InterestPointList::iterator i = ip.begin(); i != ip.end(); i++) {
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
  for (size_t i = 0; i < p.descriptor.size(); i++)
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
  for (size_t i = 0; i < size; i++)
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
  for (; iter != ip.end(); iter++)
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

  for (size_t i = 0; i < size; i++)
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

  for (size_t i = 0; i < size; i++)
    result.push_back(read_ip_record(f));
  f.close();
  return result;
}

// Routines for reading & writing interest point match files
void write_binary_match_file(std::string match_file,
                             std::vector<InterestPoint> const& ip1,
                             std::vector<InterestPoint> const& ip2) {
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
  for (; iter1 != ip1.end(); iter1++)
    write_ip_record(f, *iter1);
  for (; iter2 != ip2.end(); iter2++)
    write_ip_record(f, *iter2);
  f.close();
}

// Check if file ends with .txt (case insensitive)
bool hasTxtExtension(std::string const& filename) {
  if (filename.size() < 4)
    return false;
  std::string ext = filename.substr(filename.size() - 4);
  for (size_t i = 0; i < ext.size(); i++)
    ext[i] = std::tolower(ext[i]);
  return (ext == ".txt");
}

// Write a text file with the interest point matches. Use float precision.
void write_text_match_file(std::string match_file,
                           std::vector<InterestPoint> const& ip1,
                           std::vector<InterestPoint> const& ip2) {
  if (!hasTxtExtension(match_file))
    vw_throw(IOErr() << "Text match file must have .txt extension: " << match_file);

  vw::create_out_dir(match_file);

  size_t size1 = ip1.size();
  size_t size2 = ip2.size();
  if (size1 != size2)
    vw_throw(IOErr()
              << "The vectors of matching interest points must have the same size.\n");

  std::ofstream f(match_file.c_str());
  if (!f.is_open())
    vw_throw(IOErr() << "Failed to open \"" << match_file << "\" for writing.");

  // Must use 9 digits to avoid losing precision when writing float values
  f.precision(9);
  for (size_t i = 0; i < size1; i++) {
    f << ip1[i].x << " " << ip1[i].y << " " << ip1[i].scale << " "
      << ip2[i].x << " " << ip2[i].y << " " << ip2[i].scale << "\n";
  }
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

  for (size_t i = 0; i < size1; i++)
    ip1.push_back(read_ip_record(f));
  for (size_t i = 0; i < size2; i++)
    ip2.push_back(read_ip_record(f));
  f.close();
}

// Read a text file with the interest point matches.
void read_text_match_file(std::string match_file,
                          std::vector<InterestPoint> &ip1,
                          std::vector<InterestPoint> &ip2) {
  ip1.clear();
  ip2.clear();

  if (!hasTxtExtension(match_file))
    vw_throw(IOErr() << "Text match file must have .txt extension: " << match_file);

  std::ifstream f(match_file.c_str());
  if (!f.is_open())
    vw_throw(IOErr() << "Failed to open " << match_file << " for reading.");

  std::string line;
  int line_num = 0;
  while (std::getline(f, line)) {
    line_num++;

    // Skip lines with only whitespace
    bool only_spaces = true;
    for (size_t i = 0; i < line.size(); i++) {
      if (!std::isspace(line[i])) {
        only_spaces = false;
        break;
      }
    }
    if (only_spaces)
      continue;

    float x1, y1, unc1, x2, y2, unc2;
    char extra;
    std::istringstream iss(line);

    if (!(iss >> x1 >> y1 >> unc1 >> x2 >> y2 >> unc2))
      vw_throw(IOErr() << "Failed to parse 6 values at line " << line_num
               << " in match file: " << match_file);

    // Check for extra values
    if (iss >> extra)
      vw_out(WarningMessage) << "Warning: More than 6 values at line " << line_num
                             << " in match file: " << match_file << "\n";

    // Such a check will throw for non-positive and for nan values
    if (!(unc1 > 0))
      vw_throw(IOErr() << "Uncertainty unc1 must be positive at line " << line_num
               << " in match file: " << match_file);
    if (!(unc2 > 0))
      vw_throw(IOErr() << "Uncertainty unc2 must be positive at line " << line_num
               << " in match file: " << match_file);

    InterestPoint p1, p2;
    p1.x = x1;
    p1.y = y1;
    p1.ix = round(x1);
    p1.iy = round(y1);
    p1.scale = unc1;
    p2.x = x2;
    p2.y = y2;
    p2.ix = round(x2);
    p2.iy = round(y2);
    p2.scale = unc2;

    ip1.push_back(p1);
    ip2.push_back(p2);
  }

  f.close();
}

// Wrapper that dispatches to text or binary write based on plain_text flag
void write_match_file(std::string match_file, std::vector<InterestPoint> const& ip1,
                      std::vector<InterestPoint> const& ip2, bool plain_text) {
  if (plain_text)
    write_text_match_file(match_file, ip1, ip2);
  else
    write_binary_match_file(match_file, ip1, ip2);
}

// Wrapper that dispatches to text or binary read based on plain_text flag
void read_match_file(std::string match_file, std::vector<InterestPoint> &ip1,
                     std::vector<InterestPoint> &ip2, bool plain_text) {
  if (plain_text)
    read_text_match_file(match_file, ip1, ip2);
  else
    read_binary_match_file(match_file, ip1, ip2);
}

}} // namespace vw::ip
