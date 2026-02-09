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


/// \file Matcher.h
///
/// Classes and functions for matching image interest points.
///
#ifndef _INTERESTPOINT_MATCHER_IO_H_
#define _INTERESTPOINT_MATCHER_IO_H_

#include <vw/InterestPoint/InterestPoint.h>

#include <string>
#include <vector>


namespace vw {
namespace ip {

  /// Return the basename with no extension, and shorten it if need be.
  /// If filename starts with out_prefix followed by dash, strip those.
  std::string strip_path(std::string out_prefix, std::string filename);

  /// The name of the match file.
  std::string match_filename(std::string const& out_prefix,
                             std::string const& input_file1,
                             std::string const& input_file2,
                             bool plain_text = false);

  /// Convert match file name to clean match file name.
  std::string clean_match_filename(std::string const& match_filename,
                                   bool plain_text = false);

  /// The name of the clean match file.
  std::string clean_match_filename(std::string const& out_prefix,
                                   std::string const& input_file1,
                                   std::string const& input_file2,
                                   bool plain_text = false);
  
  /// The name of an IP file.
  std::string ip_filename(std::string const& out_prefix,
                          std::string const& input_file);

  // Routines for reading & writing interest point data files
  void write_lowe_ascii_ip_file(std::string ip_file, InterestPointList ip);
  void write_binary_ip_file    (std::string ip_file, InterestPointList ip);
  std::vector<InterestPoint> read_binary_ip_file     (std::string ip_file);
  InterestPointList          read_binary_ip_file_list(std::string ip_file);

  // Routines for reading & writing interest point match files
  void write_binary_match_file(std::string match_file, std::vector<InterestPoint> const& ip1,
                               std::vector<InterestPoint> const& ip2);
  void write_text_match_file(std::string match_file, std::vector<InterestPoint> const& ip1,
                             std::vector<InterestPoint> const& ip2);
  void read_binary_match_file(std::string match_file, std::vector<InterestPoint> &ip1,
                              std::vector<InterestPoint> &ip2);
  void read_text_match_file(std::string match_file, std::vector<InterestPoint> &ip1,
                            std::vector<InterestPoint> &ip2);

}} // namespace vw::ip

#endif // _INTEREST_POINT_MATCHER_IO_H_
