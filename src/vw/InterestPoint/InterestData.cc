// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__

/// \file InterestData.cc
/// 
/// Basic classes and structures for storing image interest points.
/// 
#include <vw/InterestPoint/InterestData.h>

namespace vw {
namespace ip {

  void write_lowe_ascii_ip_file(std::string ip_file, InterestPointList ip) {

    unsigned num_pts = ip.size();
    if (num_pts == 0) 
      vw_throw(IOErr() << "Attempted to write Lowe SIFT format interest point file with an empty list of interest points.");

    unsigned size = ip.front().descriptor.size();
    
    // Write out detected interest points to file.
    FILE *out = fopen(ip_file.c_str(), "w");
    fprintf(out, "%u %u\n", num_pts, size);
    for (InterestPointList::iterator i = ip.begin(); i != ip.end(); ++i) {
      float orientation = i->orientation;
      while (orientation > M_PI) orientation -= 2 * M_PI;
      while (orientation < -M_PI) orientation += 2 * M_PI;
      fprintf(out, "%.2f %.2f %.2f %.3f", i->y, i->x, i->scale, orientation);
      for (unsigned element = 0; element < size; ++element) {
        if (element % 20 == 0) fprintf(out, "\n");
        fprintf(out, " %u", (unsigned)(i->descriptor[element] * 255.0));
      }
      fprintf(out, "\n");
    }
    fclose(out);
  }

  void write_binary_ip_file(std::string ip_file, InterestPointList ip) {
    std::ofstream f;
    f.open(ip_file.c_str(), std::ios::binary | std::ios::out);
    InterestPointList::iterator iter = ip.begin();
    int size = ip.size();
    f.write((char*)&size, sizeof(int));
    for ( ; iter != ip.end(); ++iter) {
      f.write((char*)&(iter->x), sizeof(iter->x));
      f.write((char*)&(iter->y), sizeof(iter->y));
      f.write((char*)&(iter->ix), sizeof(iter->ix));
      f.write((char*)&(iter->iy), sizeof(iter->iy));
      f.write((char*)&(iter->orientation), sizeof(iter->orientation));
      f.write((char*)&(iter->scale), sizeof(iter->scale));
      f.write((char*)&(iter->interest), sizeof(iter->interest));
      int size = iter->size();
      f.write((char*)(&size), sizeof(int));
      for (unsigned i = 0; i < iter->descriptor.size(); ++i) 
        f.write((char*)&(iter->descriptor[i]), sizeof(iter->descriptor[i]));
    }
    f.close();
  }

  std::vector<InterestPoint> read_binary_ip_file(std::string ip_file) {
    std::vector<InterestPoint> result;
    
    std::ifstream f;
    f.open(ip_file.c_str(), std::ios::binary | std::ios::in);
    int size;
    f.read((char*)&size, sizeof(int));
    for (int i = 0; i < size; ++i) {
      InterestPoint ip;
      f.read((char*)&(ip.x), sizeof(ip.x));
      f.read((char*)&(ip.y), sizeof(ip.y));
      f.read((char*)&(ip.ix), sizeof(ip.ix));
      f.read((char*)&(ip.iy), sizeof(ip.iy));
      f.read((char*)&(ip.orientation), sizeof(ip.orientation));
      f.read((char*)&(ip.scale), sizeof(ip.scale));
      f.read((char*)&(ip.interest), sizeof(ip.interest));

      int size;
      f.read((char*)&(size), sizeof(ip.descriptor.size()));
      ip.descriptor = Vector<double>(size);
      for (int i = 0; i < size; ++i) 
        f.read((char*)&(ip.descriptor[i]), sizeof(ip.descriptor[i]));
      result.push_back(ip);
    }
    f.close();
    return result;
  }


  
}} // namespace vw::ip
