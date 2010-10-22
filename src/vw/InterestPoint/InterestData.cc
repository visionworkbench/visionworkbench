// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file InterestData.cc
///
/// Basic classes and structures for storing image interest points.
///
#include <fstream>
#include <vw/InterestPoint/InterestData.h>

namespace vw {
namespace ip {

  void write_lowe_ascii_ip_file(std::string ip_file, InterestPointList ip) {

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

    uint64 size;
    f.read((char*)&(size), sizeof(uint64));
    ip.descriptor = Vector<double>(size);
    for (size_t i = 0; i < size; ++i)
      f.read((char*)&(ip.descriptor[i]), sizeof(ip.descriptor[i]));
    return ip;
  }

  void write_binary_ip_file(std::string ip_file, InterestPointList ip) {
    std::ofstream f;
    f.open(ip_file.c_str(), std::ios::binary | std::ios::out);
    InterestPointList::iterator iter = ip.begin();
    uint64 size = ip.size();
    f.write((char*)&size, sizeof(uint64));
    for ( ; iter != ip.end(); ++iter)
      write_ip_record(f, *iter);
    f.close();
  }

  std::vector<InterestPoint> read_binary_ip_file(std::string ip_file) {
    std::vector<InterestPoint> result;

    std::ifstream f;
    f.open(ip_file.c_str(), std::ios::binary | std::ios::in);

    // Error Handling
    if ( !f.is_open() )
      vw_throw( IOErr() << "Failed to open \"" << ip_file << "\" as VWIP file." );

    uint64 size;
    f.read((char*)&size, sizeof(uint64));
    for (size_t i = 0; i < size; ++i)
      result.push_back( read_ip_record(f) );
    f.close();
    return result;
  }

  // Routines for reading & writing interest point match files
  void write_binary_match_file(std::string match_file, std::vector<InterestPoint> const& ip1, std::vector<InterestPoint> const& ip2) {
    std::ofstream f;
    f.open(match_file.c_str(), std::ios::binary | std::ios::out);
    std::vector<InterestPoint>::const_iterator iter1 = ip1.begin();
    std::vector<InterestPoint>::const_iterator iter2 = ip2.begin();
    uint64 size1 = ip1.size();
    uint64 size2 = ip2.size();
    f.write((char*)&size1, sizeof(uint64));
    f.write((char*)&size2, sizeof(uint64));
    for ( ; iter1 != ip1.end(); ++iter1)
      write_ip_record(f, *iter1);
    for ( ; iter2 != ip2.end(); ++iter2)
      write_ip_record(f, *iter2);
    f.close();
  }

  void read_binary_match_file(std::string match_file, std::vector<InterestPoint> &ip1, std::vector<InterestPoint> &ip2) {
    ip1.clear();
    ip2.clear();

    std::ifstream f;
    f.open(match_file.c_str(), std::ios::binary | std::ios::in);

    // Error Handling
    if ( !f.is_open() )
      vw_throw( IOErr() << "Failed to open \"" << match_file << "\" as Match file." );

    uint64 size1, size2;
    f.read((char*)&size1, sizeof(uint64));
    f.read((char*)&size2, sizeof(uint64));
    for (size_t i = 0; i < size1; ++i)
      ip1.push_back( read_ip_record(f) );
    for (size_t i = 0; i < size2; ++i)
      ip2.push_back( read_ip_record(f) );
    f.close();
  }

  std::vector<Vector3> iplist_to_vectorlist(std::vector<InterestPoint> const& iplist) {
    std::vector<Vector3> result(iplist.size());
    for (size_t i=0; i < iplist.size(); ++i) {
      result[i][0] = iplist[i].x;
      result[i][1] = iplist[i].y;
      result[i][2] = 1; // homogeneous vector
    }
    return result;
  }

  // You are highly discouraged in using this as all descriptor
  // information is lost.
  std::vector<InterestPoint> vectorlist_to_iplist(std::vector<Vector3> const& veclist) {
    std::vector<InterestPoint> result(veclist.size());
    for (size_t i=0; i < veclist.size(); ++i) {
      result[i].x = veclist[i][0];
      result[i].y = veclist[i][1];
    }
    return result;
  }

  /// Helpful functors
  void remove_descriptor( InterestPoint & ip ) { ip.descriptor.set_size(0); }

}} // namespace vw::ip
