// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
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

/// \file FileIO/FileMetadata.cc
/// 
/// Provides support for image file metadata.
///

#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <vw/config.h>
#include <iostream>

#include <vw/FileIO/FileMetadata.h>
#include <vw/FileIO/DiskImageResource.h>

#include <boost/algorithm/string.hpp>

namespace vw {
  
  FileMetadata* FileMetadataCollection::file_metadata(bool* is_readable /*= 0*/, FileMetadataCollectionIterator* pos /*= 0*/) const {
    std::list<std::pair<FileMetadata*, const FileMetadata*> >::const_iterator i;
    FileMetadata* m = 0;
    
    if(pos && pos->initialized)
      i = pos->i;
    else
      i = metas.begin();
      
    for(; i != metas.end(); i++) {
      if((m = (*i).first)) {
        if(is_readable)
          *is_readable = true;
        i++;
        break;
      }
    }
      
    if(pos)
    {
      pos->i = i;
      pos->initialized = true;
    }
    
    return m;
  }

  const FileMetadata* FileMetadataCollection::file_metadata_const(bool* is_readable /*= 0*/, FileMetadataCollectionIterator* pos /*= 0*/) const {
    std::list<std::pair<FileMetadata*, const FileMetadata*> >::const_iterator i;
    const FileMetadata* m = 0;
    
    if(pos && pos->initialized)
      i = pos->i;
    else
      i = metas.begin();
      
    for(; i != metas.end(); i++) {
      if((m = (*i).second)) {
        if(is_readable)
          *is_readable = ((*i).first != 0);
        i++;
        break;
      }
    }
      
    if(pos)
    {
      pos->i = i;
      pos->initialized = true;
    }
    
    return m;
  }

  bool FileMetadataCollection::contains_file_metadata(const std::string& metadata_type, bool require_readable /*= false*/) const {
    std::list<std::pair<FileMetadata*, const FileMetadata*> >::const_iterator i;
    for(i = metas.begin(); i != metas.end(); i++) {
      if(require_readable) {
        if((*i).first && (*i).first->metadata_type() == metadata_type)
          return true;
      }
      else {
        if((*i).second && (*i).second->metadata_type() == metadata_type)
          return true;
      }
    }
    return false;
  }

  void FileMetadataCollection::read_file_metadata(DiskImageResource* r) {
    std::list<std::pair<FileMetadata*, const FileMetadata*> >::iterator i;
    //std::cout << "Reading metadata from DiskImageResource" << r->type() << "." << std::endl;
    for(i = metas.begin(); i != metas.end(); i++) {
      if((*i).first) {
        if(vw::DiskImageResource::supports_metadata_type(r->type(), (*i).first->metadata_type())) {
          //std::cout << "Reading from DiskImageResource" << r->type() << ": " << (*i).first->metadata_type() << "." << std::endl;
          (*i).first->read_file_metadata(r);
        }
        else
          std::cout << "Best DiskImageResource ('" << r->type() << "') does not support metadata type '" << (*i).first->metadata_type() << "', skipping." << std::endl;
      }
    }
  }

  void FileMetadataCollection::write_file_metadata(DiskImageResource* r) const {
    std::list<std::pair<FileMetadata*, const FileMetadata*> >::const_iterator i;
    //std::cout << "Writing metadata to DiskImageResource" << r->type() << "." << std::endl;
    for(i = metas.begin(); i != metas.end(); i++) {
      if((*i).second) {
        if(vw::DiskImageResource::supports_metadata_type(r->type(), (*i).second->metadata_type())) {
          //std::cout << "Writing to DiskImageResource" << r->type() << ": " << (*i).second->metadata_type() << "." << std::endl;
          (*i).second->write_file_metadata(r);
        }
        else
          std::cout << "Best DiskImageResource ('" << r->type() << "') does not support metadata type '" << (*i).second->metadata_type() << "', skipping." << std::endl;
      }
    }
  }

} // namespace vw
