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
#ifndef __VW_FILEIO_FILE_METADATA_H__
#define __VW_FILEIO_FILE_METADATA_H__

#include <list>

// Boost
#include <boost/algorithm/string.hpp>
 
namespace vw {

  /// Forward declarations.
  class DiskImageResource;
  template<class KeyT,class DataT,class PropertyT> class PropertyMultiMap;

  /// All file metadata types (GeoReference, CameraMetadata, etc) inherit from this class.
  class FileMetadata
  {
  public:
    /// Return the type of metadata. Subclasses should implement this by
    /// by calling a static function metadata_type_static().
    virtual std::string metadata_type(void) const = 0;
    /// Read metadata from the given DiskImageResource.
    virtual void read_file_metadata(DiskImageResource* r) = 0;
    /// Write metadata to the given DiskImageResource.
    virtual void write_file_metadata(DiskImageResource* r) const = 0;

    /// Useful types for all child classes of FileMetadata.
    typedef void (*disk_image_resource_register_func)(std::string const& metadata_type);
    typedef void (*read_metadata_func)(FileMetadata* fmeta, DiskImageResource* r);
    typedef void (*write_metadata_func)(DiskImageResource* r, FileMetadata const* fmeta);
    typedef PropertyMultiMap<std::string,read_metadata_func,std::string> ReadMetadataMapType;
    typedef PropertyMultiMap<std::string,write_metadata_func,std::string> WriteMetadataMapType;
    
    /// Constructor.
    FileMetadata() {}
    /// Destructor.
    virtual ~FileMetadata() {};
  };
  
  
  /// The FileMetadataCollection class allows image IO classes to
  /// read and write GeoReference, Camera, etc, metadata.
  class FileMetadataCollection {
  public:
    /// The FileMetadataCollectionIterator class stores the current position
    /// in a FileMetadataCollection when calling FileMetadataCollection::file_metadata()
    /// or FileMetadataCollection::file_metadata_const().
    class FileMetadataCollectionIterator
    {
    public:
      std::list<std::pair<FileMetadata*, const FileMetadata*> >::const_iterator i;
      bool initialized;
      
      void reset(void)
      {
        initialized = false;
      }
      
      FileMetadataCollectionIterator() : initialized(false)
      {
      }
      
      ~FileMetadataCollectionIterator()
      {
      }
    };
  
  protected:
    std::list<std::pair<FileMetadata*, const FileMetadata*> > metas;

  public:
    /// Create an empty FileMetadataCollection.
    static FileMetadataCollection create(void) {
      FileMetadataCollection fmeta;
      return fmeta;
    }
    
    /// Return the next registered non-const FileMetadata* in the collection,
    /// and optionally whether it is readable.
    FileMetadata* file_metadata(bool* is_readable = 0, FileMetadataCollectionIterator* pos = 0) const;
    
    /// Return the next registered const or non-const FileMetadata* in the collection,
    /// and optionally whether it is readable.
    const FileMetadata* file_metadata_const(bool* is_readable = 0, FileMetadataCollectionIterator* pos = 0) const;
    
    /// Return whether the collection includes a specific metadata type, optionally
    /// only if it is readable.
    bool contains_file_metadata(const std::string& metadata_type, bool require_readable = false) const;
    
    /// Associate a FileMetadata* with the collection.
    void associate_file_metadata(FileMetadata* m) {
      metas.push_back(std::make_pair(m, m));
    }
    
    /// Associate a const FileMetadata* with the collection.
    void associate_file_metadata_const(const FileMetadata* m) {
      metas.push_back(std::make_pair((FileMetadata*)0, m));
    }
  
    /// Read all supported metadata from the given DiskImageResource.
    void read_file_metadata(DiskImageResource* r);
  
    /// Write all supported metadata to the given DiskImageResource.
    void write_file_metadata(DiskImageResource* r) const;
  
    /// Construct a default FileMetadataCollection. This FileMetadataCollection
    /// does not include any metadata (Cartography, Camera, etc).
    FileMetadataCollection() {
    }
    
    /// Destroy FileMetadataCollection.
    ~FileMetadataCollection() {
    }
    
  };
  
  inline std::ostream& operator<<(std::ostream& os, const FileMetadataCollection& fmeta) {
    const vw::FileMetadata* m;
    bool is_readable = true;
    vw::FileMetadataCollection::FileMetadataCollectionIterator i;
    os << "-- File Metadata Collection Object --\n";
    while((m = fmeta.file_metadata_const(&is_readable, &i)) != 0)
      os << "\tContains " <<  m->metadata_type() << (is_readable ? "" : " (write-only)") << " Metadata\n";
    return os;
  }


} // namespace vw

#endif // __VW_FILEIO_FILE_METADATA_H__
