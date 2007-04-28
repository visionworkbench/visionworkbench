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

/// \file JP2.h
/// 
/// Provides support for the JPEG2000 image file formats.
///
#ifndef __VW_FILEIO_JPEG2000_H__
#define __VW_FILEIO_JPEG2000_H__

#include <list>
#include <iostream>

#include <arpa/inet.h>			   // ntohl, htonl
#include <endian.h>			   // __BYTE_ORDER

#include <boost/algorithm/string.hpp>

#include <vw/config.h>
//FIXME: probably need to include something so that vw::uint8, etc are defined
//#include <vw/FileIO/FileMetadata.h>
#include <vw/Math/Matrix.h>

#if __BYTE_ORDER == __BIG_ENDIAN
#define ntohll(a) (a)
#define htonll(a) (a)
#elif __BYTE_ORDER == __LITTLE_ENDIAN
//NOTE: "ULL" is for gcc; apparently msvc does not like it (in case this gets ported to Win32)
#define ntohll(a) ((((uint64)(a) & 0xff00000000000000ULL) >> 56) | \
                   (((uint64)(a) & 0x00ff000000000000ULL) >> 40) | \
                   (((uint64)(a) & 0x0000ff0000000000ULL) >> 24) | \
                   (((uint64)(a) & 0x000000ff00000000ULL) >>  8) | \
                   (((uint64)(a) & 0x00000000ff000000ULL) <<  8) | \
                   (((uint64)(a) & 0x0000000000ff0000ULL) << 24) | \
                   (((uint64)(a) & 0x000000000000ff00ULL) << 40) | \
                   (((uint64)(a) & 0x00000000000000ffULL) << 56))
#define htonll(a) ((((uint64)(a) & 0xff00000000000000ULL) >> 56) | \
                   (((uint64)(a) & 0x00ff000000000000ULL) >> 40) | \
                   (((uint64)(a) & 0x0000ff0000000000ULL) >> 24) | \
                   (((uint64)(a) & 0x000000ff00000000ULL) >>  8) | \
                   (((uint64)(a) & 0x00000000ff000000ULL) <<  8) | \
                   (((uint64)(a) & 0x0000000000ff0000ULL) << 24) | \
                   (((uint64)(a) & 0x000000000000ff00ULL) << 40) | \
                   (((uint64)(a) & 0x00000000000000ffULL) << 56))
#else
#error "__BYTE_ORDER must be either __LITTLE_ENDIAN or __BIG_ENDIAN."
#endif


namespace vw
{

  // These classes parse jp2 and jpx files, and upgrade jp2 -> jp2-compatible jpx
  // References: jp2 file format: http://www.jpeg.org/public/15444-1annexi.pdf
  //             jpx file format: http://www.jpeg.org/public/15444-2annexm.pdf
  
  // Forward declaration.
  class JP2Box;
  
  typedef std::list<JP2Box*> JP2BoxList;
  typedef std::list<std::pair<uint16,bool> > JP2ReaderRequirementsList;
  
  #if 0
  template <class MaskT>
  class JP2ReaderRequirements
  {
  protected:
    uint8* ML;
    uint16* NSF;
    uint16* NVF;
  
  public:
    MaskT* FUAM;
    MaskT* DCM;
    uint16** SF;
    MaskT** SM;
    uint16** VF;
    MaskT** VM;
    uint8* dat;
    uint64 bytes;
    bool free_dat;

    void print(void)
    {
      MaskT m;
      bool first, firstfirst;
      int i, j;

      std::cout << "standard features";
      firstfirst = true;
      for(i = 0, m = 1; i < (8 * (*ML)); i++, m <<= 1)
      {
        first = true;
        for(j = 0; j < *NSF; j++)
        {
          if(*(SM[i]) & m)
          {
            if(firstfirst)
            {
              std::cout << " (";
              firstfirst = false;
              first = false;
            }
            else if(first)
            {
              std::cout << " & (";
              first = false;
            }
            else
              std::cout << " |";
            std::cout << " " << *(SF[i]);
          }
        }
        if(!first)
          std::cout << ")";
      }

      std::cout << ", vendor features";
      firstfirst = true;
      for(i = 0, m = 1; i < (8 * (*ML)); i++, m <<= 1)
      {
        first = true;
        for(j = 0; j < *NVF; j++)
        {
          if(*(VM[i]) & m)
          {
            if(firstfirst)
            {
              std::cout << " (";
              firstfirst = false;
              first = false;
            }
            else if(first)
            {
              std::cout << " & (";
              first = false;
            }
            else
              std::cout << " |";
            std::cout << " " << *(VF[i]);
          }
        }
        if(!first)
          std::cout << ")";
      }
    }
    
    // Construct an empty Reader Requirements box with NSF_ standard features and NVF_ vendor features
    // 'free_dat_' is whether to free the data buffer when this object is destroyed
    JP2ReaderRequirements(uint8 NSF_, uint8 NVF_, bool free_dat_ = true) : free_dat(free_dat_)
    {
      uint8* d;
      uint8 i;
      
      bytes = sizeof(uint8) + (2 + NSF_ + NVF_) * sizeof(uint16) + (2 + NSF_ + NVF_) * sizeof(MaskT);
      dat = new uint8[bytes];
      if(NSF_ > 0)
      {
        SF = new uint16*[NSF_];
        SM = new MaskT*[NSF_];
      }
      else
      {
        SF = 0;
        SM = 0;
      }
      if(NVF_ > 0)
      {
        VF = new uint16*[NVF_];
        VM = new MaskT*[NVF_];
      }
      else
      {
        VF = 0;
        VM = 0;
      }
      d = dat;
      ML = d; *ML = sizeof(MaskT); d++;
      FUAM = (MaskT*)d; *FUAM = 0; d += sizeof(MaskT);
      DCM = (MaskT*)d; *DCM = 0; d += sizeof(MaskT);
      NSF = (uint16*)d; *NSF = NSF_; d += sizeof(uint16);
      for(i = 0; i < *NSF; i++)
      {
        SF[i] = (uint16*)d; *(SF[i]) = 0; d += sizeof(uint16);
        SM[i] = (MaskT*)d; *(SM[i]) = 0; d += sizeof(MaskT);
      }
      NVF = (uint16*)d; *NVF = NVF_; d += sizeof(uint16);
      for(i = 0; i < *NVF; i++)
      {
        VF[i] = (uint16*)d; *(VF[i]) = 0; d += sizeof(uint16);
        VM[i] = (MaskT*)d; *(VM[i]) = 0; d += sizeof(MaskT);
      }
    }

    // Construct a Reader Requirements box from data buffer dat_ of length bytes_
    // 'free_dat_' is whether to free the data buffer when this object is destroyed
    JP2ReaderRequirements(uint8* dat_, uint64 bytes_, bool free_dat_ = false) : dat(dat_), free_dat(free_dat_)
    {
      uint8* d = dat_;
      int i;

      ML = d; d++;
      FUAM = (MaskT*)d; d += sizeof(MaskT);
      DCM = (MaskT*)d; d += sizeof(MaskT);
      NSF = (uint16*)d; d += sizeof(uint16);
      for(i = 0; i < *NSF; i++)
      {
        SF[i] = (uint16*)d; d += sizeof(uint16);
        SM[i] = (MaskT*)d; d += sizeof(MaskT);
      }
      NVF = (uint16*)d; d += sizeof(uint16);
      for(i = 0; i < *NVF; i++)
      {
        VF[i] = (uint16*)d; d += sizeof(uint16);
        VM[i] = (MaskT*)d; d += sizeof(MaskT);
      }
    }
    
    // Destructor
    ~JP2ReaderRequirements()
    {
      if(SF)
        delete[] SF;
      if(SM)
        delete[] SM;
      if(VF)
        delete[] VF;
      if(VM)
        delete[] VM;
      if(free_dat)
        delete[] dat;
    }
  };
  #endif

  // Generic jp2/jpx box
  class JP2Box
  {
  protected:
    uint32 TBox;
  
  public:
    // Get box type
    uint32 box_type(void)
    {
      return TBox;
    }
    
    // Get number of bytes in the body of the box (DBox)
    virtual uint64 bytes_dbox(void) = 0;
    
    // Get number of bytes including the header
    virtual uint64 bytes(void) = 0;

    // Print with d intentations
    virtual void print(int d = 0) = 0;
    
    // Serialize into buffer d, which must be large enough
    // (call bytes() and allocate buffer first)
    virtual void serialize(uint8* d) = 0;

    // Construct a JP2Box with only a box type
    JP2Box(uint32 type) : TBox(type)
    {
    }

    // Destructor
    ~JP2Box()
    {
    }
  
  protected:
    // Determine whether box with type 'type' is a superbox
    static bool is_superbox(uint32 type)
    {
      bool retval = false;
    
      switch(type)
      {
      case 0x6A703268: // "jp2h"
      case 0x6A706368: // "jpch"
      case 0x6A706C68: // "jplh"
      case 0x63677270: // "cgrp"
      case 0x6674626C: // "ftbl"
      case 0x636F6D70: // "comp"
      case 0x61736F63: // "asoc"
      case 0x64726570: // "drep"
      case 0x00000000: // fake box type for root JP2File "box"
        retval = true;
        break;
      default:
        break;
      }
      
      return retval;
    }

    // Retrieve type from box data d (d includes the box header)
    static uint32 interp_type(uint8* d)
    {
      return ntohl(*((uint32*)&d[4]));
    }

    // Reverse the bytes in TBox_ if necessary
    static uint32 interp_type(uint32 TBox_, bool byte_order_converted)
    {
      return byte_order_converted ? TBox_ : ntohl(TBox_);
    }

    // Determine how many bytes (including the header) are in a box with the given LBox_, XLBox_,
    // and number of bytes remaining in the file. If necessary, reverses the bytes of LBox_ and XLBox_.
    //NOTE: nbytes is assumed to always have had its byte order coverted
    static uint64 interp_bytes(uint32 LBox_, uint64 XLBox_, uint64 nbytes, bool byte_order_converted = true)
    {
      uint32 LBox = byte_order_converted ? LBox_ : ntohl(LBox_);
      uint64 retval = 0;

      switch(LBox)
      {
      case 0:
        retval = nbytes;
        break;
      case 1:
        retval = byte_order_converted ? XLBox_ : ntohll(XLBox_);
        break;
      default:
        retval = LBox;
        break;
      }

      return retval;
    }

    // Determine how many bytes (including the header) are in the box d with nbytes remaining in the file
    // (d includes the box header)
    //NOTE: nbytes is assumed to always have had its byte order coverted
    static uint64 interp_bytes(uint8* d, uint64 nbytes)
    {
      uint32 LBox = ntohl(*((uint32*)d));
      uint64 retval = 0;

      switch(LBox)
      {
      case 0:
        retval = nbytes;
        break;
      case 1:
        retval = ntohll(*((uint64*)&d[8]));
        break;
      default:
        retval = LBox;
        break;
      }

      return retval;
    }

    // Determine how many header bytes to add to a box with b bytes of payload
    static uint64 header_bytes_to_add(uint64 b)
    {
      return ((b + 8) > 0x00000000ffffffffULL ? 16 : 8);
    }

    // Add the appropriate number of header bytes to a box with b bytes of payload
    static uint64 add_header_bytes(uint64 b)
    {
      return b + header_bytes_to_add(b);
    }

    // Determine how many header bytes to remove from a box with b total bytes (b includes the header)
    static uint64 header_bytes_to_remove(uint64 b)
    {
      return (b > 0x00000000ffffffffULL ? 16 : 8);
    }

    // Remove the appropriate number of header bytes from a box with b total bytes (b includes the header)
    static uint64 remove_header_bytes(uint64 b)
    {
      return b - header_bytes_to_remove(b);
    }

    // Print the box type string corresponding to numeric box type t
    static void print_type(uint32 t)
    {
      uint8* c;
      int j;
      
      t = htonl(t);
      c = (uint8*)&t;
      for(j = 0; j < 4; j++)
        std::cout << (char)c[j];
    }
    
    // Basic common print functionality for print()
    void print_basic(void)
    {
      print_type(TBox);
      std::cout << ": " << this->bytes() << "/" << this->bytes_dbox() << " bytes";
    }
    
    // Basic serialization functionality for serialize()
    virtual void serialize_basic(uint8** d_)
    {
      uint64 s = this->bytes();
      uint8* d = *d_;
    
      if(s > 0x00000000ffffffffULL)
      {
        *((uint32*)d) = htonl(1); d += sizeof(uint32);
        *((uint32*)d) = htonl(TBox); d += sizeof(uint32);
        *((uint32*)d) = htonll(s); d += sizeof(uint64);
      }
      else
      {
        *((uint32*)d) = htonl(s); d += sizeof(uint32);
        *((uint32*)d) = htonl(TBox); d += sizeof(uint32);
      }
      
      *d_ = d;
    }
  };

  // Reader Requirements box
  class JP2ReaderRequirementsBox : public JP2Box
  {
  protected:
    uint16 NSF;
    uint16 NVF;
    Matrix<uint16> requirements;

  public:

    // Get number of bytes in the body of the box (DBox)
    virtual uint64 bytes_dbox(void)
    {
      uint8 ML = mask_length(requirements.cols() - 2);
      return (NSF + NVF + 2) * (sizeof(uint16) + ML) + 1;
    }
    
    // Get number of bytes including the header
    virtual uint64 bytes(void)
    {
      return add_header_bytes(bytes_dbox());
    }

    // Print with d intentations
    virtual void print(int d = 0)
    {
      int terms;
      bool first, firstfirst;
      int i, j, k;

      for(j = 0; j < d; j++)
        std::cout << "  ";
      print_basic();

      std::cout << " (";
      for(k = 0; k <= 1; k++)
      {
        if(k == 0)
          std::cout << "fully understood:";
        else
          std::cout << ", display:";
        firstfirst = true;
        for(j = 2; j < requirements.cols(); j++)
        {
          if(!requirements(k, j))
            continue;
          terms = 0;
          for(i = 2; i < requirements.rows(); i++)
          {
            if(requirements(i, j))
              terms++;
          }
          first = true;
          for(i = 2; i < requirements.rows(); i++)
          {
            if(!requirements(i, j))
              continue;
            if(firstfirst)
            {
              std::cout << " ";
              if(terms > 1)
                std::cout << "(";
              firstfirst = false;
              first = false;
            }
            else if(first)
            {
              std::cout << " & ";
              if(terms > 1)
                std::cout << "(";
              first = false;
            }
            else
              std::cout << " | ";
            std::cout << requirements(i, 0) << (requirements(i, 1) ? "v" : "");
          }
          if(!first && terms > 1)
            std::cout << ")";
        }
      }
      std::cout << ")" << std::endl;
    }

    // Serialize into buffer d, which must be large enough
    // (call bytes() and allocate buffer first)
    virtual void serialize(uint8* d)
    {
      uint8 ML = mask_length(requirements.cols() - 2);
      int i;

      serialize_basic(&d);

      *d = ML; d++;
      serialize_mask(&d, ML, fuam());
      serialize_mask(&d, ML, dcm());
      *((uint16*)d) = htons(NSF); d += sizeof(uint16);
      for(i = 2; i < (NSF + NVF + 2); i++)
      {
        if(requirements(i, 1))
          continue;
        *((uint16*)d) = htons(requirements(i, 0)); d += sizeof(uint16);
        serialize_mask(&d, ML, requirements_row_(i));
      }
      *((uint16*)d) = htons(NVF); d += sizeof(uint16);
      for(i = 2; i < (NSF + NVF + 2); i++)
      {
        if(!requirements(i, 1))
          continue;
        *((uint16*)d) = htons(requirements(i, 0)); d += sizeof(uint16);
        serialize_mask(&d, ML, requirements_row_(i));
      }
    }

    // Number of rows in requirements table
    int requirements_rows(void)
    {
      return requirements.rows() - 2;
    }

    // Number of columns in requirements table
    int requirements_cols(void)
    {
      return requirements.cols() - 2;
    }

    // Add row to reader requirements table
    void add_requirements_row(uint16 req, uint64 mask = 0, bool vendor = false)
    {
      int row;

      requirements.set_size(requirements.rows() + 1, requirements.cols(), true);
      row = requirements.rows() - 1;
      requirements(row, 0) = req;
      if(vendor)
      {
        requirements(row, 1) = 1;
        NVF++;
      }
      else
      {
        requirements(row, 1) = 0;
        NSF++;
      }
      set_requirements_row_(row, mask);
    }

    // Add column to reader requirements table
    void add_requirements_col(uint64 mask = 0, bool in_dcm = true)
    {
      int col;

      requirements.set_size(requirements.rows(), requirements.cols() + 1, true);
      col = requirements.cols() - 1;
      requirements(0, col) = 1;
      requirements(1, col) = (in_dcm ? 1 : 0);
      set_requirements_col_(col, mask);
    }

    // Set row of reader requirements table to mask
    void set_requirements_row(int row, uint64 mask)
    {
      set_requirements_row_(row + 2, mask);
    }

    // Get row of reader requirements table
    uint64 requirements_row(int row)
    {
      return requirements_row_(row + 2);
    }

    // Get requirement number of row 'row'; optionally return whether it is vendor-specific
    uint16 requirement(int row, bool* vendor)
    {
      if(vendor)
        *vendor = requirements(row + 2, 1);
      return requirements(row + 2, 0);
    }

    // Set col of reader requirements table to mask
    void set_requirements_col(int col, uint64 mask)
    {
      set_requirements_col_(col + 2, mask);
    }

    // Get row of reader requirements table
    uint64 requirements_col(int col)
    {
      return requirements_col_(col + 2);
    }

    // Set FUAM of reader requirements table to mask
    void set_fuam(uint64 mask)
    {
      set_requirements_row_(0, mask);
    }

    // Get FUAM of reader requirements table
    uint64 fuam(void)
    {
      return requirements_row_(0);
    }

    // Set DCM of reader requirements table to mask
    void set_dcm(uint64 mask)
    {
      set_requirements_row_(1, mask);
    }

    // Get DCM of reader requirements table
    uint64 dcm(void)
    {
      return requirements_row_(1);
    }
    
    // Construct an empty Reader Requirements box
    JP2ReaderRequirementsBox() : JP2Box(0x72726571 /*"rreq"*/), NSF(0), NVF(0), requirements(2, 2)
    {
    }

    // Construct JP2ReaderRequirementsBox from buffer d containing a serialized box (d includes the header)
    // Buffer d is not deleted when this box is destroyed
    JP2ReaderRequirementsBox(uint8* d, uint64 nbytes) : JP2Box(0x72726571 /*"rreq"*/)
    {
      uint8 ML;
      uint64 FUAM;
      uint64 DCM;
      uint64 mask;
      int i;

      d += header_bytes_to_remove(interp_bytes(d, nbytes));

      ML = *d; d++;
      FUAM = interp_mask(d, ML); d += ML;
      DCM = interp_mask(d, ML); d += ML;
      NSF = ntohs(*((uint16*)d)); d += sizeof(uint16);
      requirements.set_size(NSF + 2, std::max(mask_used_bits(FUAM), mask_used_bits(DCM)) + 2);
      set_fuam(FUAM);
      set_dcm(DCM);
      for(i = 2; i < (NSF + 2); i++)
      {
        requirements(i, 0) = ntohs(*((uint16*)d)); d += sizeof(uint16);
        requirements(i, 1) = 0;
        mask = interp_mask(d, ML); d += ML;
        set_requirements_row_(i, mask);
      }
      NVF = ntohs(*((uint16*)d)); d += sizeof(uint16);
      requirements.set_size(NSF + NVF + 2, requirements.cols(), true);
      for(; i < (NSF + NVF + 2); i++)
      {
        requirements(i, 0) = ntohs(*((uint16*)d)); d += sizeof(uint16);
        requirements(i, 1) = 1;
        mask = interp_mask(d, ML); d += ML;
        set_requirements_row_(i, mask);
      }
    }
    
    // Destructor
    ~JP2ReaderRequirementsBox()
    {
    }
    
  protected:
  
    // Find (64 - number of leading zeros).
    static uint8 mask_used_bits(uint64 mask)
    {
      uint8 i;
      uint64 m;

      for(i = 64, m = 0x8000000000000000ULL; i > 0 && !(mask & m); i--, m >>= 1);

      return i;
    }

    // Find mask length in bytes.
    static uint8 mask_length(uint8 num_expressions)
    {
      return (uint8)ceil((double)num_expressions / 8.0); 
    }
    
    // Serialize mask into buffer *d_.
    //NOTE: returned mask is in host order
    static uint64 interp_mask(uint8* d, uint8 nbytes)
    {
      uint64 mask = 0;
    
      switch(nbytes)
      {
      case 1:
        mask = (uint64)(((uint8*)d)[0]);
        break;
      case 2:
        mask = (uint64)ntohs(((uint16*)d)[0]);
        break;
      case 4:
        mask = (uint64)ntohl(((uint32*)d)[0]);
        break;
      case 8:
        mask = ntohll(((uint64*)d)[0]);
        break;
      default:
        //FIXME: do something ugly
        break;
      }

      return mask;
    }
    
    // Serialize mask of length nbytes into buffer *d_.
    //NOTE: mask is in host order
    static void serialize_mask(uint8** d_, uint8 nbytes, uint64 mask)
    {
      uint8* d = *d_;
      
      switch(nbytes)
      {
      case 1:
        ((uint8*)d)[0] = (uint8)mask;
        d += 1;
        break;
      case 2:
        ((uint16*)d)[0] = htons((uint16)mask);
        d += 2;
        break;
      case 4:
        ((uint32*)d)[0] = htonl((uint32)mask);
        d += 4;
        break;
      case 8:
        ((uint64*)d)[0] = htonll(mask);
        d += 8;
        break;
      default:
        //FIXME: do something ugly
        break;
      }
      
      *d_ = d;
    }

    void set_requirements_row_(int row, uint64 mask)
    {
      int j;
      uint64 m;

      for(j = 2, m = 1; j < requirements.cols(); j++, m <<= 1)
        requirements(row, j) = ((mask & m) ? 1 : 0);
    }

    uint64 requirements_row_(int row)
    {
      int j;
      uint64 mask = 0;

      for(j = 2; j < requirements.cols(); j++)
        mask |= ((uint64)requirements(row, j) << (j - 2));

      return mask;
    }

    void set_requirements_col_(int col, uint64 mask)
    {
      int i;
      uint64 m;

      for(i = 2, m = 1; i < requirements.rows(); i++, m <<= 1)
        requirements(i, col) = ((mask & m) ? 1 : 0);
    }

    uint64 requirements_col_(int col)
    {
      int i;
      uint64 mask = 0;

      for(i = 2; i < requirements.rows(); i++)
        mask |= ((uint64)requirements(i, col) << (i - 2));

      return mask;
    }
  };

  // Generic jp2/jpx data box (non-superbox)
  class JP2DataBox : public JP2Box
  {
    //friend class JP2File;
    
  protected:
    uint8* DBox;
    bool dbox_allocated;
    uint64 dbox_bytes;
    uint8* dbox_addon;
    bool dbox_addon_allocated;
    uint64 dbox_addon_bytes;
    
  public:
    
    // Get number of bytes in the body of the box (DBox)
    virtual uint64 bytes_dbox(void)
    {
      return dbox_bytes + dbox_addon_bytes;
    }
    
    // Get number of bytes including the header
    virtual uint64 bytes(void)
    {
      return add_header_bytes(bytes_dbox());
    }

    // Print with d intentations
    virtual void print(int d = 0)
    {
      int j;
      
      for(j = 0; j < d; j++)
        std::cout << "  ";
      print_basic();
      switch(TBox)
      {
      case 0x6C626C20: // "lbl\040"
        std::cout << " (label is \"";
        for(j = 0; j < dbox_bytes; j++)
          std::cout << (char)DBox[j];
        std::cout << "\")";
        break;
      case 0x66747970: // "ftyp"
        std::cout << " (filetype is \"";
        for(j = 0; j < 4; j++)
          std::cout << (char)DBox[j];
        std::cout << "\"";
        std::cout << ", compatibility list is \"";
        for(j = 8; j < dbox_bytes; j++)
          std::cout << (char)DBox[j];
        std::cout << "\")";
        break;
      case 0x75756964: // "uuid"
        std::cout << " (id is \"";
        for(j = 0; j < 16; j++)
        {
          std::cout << (j == 0 ? "" : " ") << "0x";
          printf("%02hhx", DBox[j]); //FIXME: better to do this with cout
        }
        std::cout << "\")";
        break;
      default:
        break;
      }
      std::cout << std::endl;
    }
    
    // Serialize into buffer d, which must be large enough
    // (call bytes() and allocate buffer first)
    virtual void serialize(uint8* d)
    {
      serialize_basic(&d);
     
      if(dbox_bytes > 0)
      {
        memcpy(d, DBox, dbox_bytes);
        d += dbox_bytes;
      }
      if(dbox_addon_bytes > 0)
      {
        memcpy(d, dbox_addon, dbox_addon_bytes);
        d += dbox_addon_bytes;
      }
    }

    // Get DBox data
    //NOTE: this does NOT include any dbox_addon that might have been added
    uint8* data(void)
    {
      return DBox;
    }
    
    // Add some data on to the payload
    // 'allocated' is whether to delete the buffer when this box is destroyed
    void addon_dbox(uint8* d, uint64 nbytes, bool allocated)
    {
      if(dbox_addon_bytes > 0 && dbox_addon_allocated)
        delete[] dbox_addon;
      dbox_addon = d;
      dbox_addon_bytes = nbytes;
      dbox_addon_allocated = allocated;
    }
    
    // Construct JP2DataBox from buffer d containing a serialized box (d includes the header)
    // Buffer d is not deleted when this box is destroyed
    JP2DataBox(uint8* d, uint64 nbytes) : JP2Box(interp_type(d)), dbox_addon(0), dbox_addon_allocated(false), dbox_addon_bytes(0)
    {
      uint64 s, r;

      s = interp_bytes(d, nbytes);
      r = header_bytes_to_remove(s);
      DBox = &d[r];
      dbox_allocated = false;
      dbox_bytes = s - r;
    }

    // Construct JP2DataBox from box type, DBox_ payload buffer, and number of payload bytes
    // 'dbox_allocated' is whether to delete the buffer when this box is destroyed
    JP2DataBox(uint32 type, uint8* DBox_, uint64 nbytes, bool dbox_allocated_ = false) : JP2Box(type), dbox_addon(0), dbox_addon_allocated(false), dbox_addon_bytes(0)
    {
      DBox = DBox_;
      dbox_allocated = dbox_allocated_;
      dbox_bytes = nbytes;
    }

    // Construct JP2DataBox from box size {LBox_, XLBox_}, box type TBox_, payload buffer TBox_, and number
    // of bytes remaining in file nbytes
    // 'dbox_allocated' is whether to delete the buffer when this box is destroyed
    // 'byte_order_converted' is whether the byte order of LBox_, TBox_, and XLBox_ have already been converted
    JP2DataBox(uint32 LBox_, uint32 TBox_, uint64 XLBox_, uint8* DBox_, uint64 nbytes, bool dbox_allocated_ = false, bool byte_order_converted = true) : JP2Box(interp_type(TBox_, byte_order_converted)), dbox_addon(0), dbox_addon_allocated(false), dbox_addon_bytes(0)
    {
      uint64 s, r;

      s = interp_bytes(LBox_, XLBox_, nbytes, byte_order_converted);
      r = header_bytes_to_remove(s);
      DBox = DBox_;
      dbox_allocated = dbox_allocated_;
      dbox_bytes = s - r;
    }
    
    // Destructor
    ~JP2DataBox()
    {
      if(DBox && dbox_allocated)
        delete[] DBox;
      if(dbox_addon && dbox_addon_allocated)
        delete[] dbox_addon;
    }
  };
  
  // Generic jp2/jpx superbox
  class JP2SuperBox : public JP2Box
  {
    //friend class JP2File;
    
  public:
    class JP2BoxIterator
    {
    public:
      JP2BoxList::iterator i;
      bool initialized;
      
      void reset(void)
      {
        initialized = false;
      }
      
      JP2BoxIterator() : initialized(false)
      {
      }
      
      ~JP2BoxIterator()
      {
      }
    };
    
  protected:
    JP2BoxList sub_boxes;
    
  public:
    // Get number of bytes in the body of the box (DBox)
    virtual uint64 bytes_dbox(void)
    {
      JP2BoxList::iterator i;
      uint64 retval = 0;
      
      for(i = sub_boxes.begin(); i != sub_boxes.end(); i++)
        retval += (*i)->bytes();
        
      return retval;
    }
    
    // Get number of bytes including the header
    virtual uint64 bytes(void)
    {
      return add_header_bytes(bytes_dbox());
    }

    // Print with d intentations
    virtual void print(int d = 0)
    {
      JP2BoxList::iterator i;
      int j;
      
      for(j = 0; j < d; j++)
        std::cout << "  ";
      print_basic();
      std::cout << std::endl;
      for(i = sub_boxes.begin(); i != sub_boxes.end(); i++)
        (*i)->print(d + 1);
    }
    
    // Serialize into buffer d, which must be large enough
    // (call bytes() and allocate buffer first)
    virtual void serialize(uint8* d)
    {
      JP2Box* b;
      JP2BoxList::iterator i;
      
      serialize_basic(&d);
      
      for(i = sub_boxes.begin(); i != sub_boxes.end(); i++)
      {
        b = *i;
        b->serialize(d);
        d += b->bytes();
      }
    }
    
    // Find the next box with the given type (type = 0 for any box),
    // beginning at position pos (pos = beginning of file by default).
    JP2Box* find_box(uint32 type, JP2BoxIterator* pos = 0)
    {
      JP2BoxList::iterator i;
      JP2Box* b = 0;
      
      if(pos && pos->initialized)
        i = pos->i;
      else
        i = sub_boxes.begin();
        
      for( ; i != sub_boxes.end(); i++)
      {
        if(type == 0 || (*i)->box_type() == type)
        {
          b = *i;
          i++;
          break;
        }
      }
      
      if(pos)
      {
        pos->i = i;
        pos->initialized = true;
      }
        
      return b;
    }
    
    // Insert box b before all other sub-boxes.
    void insert_box_first(JP2Box* b)
    {
      sub_boxes.push_front(b);
    }

    // Insert box b after all other sub-boxes.
    void insert_box_last(JP2Box* b)
    {
      sub_boxes.push_back(b);
    }

    // Insert box b before position pos.
    void insert_box_before(JP2BoxIterator& pos, JP2Box* b)
    {
      sub_boxes.insert(pos.i, b);
    }
    
    // Insert box b after position pos.
    void insert_box_after(JP2BoxIterator& pos, JP2Box* b)
    {
      JP2BoxList::iterator i = pos.i;
      i++;
      sub_boxes.insert(i, b);
    }

    // Construct JP2SuperBox from buffer d containing a serialized box (d includes the header)
    // Buffer d is not deleted when this box is destroyed
    JP2SuperBox(uint8* d, uint64 nbytes) : JP2Box(interp_type(d))
    {
      uint64 s, r;
      JP2Box* b;

      s = interp_bytes(d, nbytes);
      r = header_bytes_to_remove(s);
      find_child_boxes(&d[r], s - r);
    }

    // Construct a JP2SuperBox with only a box type
    JP2SuperBox(uint32 type) : JP2Box(type)
    {
    }
    
    // Destructor
    ~JP2SuperBox()
    {
      JP2BoxList::iterator i;
      
      for(i = sub_boxes.begin(); i != sub_boxes.end(); i++)
        delete (*i);
      sub_boxes.clear();
    }

  protected:
    // Find child boxes of this superbox; d is the beginning of the payload, and there are nbytes
    // bytes remaining in the file
    void find_child_boxes(uint8* d, uint64 nbytes)
    {
      JP2Box* b;
      uint32 b_type;
      uint64 p;
    
      for(p = 0; p < nbytes; p += b->bytes())
      {
        b_type = interp_type(&d[p]);
        switch(b_type)
        {
        case 0x00000000:
           //FIXME: shouldn't do it like this
           fprintf(stderr, "ERROR: Trying to create an illegal JP2 box.\n");
           exit(EXIT_FAILURE);
           break;
        case 0x72726571: // "rreq"
           b = new JP2ReaderRequirementsBox(&d[p], nbytes - p);
           break;
        default:
          if(is_superbox(b_type))
            b = new JP2SuperBox(&d[p], nbytes - p);
          else
            b = new JP2DataBox(&d[p], nbytes - p);
          break;
        }
        sub_boxes.push_back(b);
      }
    }
  };
  
  // jp2/jpx file
  class JP2File : public JP2SuperBox
  {
  public:
    // Get number of bytes including the header
    virtual uint64 bytes(void)
    {
      return bytes_dbox();
    }

    // Print with d intentations
    virtual void print(int d = 0)
    {
      JP2BoxList::iterator i;
      int j;
    
      for(j = 0; j < d; j++)
        std::cout << "  ";
      std::cout << "Boxes:" << std::endl;
      for(i = sub_boxes.begin(); i != sub_boxes.end(); i++)
        (*i)->print(d + 1);
        
      for(j = 0; j < d; j++)
        std::cout << "  ";
      std::cout << "Total " << this->bytes() << " bytes" << std::endl;
    }
    
    // Serialize into buffer d, which must be large enough
    // (call bytes() and allocate buffer first)
    virtual void serialize(uint8* d)
    {
      JP2Box* b;
      JP2BoxList::iterator i;
      uint64 p = 0;
    
      for(i = sub_boxes.begin(); i != sub_boxes.end(); i++)
      {
        b = *i;
        b->serialize(&d[p]);
        p += b->bytes();
      }
    }
    
    // Convert a jp2 file into a jp2-compatible jpx file
    int convert_to_jpx()
    {
      JP2Box* b;
      JP2Box* b2;
      JP2ReaderRequirementsBox* rr;
      uint32* d32;
      char* tmpstr;
      int tmpstr_len;
      JP2BoxIterator pos;
      bool found_jp2, found_jpx, found_jpxb;
      int i;
      
      // find File Type box
      b = find_box(0x66747970); // "ftyp"
      if(!b)
        return -1;
      // BR field of File Type box should be "jpx\040"
      d32 = (uint32*)(((JP2DataBox*)b)->data());
      d32[0] = htonl(0x6A707820); // "jpx\040"
      // both "jpx\040" and "jp2\040" must be in compatibility list; this is also "jpxb"-compatible
      found_jp2 = false;
      found_jpx = false;
      found_jpxb = false;
      for(i = 2; 4 * i < b->bytes_dbox(); i++)
      {
        if(d32[i] == htonl(0x6A703220)) // "jp2\040"
        {
          found_jp2 = true;
          break;
        }
        else if(d32[i] == htonl(0x6A707820)) // "jpx\040"
        {
          found_jpx = true;
          break;
        }
        if(d32[i] == htonl(0x6A706220)) // "jpxb"
        {
          found_jpxb = true;
          break;
        }
      }
      if(!found_jp2 || !found_jpx || !found_jpxb)
      {
        tmpstr = new char[13];
        tmpstr_len = 0;
        tmpstr[0] = '\0';
        if(!found_jp2)
        {
          strcat(tmpstr, "jp2 ");
          tmpstr_len += 4;
        }
        if(!found_jpx)
        {
          strcat(tmpstr, "jpx ");
          tmpstr_len += 4;
        }
        if(!found_jpxb)
        {
          strcat(tmpstr, "jpxb");
          tmpstr_len += 4;
        }
        ((JP2DataBox*)b)->addon_dbox((uint8*)tmpstr, tmpstr_len, true);
        //NOTE: tmpstr will be deleted by JP2DataBox
      }
        
      // Find or add Reader Requirements box
      b = find_box(0x72726571); // "rreq"
      if(!b)
      {
        // Fully Understood Aspects: 1 & 5 & 8 & 12 & 18 & 24 & 31
        // Decode Completely (Display): 1 & 5 & 8 & 12 & 18 & 24 & 31
        rr = new JP2ReaderRequirementsBox;
        rr->add_requirements_row(1); // codestream contains no extensions
        rr->add_requirements_row(5); // codestream is JPEG 2000 as defined by ITU-T Rec. T.800 | ISO/IEC 15444-1
        //FIXME: if have opacity, want to take 8 out and replace it with 9 or 10
        rr->add_requirements_row(8); // no opacity
        rr->add_requirements_row(12); // contiguous codestream
        rr->add_requirements_row(18); // no compositing layers
        rr->add_requirements_row(24); // no animation
        rr->add_requirements_row(31); // no scaling
        rr->add_requirements_col(0x01);
        rr->add_requirements_col(0x02);
        rr->add_requirements_col(0x04);
        rr->add_requirements_col(0x08);
        rr->add_requirements_col(0x10);
        rr->add_requirements_col(0x20);
        rr->add_requirements_col(0x40);
        // find File Type box (Reader Requirements box will go immediately after File Type box)
        pos.reset();
        b = find_box(0x66747970, &pos); // "ftyp"
        if(!b)
          return -1;
        // insert Reader Requirements box immediately after File Type box
        insert_box_after(pos, rr);
        //NOTE: rr will be deleted from sub_boxes list
      }
      
      // find JP2 Header box
      b = find_box(0x6A703268); // "jp2h"
      if(!b)
        return -1;
      // JP2 Header box must contain a Label box (label for the codestream)
      // as its first sub-box
      b2 = ((JP2SuperBox*)b)->find_box(0);
      if(!b2 || b2->box_type() != 0x6C626C20 /*"lbl\040"*/)
      {
        tmpstr = new char[4];
        strcpy(tmpstr, "img");
        b2 = new JP2DataBox(0x6C626C20 /*"lbl\040"*/, (uint8*)tmpstr, 3, true);
        //NOTE: tmpstr will be deleted by JP2DataBox
        ((JP2SuperBox*)b)->insert_box_first(b2);
        //NOTE: b2 will be deleted by JP2SuperBox
      }
      
      return 0;
    }
    
    // Add boxes immediately before the Contiguous Codestream box.
    int add_boxes(const JP2BoxList& boxes)
    {
      JP2DataBox* b;
      JP2BoxIterator pos;
      JP2BoxList::const_iterator i;
      
      // find Contiguous Codestream box
      pos.reset();
      b = (JP2DataBox*)find_box(0x6A703263, &pos); // "jp2c"
      if(!b)
        return -1;
      
      // insert boxes immediately before Contiguous Codestream box
      for(i = boxes.begin(); i != boxes.end(); i++)
        insert_box_before(pos, *i);
      
      return 0;
    }
    
    // Add (fully understood, non-compound) requirements to the
    // Reader Requirements box.
    int add_requirements(const JP2ReaderRequirementsList& req)
    {
      JP2ReaderRequirementsBox* rr;
      JP2ReaderRequirementsList::const_iterator i;
      uint64 mask;
      bool vendor = false;
      bool foundit;
      int j;

      // find Reader Requirements box
      rr = (JP2ReaderRequirementsBox*)find_box(0x72726571); // "rreq"
      if(!rr)
        return -1;
      
      // add requirements to Reader Requirements box
      for(i = req.begin(); i != req.end(); i++)
      {
        // ensure that requirement does not already exist
        foundit = false;
        for(j = 0; j < rr->requirements_rows(); j++)
        {
          if(rr->requirement(j, &vendor) == (*i).first && vendor == (*i).second)
          {
            foundit = true;
            break;
          }
        }
        // if requirement did not already exist, then add it
        if(!foundit)
        {
          mask = 1;
          mask <<= rr->requirements_rows();
          rr->add_requirements_row((*i).first, 0x00, (*i).second);
          rr->add_requirements_col(mask, false);
        }
      }
      
      return 0;
    }
  
    // Construct JP2File from buffer d containing a serialized file
    // Buffer d is not deleted when this object is destroyed
    JP2File(uint8* d, uint64 nbytes) : JP2SuperBox(0x00000000 /*fake box type for root JP2File "box"*/)
    {
      find_child_boxes(d, nbytes);
    }
    
    // Destructor
    ~JP2File()
    {
     // everything that needs to be done is done in the JP2SuperBlock destructor
    }
  };

} // namespace vw

#endif // __VW_FILEIO_JPEG2000_H__
