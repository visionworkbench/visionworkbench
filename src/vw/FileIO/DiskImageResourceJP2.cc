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

/// \file DiskImageResourceJP2.cc
/// 
#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <list>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>   
#include <sstream>   

#include <values.h>			   // for BITSPERBYTE
#include <strings.h>			   // for strncasecmp
#include <arpa/inet.h>			   // ntohl, htonl
#include <endian.h>			   // __BYTE_ORDER

#include <boost/algorithm/string.hpp>

#include <openjpeg.h>
#include <xmlParser.h>

#include <vw/Core/Exception.h>
#include <vw/Cartography/GeoReference.h> //FIXME: should this whole thing go in Cartography?
#include <vw/FileIO/DiskImageResourceJP2.h>

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

using namespace std;
using namespace boost;

extern "C" {
int cio_numbytesleft(opj_cio_t *cio);
}

namespace vw
{

  // These classes parse jp2 and jpx files, upgrade jp2 -> jp2-compatible jpx,
  // and insert GML geospatial data into jpx files.
  // References: jp2 file format: http://www.jpeg.org/public/15444-1annexi.pdf
  //             jpx file format: http://www.jpeg.org/public/15444-2annexm.pdf
  //             GML in jpx: http://www.opengeospatial.org/standards/gmljp2

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

    // Constructor
    JP2Box()
    {
    }

    // Destructor
    ~JP2Box()
    {
    }
  
  protected:
    // Determine whether box with type 'type' is a superbox
    bool is_superbox(uint32 type)
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
    uint32 interp_type(uint8* d)
    {
      return ntohl(*((uint32*)&d[4]));
    }

    // Reverse the bytes in TBox_ if necessary
    uint32 interp_type(uint32 TBox_, bool byte_order_converted)
    {
      return byte_order_converted ? TBox_ : ntohl(TBox_);
    }

    // Determine how many bytes (including the header) are in a box with the given LBox_, XLBox_,
    // and number of bytes remaining in the file. If necessary, reverses the bytes of LBox_ and XLBox_.
    //NOTE: nbytes is assumed to always have had its byte order coverted
    uint64 interp_bytes(uint32 LBox_, uint64 XLBox_, uint64 nbytes, bool byte_order_converted = true)
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
    uint64 interp_bytes(uint8* d, uint64 nbytes)
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
    uint64 header_bytes_to_add(uint64 b)
    {
      return ((b + 8) > 0x00000000ffffffffULL ? 16 : 8);
    }

    // Add the appropriate number of header bytes to a box with b bytes of payload
    uint64 add_header_bytes(uint64 b)
    {
      return b + header_bytes_to_add(b);
    }

    // Determine how many header bytes to remove from a box with b total bytes (b includes the header)
    uint64 header_bytes_to_remove(uint64 b)
    {
      return (b > 0x00000000ffffffffULL ? 16 : 8);
    }

    // Remove the appropriate number of header bytes from a box with b total bytes (b includes the header)
    uint64 remove_header_bytes(uint64 b)
    {
      return b - header_bytes_to_remove(b);
    }

    // Print the box type string corresponding to numeric box type t
    void print_type(uint32 t)
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
    void serialize_basic(uint8** d_)
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

  // Generic jp2/jpx data box (non-superbox)
  class JP2DataBox : public JP2Box
  {
    friend class JP2File;
    
  protected:
    uint8* DBox;
    bool dbox_allocated;
    uint64 dbox_bytes;
    uint8* dbox_addon;
    bool dbox_addon_allocated;
    uint64 dbox_addon_bytes;
    
  public:
    
    // Get number of bytes in the body of the box (DBox)
    uint64 bytes_dbox(void)
    {
      return dbox_bytes + dbox_addon_bytes;
    }
    
    // Get number of bytes including the header
    uint64 bytes(void)
    {
      return add_header_bytes(bytes_dbox());
    }

    // Print with d intentations
    void print(int d = 0)
    {
      int j;
      
      for(j = 0; j < d; j++)
        std::cout << "  ";
      print_basic();
      if(TBox == 0x6C626C20) // "lbl\040"
      {
        std::cout << " (label is \"";
        for(j = 0; j < dbox_bytes; j++)
          std::cout << (char)DBox[j];
        std::cout << "\")";
      }
      std::cout << std::endl;
    }
    
    // Serialize into buffer d, which must be large enough
    // (call bytes() and allocate buffer first)
    void serialize(uint8* d)
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
    
    // Create JP2DataBox from buffer d containing a serialized box (d includes the header)
    // Buffer d is not deleted when this box is destroyed
    JP2DataBox(uint8* d, uint64 nbytes) : dbox_addon(0), dbox_addon_allocated(false), dbox_addon_bytes(0)
    {
      uint64 s, r;

      TBox = interp_type(d);

      s = interp_bytes(d, nbytes);
      r = header_bytes_to_remove(s);
      DBox = &d[r];
      dbox_allocated = false;
      dbox_bytes = s - r;
    }

    // Create JP2DataBox from box type, DBox_ payload buffer, and number of payload bytes
    // 'dbox_allocated' is whether to delete the buffer when this box is destroyed
    JP2DataBox(uint32 type, uint8* DBox_, uint64 nbytes, bool dbox_allocated_ = false) : dbox_addon(0), dbox_addon_allocated(false), dbox_addon_bytes(0)
    {
      TBox = type;

      DBox = DBox_;
      dbox_allocated = dbox_allocated_;
      dbox_bytes = nbytes;
    }

    // Create JP2DataBox from box size {LBox_, XLBox_}, box type TBox_, payload buffer TBox_, and number
    // of bytes remaining in file nbytes
    // 'dbox_allocated' is whether to delete the buffer when this box is destroyed
    // 'byte_order_converted' is whether the byte order of LBox_, TBox_, and XLBox_ have already been converted
    JP2DataBox(uint32 LBox_, uint32 TBox_, uint64 XLBox_, uint8* DBox_, uint64 nbytes, bool dbox_allocated_ = false, bool byte_order_converted = true) : dbox_addon(0), dbox_addon_allocated(false), dbox_addon_bytes(0)
    {
      uint64 s, r;

      TBox = interp_type(TBox_, byte_order_converted);

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
    friend class JP2File;
    
  protected:
    std::list<JP2Box*> sub_boxes;
    
  public:
    // Get number of bytes in the body of the box (DBox)
    uint64 bytes_dbox(void)
    {
      std::list<JP2Box*>::iterator i;
      uint64 retval = 0;
      
      for(i = sub_boxes.begin(); i != sub_boxes.end(); i++)
        retval += (*i)->bytes();
        
      return retval;
    }
    
    // Get number of bytes including the header
    uint64 bytes(void)
    {
      return add_header_bytes(bytes_dbox());
    }

    // Print with d intentations
    void print(int d = 0)
    {
      std::list<JP2Box*>::iterator i;
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
    void serialize(uint8* d)
    {
      JP2Box* b;
      std::list<JP2Box*>::iterator i;
      
      serialize_basic(&d);
      
      for(i = sub_boxes.begin(); i != sub_boxes.end(); i++)
      {
        b = *i;
        b->serialize(d);
        d += b->bytes();
      }
    }

    // Non-initializing constructor
    //NOTE: this should not be used except by the JP2File constructor
    JP2SuperBox()
    {
    }

    // Create JP2SuperBox from buffer d containing a serialized box (d includes the header)
    // Buffer d is not deleted when this box is destroyed
    JP2SuperBox(uint8* d, uint64 nbytes)
    {
      uint64 s, r;
      JP2Box* b;

      TBox = interp_type(d);

      s = interp_bytes(d, nbytes);
      r = header_bytes_to_remove(s);
      find_child_boxes(&d[r], s - r);
    }

    // Construct a JP2SuperBox with only a box type
    JP2SuperBox(uint32 type)
    {
      TBox = type;
    }
    
    // Destructor
    ~JP2SuperBox()
    {
      std::list<JP2Box*>::iterator i;
      
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
      uint64 p;
    
      for(p = 0; p < nbytes; p += b->bytes())
      {
        if(is_superbox(interp_type(&d[p])))
          b = new JP2SuperBox(&d[p], nbytes - p);
        else
          b = new JP2DataBox(&d[p], nbytes - p);
        sub_boxes.push_back(b);
      }
    }
  };
  
  // Reader Requirements box
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
      for(i = 0; i < NSF_; i++)
      {
        SF[i] = (uint16*)d; *(SF[i]) = 0; d += sizeof(uint16);
        SM[i] = (MaskT*)d; *(SM[i]) = 0; d += sizeof(MaskT);
      }
      NVF = (uint16*)d; *NVF = NVF_; d += sizeof(uint16);
      for(i = 0; i < NVF_; i++)
      {
        VF[i] = (uint16*)d; *(VF[i]) = 0; d += sizeof(uint16);
        VM[i] = (MaskT*)d; *(VM[i]) = 0; d += sizeof(MaskT);
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
  
  // jp2/jpx file
  class JP2File : public JP2SuperBox
  {
  public:
    // Get number of bytes including the header
    uint64 bytes(void)
    {
      return bytes_dbox();
    }

    // Print with d intentations
    void print(int d = 0)
    {
      std::list<JP2Box*>::iterator i;
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
    void serialize(uint8* d)
    {
      JP2Box* b;
      std::list<JP2Box*>::iterator i;
    
      for(i = sub_boxes.begin(); i != sub_boxes.end(); i++)
      {
        b = *i;
        b->serialize(d);
        d += b->bytes();
      }
    }
    
    // Convert a jp2 file into a jp2-compatible jpx file
    // 'contains_gml' is whether the jpx file will contain GML
    int convert_to_jpx(bool contains_gml)
    {
      JP2Box* b;
      JP2Box* b2;
      uint32* d32;
      char* tmpstr;
      std::list<JP2Box*>::iterator i;
      int j;
      bool foundit;
      
      // find File Type box
      foundit = false;
      for(b = 0, i = sub_boxes.begin(); i != sub_boxes.end(); i++)
      {
        b = *i;
        if(b->box_type() == 0x66747970) // "ftyp"
        {
          foundit = true;
          break;
        }
      }
      if(!foundit)
        return -1;
      // BR field of File Type box should be "jpx\040"
      d32 = (uint32*)(((JP2DataBox*)b)->data());
      d32[0] = htonl(0x6A707820); // "jpx\040"
      // both "jpx\040" and "jp2\040" must be in compatibility list; this is also "jpxb"-compatible
      foundit = false;
      for(j = 2; 4 * j < b->bytes_dbox(); j++)
      {
        if(d32[j] == htonl(0x6A703220)) // "jp2\040"
        {
          foundit = true;
          break;
        }
      }
      if(foundit)
      {
        tmpstr = new char[9];
        strcpy(tmpstr, "jpx jpxb");
        ((JP2DataBox*)b)->addon_dbox((uint8*)tmpstr, 8, true);
        //NOTE: tmpstr will be deleted by JP2DataBox
      }
      else
      {
        tmpstr = new char[13];
        strcpy(tmpstr, "jpx jpxbjp2 ");
        ((JP2DataBox*)b)->addon_dbox((uint8*)tmpstr, 12, true);
        //NOTE: tmpstr will be deleted by JP2DataBox
      }
        
      // add Reader Requirements box
      // Fully Understood Aspects: 1 & 5 & 8 & 12 & 18 & 24 & 31 (& 67)
      // Decode Completely (Display): 1 & 5 & 8 & 12 & 18 & 24 & 31
      JP2ReaderRequirements<uint8> rr(contains_gml ? 8 : 7, 0, false);
      //NOTE: if we have to upgrade these masks to uint16, we need some htons()
      *(rr.FUAM) = (uint8)0xFF;
      *(rr.DCM)  = (uint8)0xFE;
      *(rr.SF[0]) = htons(1); // codestream contains no extensions
      *(rr.SM[0]) = (uint8)0x80;
      *(rr.SF[1]) = htons(5); // codestream is JPEG 2000 as defined by ITU-T Rec. T.800 | ISO/IEC 15444-1
      *(rr.SM[1]) = (uint8)0x40;
      //FIXME: if have opacity, want to take 8 out and replace it with 9 or 10
      *(rr.SF[2]) = htons(8); // no opacity
      *(rr.SM[2]) = (uint8)0x20;
      *(rr.SF[3]) = htons(12); // contiguous codestream
      *(rr.SM[3]) = (uint8)0x10;
      *(rr.SF[4]) = htons(18); // no compositing layers
      *(rr.SM[4]) = (uint8)0x08;
      *(rr.SF[5]) = htons(24); // no animation
      *(rr.SM[5]) = (uint8)0x04;
      *(rr.SF[6]) = htons(31); // no scaling
      *(rr.SM[6]) = (uint8)0x02;
      if(contains_gml)
      {
        *(rr.SF[7]) = htons(67); // GML
        *(rr.SM[7]) = (uint8)0x01;
      }
      // find File Type box (Reader Requirements box will go immediately after File Type box)
      foundit = false;
      for(b = 0, i = sub_boxes.begin(); i != sub_boxes.end(); i++)
      {
        b = *i;
        if(b->box_type() == 0x66747970) // "ftyp"
        {
          foundit = true;
          break;
        }
      }
      if(!foundit)
        return -1;
      i++;
      // create Reader Requirements box
      b = new JP2DataBox(0x72726571 /*"rreq"*/, rr.dat, rr.bytes, true);
      //NOTE: rr.dat will be deleted by JP2DataBox
      // insert Reader Requirements box immediately after File Type box
      sub_boxes.insert(i, b);
      //NOTE: b will be deleted from sub_boxes list
      
      // find JP2 Header box
      foundit = false;
      for(b = 0, i = sub_boxes.begin(); i != sub_boxes.end(); i++)
      {
        b = *i;
        if(b->box_type() == 0x6A703268) // "jp2h"
        {
          foundit = true;
          break;
        }
      }
      if(!foundit)
        return -1;
      // JP2 Header box must contain a Label box (label for the codestream)
      tmpstr = new char[4];
      strcpy(tmpstr, "img");
      b2 = new JP2DataBox(0x6C626C20 /*"lbl\040"*/, (uint8*)tmpstr, 3, true);
      //NOTE: tmpstr will be deleted by JP2DataBox
      ((JP2SuperBox*)b)->sub_boxes.push_front(b2);
      //NOTE: b2 will be deleted by JP2SuperBox
      
      return 0;
    }

    int add_gml(uint8* d, uint64 nbytes, bool allocated)
    {
      JP2DataBox* b;
      JP2SuperBox* a;
      JP2SuperBox* a2;
      //JP2Box* b2;
      //uint32* d32;
      char* tmpstr;
      std::list<JP2Box*>::iterator i;
      //int j;
      bool foundit;
      
      // find Contiguous Codestream box
      foundit = false;
      for(i = sub_boxes.begin(); i != sub_boxes.end(); i++)
      {
        if((*i)->box_type() == 0x6A703263) // "jp2c"
        {
          foundit = true;
          break;
        }
      }
      if(!foundit)
        return -1;

      // add outer Association box before Contiguous Codestream box
      a = new JP2SuperBox(0x61736F63 /*"asoc"*/);
      sub_boxes.insert(i, a);
      //NOTE: a will be deleted from sub_boxes list

      // add Label box with "gml.data" to outer Association box
      tmpstr = new char[9];
      strcpy(tmpstr, "gml.data");
      b = new JP2DataBox(0x6C626C20 /*"lbl\040"*/, (uint8*)tmpstr, 8, true);
      //NOTE: tmpstr will be deleted by JP2DataBox
      a->sub_boxes.push_back(b);
      //NOTE: b will be deleted by JP2SuperBox

      // add inner Association box to outer Association box
      a2 = new JP2SuperBox(0x61736F63 /*"asoc"*/);
      a->sub_boxes.push_back(a2);
      //NOTE: a2 will be deleted by JP2SuperBox

      // add Label box with "gml.root-instance" to inner Association box
      tmpstr = new char[18];
      strcpy(tmpstr, "gml.root-instance");
      b = new JP2DataBox(0x6C626C20 /*"lbl\040"*/, (uint8*)tmpstr, 17, true);
      //NOTE: tmpstr will be deleted by JP2DataBox
      a2->sub_boxes.push_back(b);
      //NOTE: b will be deleted by JP2SuperBox

      // add XML box with GML root-instance data to inner Association box
      b = new JP2DataBox(0x786D6C20 /*"xml\040"*/, d, nbytes, allocated);
      //NOTE: d will be deleted by JP2DataBox iff allocated
      a2->sub_boxes.push_back(b);
      //NOTE: b will be deleted by JP2SuperBox

      return 0;
    }
  
    // Create JP2DataBox from buffer d containing a serialized file
    // Buffer d is not deleted when this object is destroyed
    JP2File(uint8* d, uint64 nbytes) : JP2SuperBox()
    {
      TBox = 0x00000000; // fake box type for root JP2File "box"
      find_child_boxes(d, nbytes);
    }
    
    // Destructor
    ~JP2File()
    {
     // everything that needs to be done is done in the JP2SuperBlock destructor
    }
  };
  
  
  // Make GML from a GeoReference
  //NOTE: caller must free returned string
  //FIXME: needs to take a GeoReference instead of just making boilerplate GML
  char* make_gml(/*cartography::GeoReference const& georef*/void)
  {
    char* s;
    XMLNode t, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10;
    
    t = XMLNode::createXMLTopNode("xml", TRUE);
    t.addAttribute("version", "1.0");
    t.addAttribute("encoding", "UTF-8");
    
    //NOTE: this was initially generated from a sample gml file using gml_print_code()
    
    n1 = t.addChild("gml:FeatureCollection");
    n1.addAttribute("xmlns", "http://www.opengis.net/gml");
    n1.addAttribute("xmlns:gml", "http://www.opengis.net/gml");
    n1.addAttribute("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
    n1.addAttribute("xsi:schemaLocation", "http://www.opengis.net/gml gmlJP2Profile.xsd"); //FIXME: see p. 41 of GMLJP2 standard, which appears to say that gmlJP2Profile.xsd must be included in the jp2 file and cannot be linked like in the commented-out line below
    //n1.addAttribute("xsi:schemaLocation", "http://www.opengis.net/gml http://schemas.opengis.net/gml/3.1.1/profiles/gmlJP2Profile/1.0.0/gmlJP2Profile.xsd");
    
    n2 = n1.addChild("gml:boundedBy");
    n3 = n2.addChild("gml:Envelope");
    n4 = n3.addChild("gml:lowerCorner");
    n4.addText("270379.500 3942462.000");
    n4 = n3.addChild("gml:upperCorner");
    n4.addText("518842.500 3942462.000");
    
    n2 = n1.addChild("gml:featureMember");
    n3 = n2.addChild("gml:FeatureCollection");
    
    n4 = n3.addChild("gml:boundedBy");
    n5 = n4.addChild("gml:Envelope");
    n6 = n5.addChild("gml:lowerCorner");
    n6.addText("270379.500 3942462.000");
    n6 = n5.addChild("gml:upperCorner");
    n6.addText("518842.500 3942462.000");
    
    n4 = n3.addChild("gml:featureMember");
    n5 = n4.addChild("gml:RectifiedGridCoverage");
    n5.addAttribute("dimension", "2");
    n5.addAttribute("gml:id", "RGC0001");
    n6 = n5.addChild("gml:description");
    n6.addText("This GMLJP2 Minimal Root Instance contains a GML Rectified Grid. The rectified grid is embedded in a RectifiedGridCoverage with generic range parameters (to be ignored).");
    n6 = n5.addChild("gml:rectifiedGridDomain");
    n7 = n6.addChild("gml:RectifiedGrid");
    n7.addAttribute("dimension", "2");
    
    n8 = n7.addChild("gml:limits");
    n9 = n8.addChild("gml:GridEnvelope");
    n10 = n9.addChild("gml:low");
    n10.addText("0 0");
    n10 = n9.addChild("gml:high");
    n10.addText("8718 7812");
    
    n8 = n7.addChild("gml:axisName");
    n8.addText("x");
    n8 = n7.addChild("gml:axisName");
    n8.addText("y");
    
    n8 = n7.addChild("gml:origin");
    n9 = n8.addChild("gml:Point");
    n9.addAttribute("gml:id", "Pt001");
    n9.addAttribute("srsName", "urn:ogc:def:crs:EPSG:6.6:32612");
    n10 = n9.addChild("gml:description");
    n10.addText("\"Upper-left\" image origin");
    n10 = n9.addChild("gml:coordinates");
    n10.addText("270379.500000, 3942462.000000");
    
    n8 = n7.addChild("gml:offsetVector");
    n8.addAttribute("srsName", "urn:ogc:def:crs:EPSG:6.6:32612");
    n8.addText("28.5 0");
    n8 = n7.addChild("gml:offsetVector");
    n8.addAttribute("srsName", "urn:ogc:def:crs:EPSG:6.6:32612");
    n8.addText("0 28.5");
    
    n6 = n5.addChild("gml:rangeSet");
    n7 = n6.addChild("gml:File");
    n8 = n7.addChild("gml:rangeParameters");
    n9 = n8.addChild("gml:QuantityList");
    n9.addAttribute("uom", "urn:ogc:def:crs:EPSG:6.6:32612");
    n9.addText("inapplicable");
    n8 = n7.addChild("gml:fileName");
    n8.addText("Not Applicable");
    n8 = n7.addChild("gml:fileStructure");
    n8.addText("Record Interleaved");
    
    n6 = n5.addChild("gml:coverageFunction");
    n7 = n6.addChild("gml:GridFunction");
    n8 = n7.addChild("gml:sequenceRule");
    n8.addAttribute("order", "+x+y");
    n8.addText("Linear");
    n8 = n7.addChild("gml:startPoint");
    n8.addText("0 0");
    
    s = t.createXMLString(false);
    return s;
  }
  
  // Recursive implementation of gml_print_code()
  void gml_print_code_recursive(XMLNode const& n, int d)
  {
    int num_attributes, num_text, num_children;
    int i;
    
    std::cout << "n" << d << " = n" << (d - 1) << ".addChild(\"" << (strchr(n.getName(), ':') ? "" : "gml:") << n.getName() << "\");" << std::endl;
    
    num_attributes = n.nAttribute();
    for(i = 0; i < num_attributes; i++)
      std::cout << "n" << d << ".addAttribute(\"" << n.getAttributeName(i) << "\", \"" << n.getAttributeValue(i) << "\");" << std::endl;
      
    num_text = n.nText();
    for(i = 0; i < num_text; i++)
      std::cout << "n" << d << ".addText(\"" << n.getText(i) << "\");" << std::endl;
  
    num_children = n.nChildNode();
    for(i = 0; i < num_children; i++)
      gml_print_code_recursive(n.getChildNode(i), d + 1);
  }
  
  // Print C++ code to create the XML file fn
  void gml_print_code(const char* fn)
  {
    XMLNode t;
     
    t = XMLNode::openFileHelper(fn, "FeatureCollection");
    gml_print_code_recursive(t, 1);
  }
  
  

  /// Close the file when the object is destroyed
  DiskImageResourceJP2::~DiskImageResourceJP2()
  {
    this->flush();
  }
 
  /// Flush the buffered data to disk
  void
  DiskImageResourceJP2::flush()
  {
  }

#if 0
  static int
  get_file_format(char *filename)
  {
    static const char *extension[] = { "j2k", "jp2", "jpt" };
    static const int format[] = { J2K_CFMT, JP2_CFMT, JPT_CFMT };
    char * ext = strrchr(filename, '.') + 1;

    if (ext != 0)
    {
      for (int i = 0; i < sizeof(format); i++)
      {
	if (strncasecmp(ext, extension[i], 3) == 0)
	  return format[i];
      }
    }

    return -1;
  }
#endif
  // sample error callback expecting a FILE* client object
  static void
  error_callback(const char *msg, void *client_data)
  {
    FILE *stream = (FILE*)client_data;
    fprintf(stream, "[ERROR] %s", msg);
  }
  // sample warning callback expecting a FILE* client object
  void
  warning_callback(const char *msg, void *client_data)
  {
    FILE *stream = (FILE*)client_data;
    fprintf(stream, "[WARNING] %s", msg);
  }
  // sample debug callback expecting a FILE* client object
  void
  info_callback(const char *msg, void *client_data)
  {
    FILE *stream = (FILE*)client_data;
    fprintf(stream, "[INFO] %s", msg);
  }

  static void
  setup_event_handling(opj_common_ptr cinfo)
  {
    opj_event_mgr_t event_manager;

    // configure the event callbacks (not required)
    // setting of each callback is optionnal 
    memset(&event_manager, 0, sizeof(opj_event_mgr_t));
    event_manager.error_handler = error_callback;
    event_manager.warning_handler = warning_callback;
    event_manager.info_handler = info_callback;
    // catch events using our callbacks and give a local context
    opj_set_event_mgr((opj_common_ptr)cinfo, &event_manager, stderr);
  }

  static opj_image_t *
  init_jpeg2k_image(opj_cparameters_t &compress_parameters,
		    int width, int height, int num_components,
		    OPJ_COLOR_SPACE color_space,
		    ImageBuffer &dst)
  {
    cout << "-- init_jpeg2k_image(): start" << endl;
    opj_image_t *new_image = new opj_image_t;

    if (new_image != 0)
    {
      new_image->color_space = color_space;
      new_image->numcomps = num_components;

      // allocate memory for the per-component information
      new_image->comps = new opj_image_comp_t [new_image->numcomps];
	
      if (new_image->comps == 0)
      {
	delete new_image;
	return 0;
      }

      int dx = compress_parameters.subsampling_dx;
      int dy = compress_parameters.subsampling_dy;

      // init the individual image components
      for (int i = 0; i < num_components; i++)
      {
	opj_image_comp_t *component = &new_image->comps[i];

	component->dx = dx;
	component->dy = dy;
	component->w = width;
	component->h = height;
	component->x0 = 0;
	component->y0 = 0;
	component->prec = 8;
	component->bpp = 8;
	component->sgnd = 0;
	component->resno_decoded = 0;
	component->factor = 0;

	component->data = (int*) calloc(width * height, sizeof(int));
	if (component->data == 0)
	{
	  cout << "init_jpeg2k_image(): could not allocate memory for data"
	       << endl;
	  delete [] new_image->comps;
	  delete new_image;
	  return 0;
	}
	// maybe we'll eventaully be able to just point it at the data...
// 	component->data = (int *) dst(0, 0, i);
	cout << "init_jpeg2k_image(): initializing data" << endl;
	// Only have had luck with 8 bit images so far...
	uint8 *values = (uint8 *) dst(0, 0, i);
	int minval = values[0];
	int maxval = minval;
	cout << "values[0] = " << int(values[0]) << endl;
	for (int index = 0; index < width * height; index++)
	{
	  component->data[index] = (int) values[index];
	  if (minval > component->data[index])
	    minval = component->data[index];
	  if (maxval < component->data[index])
	    maxval = component->data[index];;
	}
	cout << "done" << endl;
	cout << "minval = " << minval << endl;
	cout << "maxval = " << maxval << endl;
      }

      new_image->x0 = compress_parameters.image_offset_x0;
      new_image->y0 = compress_parameters.image_offset_y0;
      new_image->x1 = new_image->x0 + ((width - 1) * dx) + 1;
      new_image->y1 = new_image->y0 + ((height - 1) * dy) + 1;
    }

    cout << "-- init_jpeg2k_image(): end" << endl;

    return new_image;
  }

  /// Bind the resource to a file for reading.  Confirm that we can open
  /// the file and that it has a sane pixel format.  
  void
  DiskImageResourceJP2::open(std::string const& filename)
  {
    vw_throw(NoImplErr() << "DiskImageResourceJP2 (read) Error: "
	     "JPEG2000 file formats are not supported yet!");
#if 0
    ifstream image_file(filename.c_str(), ios::in | ios::binary);

    if (!image_file)
      throw IOErr() << "DiskImageResourceJP2::open(): could not open \""
		    << filename << "\".";

    // - - - - - - - - - - - -
    // decode the code-stream
    //
    // Note: A JPEG 2000 file (extension .jp2 and magic number
    // 0000000C6A5020200D0A870A) consists of a number of "boxes" that
    // contain metadata or code streams (i.e., encoded image data).
    //
    // A file with a .j2k extension contains only a "raw" code stream,
    // which actually has a header:
    //
    // FF4F -- Start of Codestream (SOC) marker
    // FF51 -- Size (SIZ) marker
    // #### -- Length of SIZ segment (Lsiz)
    // ...  -- SIZ parameters
    // FF5C -- Quantization Default (QCD) marker
    // #### -- Length of QCD segment (Lqcd)
    // ...  -- QCD parameters
    // FF52 -- Coding style Default (COD) marker
    // #### -- Length of COD segment (Lcod)
    // ...  -- COD parameters
    // Optional marker segments COC (FF53), QCC (FF5D), RGN (FF5E),
    // POD (FF5F), PPM (FF60), PLM (FF57), TLM (FF55), CME (Comment
    // and Extension info, FF64)
    //
    // followed by encoded image data:
    //
    // FF90 -- Start of Tile (SOT) marker
    // ...
    // FF93 -- Start of data (SOD) marker
    // ...
    // FFD9 -- End of Codestream (EOC) marker
    //
    // - - - - - - - - - - - -

    opj_dinfo_t* dinfo = 0;		   // handle to a decompressor
    opj_cio_t *cio = 0;
    opj_image_t *image = 0;

    int file_format = get_file_format(file_format.c_str();

    // We can handle *.j2k, *.jp2, *.jpc or *.jpt files
    switch (file_format)
    {
    case J2K_CFMT:			   // JPEG-2000 codestream
      {
	// get a decoder handle
	dinfo = opj_create_decompress(CODEC_J2K);
	
	// catch events using our callbacks and give a local context
	opj_set_event_mgr((opj_common_ptr)dinfo, &event_mgr, stderr);

	// setup the decoder decoding parameters using user parameters
	opj_setup_decoder(dinfo, &parameters);

	// open a byte stream
	cio = opj_cio_open((opj_common_ptr)dinfo, src, file_length);

	// decode the stream and fill the image structure
	image = opj_decode(dinfo, cio);
	if (!image)
	{
	  cerr << "ERROR in DiskImageResourceJP2::open(): "
	       << "failed to decode image!" << endl;
	  opj_destroy_decompress(dinfo);
	  opj_cio_close(cio);
	  return 1;
	}
	
	// close the byte stream
	opj_cio_close(cio);
      }
      break;
    case JP2_CFMT:			   // JPEG 2000 compressed image data
      {
	// get a decoder handle
	dinfo = opj_create_decompress(CODEC_JP2);
	
	// catch events using our callbacks and give a local context
	opj_set_event_mgr((opj_common_ptr)dinfo, &event_mgr, stderr);

	// setup the decoder decoding parameters using the current
	// image and using user parameters
	opj_setup_decoder(dinfo, &parameters);

	// open a byte stream
	cio = opj_cio_open((opj_common_ptr)dinfo, src, file_length);

	// decode the stream and fill the image structure
	image = opj_decode(dinfo, cio);
	if(!image)
	{
	  cerr << "ERROR in DiskImageResourceJP2::open(): "
	       << "failed to decode image!" << endl;
	  opj_destroy_decompress(dinfo);
	  opj_cio_close(cio);
	  return 1;
	}

	// close the byte stream
	opj_cio_close(cio);
      }
      break;
    case JPT_CFMT:			   // JPEG 2000, JPIP
      {
	// get a decoder handle
	dinfo = opj_create_decompress(CODEC_JPT);
	
	// catch events using our callbacks and give a local context
	opj_set_event_mgr((opj_common_ptr)dinfo, &event_mgr, stderr);

	// setup the decoder decoding parameters using user parameters
	opj_setup_decoder(dinfo, &parameters);

	// open a byte stream
	cio = opj_cio_open((opj_common_ptr)dinfo, src, file_length);

	// decode the stream and fill the image structure
	image = opj_decode(dinfo, cio);
	if (!image)
	{
	  cerr << "ERROR in DiskImageResourceJP2::open(): "
	       << "failed to decode image!" << endl;
	  opj_destroy_decompress(dinfo);
	  opj_cio_close(cio);
	  return 1;
	}

	// close the byte stream
	opj_cio_close(cio);
      }
      break;
    default:
      cerr << "ERROR in DiskImageResourceJP2::open(): "
	   << "unknown JPEG2000 image format variant" << endl;
      return 1;
      break;
    }

    m_format.planes = 1;
    m_format.pixel_format = VW_PIXEL_GRAY;
    m_format.cols = header.bytesPerScanline / bytes_per_pixel;
    m_format.rows = header.numScanlines;
 
    image_file.close();
#endif
  }

  /// Bind the resource to a file for writing.
  void
  DiskImageResourceJP2::create(std::string const& filename,
			       ImageFormat const& format)
  {
    cout << "DiskImageResourceJP2::create(): enter" << endl;

    VW_ASSERT(format.planes == 1 || format.pixel_format == VW_PIXEL_SCALAR,
	      NoImplErr() << "DiskImageResourceOpenJP2: Cannot create "
	      << filename << "\n\tThe image cannot have both multiple "
	      << "channels and multiple planes.\n");

    m_filename = filename;
    m_format = format;

    std::cout << "DiskImageResourceJP2::create(): channel type = " << std::flush;

    switch (format.channel_type)
    {
    case VW_CHANNEL_UNKNOWN:
      cout << "VW_CHANNEL_UNKNOWN" << endl;
      break;
    case VW_CHANNEL_INT8:
      cout << "VW_CHANNEL_INT8" << endl;
      break;
    case VW_CHANNEL_UINT8:
      cout << "VW_CHANNEL_UINT8" << endl;
      break;
    case VW_CHANNEL_INT16:
      cout << "VW_CHANNEL_INT16" << endl;
      break;
    case VW_CHANNEL_UINT16:
      cout << "VW_CHANNEL_UINT16" << endl;
      break;
    case VW_CHANNEL_INT32:
      cout << "VW_CHANNEL_INT32" << endl;
      break;
    case VW_CHANNEL_UINT32:
      cout << "VW_CHANNEL_UINT32" << endl;
      break;
    case VW_CHANNEL_INT64:
      cout << "VW_CHANNEL_INT64" << endl;
      break;
    case VW_CHANNEL_UINT64:
      cout << "VW_CHANNEL_UINT64" << endl;
      break;
    case VW_CHANNEL_FLOAT16:
      cout << "VW_CHANNEL_FLOAT16" << endl;
      break;
    case VW_CHANNEL_FLOAT32:
      cout << "VW_CHANNEL_FLOAT32" << endl;
      break;
    case VW_CHANNEL_FLOAT64:
      cout << "VW_CHANNEL_FLOAT64" << endl;
      break;
    case VW_CHANNEL_BOOL:
      cout << "VW_CHANNEL_BOOL" << endl;
      break;
    case VW_CHANNEL_USER:
      cout << "VW_CHANNEL_USER" << endl;
      break;
    }

    // The open JPEG code only seems to work with 8 bit images right now...
    m_format.channel_type = VW_CHANNEL_UINT8;
    m_format.planes = std::max(format.planes,
			       num_channels(format.pixel_format));

    cout << "DiskImageResourceJP2::create(): return" << endl;
  }

  /// Read the disk image into the given buffer.
  void
    DiskImageResourceJP2::read(ImageBuffer const& dest,
			       BBox2i const& bbox) const
  {
    VW_ASSERT((bbox.width() == int(cols())) && (bbox.height()==int(rows())),
	      NoImplErr() << "DiskImageResourceOpenJP2 does not support"
	      " partial reads." );
#if 0
    VW_ASSERT((dest.format.cols == cols()) && (dest.format.rows == rows()),
	      IOErr() << "Buffer has wrong dimensions in JP2 read.");
  
    ifstream image_file(m_filename.c_str(), ios::in | ios::binary);

    if (!image_file)
      vw_throw(IOErr() << "DiskImageResourceJP2::read(): "
	       "failed to open \"" << m_filename << "\" for reading!");

    // Set the file offset to the position of the first image byte
    // (read in from the previously opened JP2 header)
    image_file.seekg(m_image_data_offset, ios::beg);

    image_file.read((char *) image_data, bytes_per_pixel * total_pixels);

    if (image_file.bad())
      vw_throw(IOErr() << "DiskImageResourceJP2::read(): "
	       "an unrecoverable error occured while reading the image data.");

    // Create image buffer from the JP2 data.
    ImageBuffer src;
    src.data = image_data;
    src.format = m_format;
    src.cstride = bytes_per_pixel;
    src.rstride = bytes_per_pixel * m_format.cols;
    src.pstride = bytes_per_pixel * m_format.cols * m_format.rows;
  
    convert(dest, src);

    delete[] image_data;

    image_file.close();
#endif
  }

  // Write the given buffer into the disk image.
  void
    DiskImageResourceJP2::write(ImageBuffer const& src, BBox2i const& bbox)
  {
    VW_ASSERT((bbox.width()==int(cols())) && (bbox.height()==int(rows())),
	      NoImplErr() << "DiskImageResourceOpenJP2 does not support"
	      " partial writes." );

    VW_ASSERT((src.format.cols == cols()) && (src.format.rows == rows()),
	      IOErr() << "Buffer has wrong dimensions in "
	      "DiskImageResourceJP2::write().");

    VW_ASSERT(src.format.planes == 1, IOErr() << "Buffer has multiple planes "
	      "in DiskImageResourceJP2::write(). Only gray scale "
	      "images currently handled.");
  
    // This is pretty simple since we always write 32 bit integer
    // files.  Note that we handle multi-channel images with
    // interleaved planes.
    // We've already ensured that either planes == 1 or channels == 1.
    // Can't seem to get anything other than 8 bit JP2s working right now...
    ImageView<uint8> jp2_image_view(m_format.cols, m_format.rows,
				    m_format.planes);
    ImageBuffer dst = jp2_image_view.buffer();

    convert(dst, src);

    cout << "DiskImageResourceJP2::write(): "
	 << "initializing encoder and parameters" << endl;


    // For the moment we only handle grayscale images
    OPJ_COLOR_SPACE color_space = CLRSPC_GRAY;
    int num_components = 1;
    opj_cparameters_t coding_parameters;

    opj_set_default_encoder_parameters(&coding_parameters);
    coding_parameters.cp_comment = "Created by OpenJPEG version 0.9";
    // if no rate entered, lossless by default
    if (coding_parameters.tcp_numlayers == 0)
    {
      cout << "DiskImageResourceJP2::write(): encoding lossless!" << endl;
      coding_parameters.tcp_rates[0] = 0;	// MOD antonin : losslessbug
      coding_parameters.tcp_numlayers++;
      coding_parameters.cp_disto_alloc = 1;
    }
    // experiment with tiles -- LJE
//     coding_parameters.cp_tdx = 128;
//     coding_parameters.cp_tdy = 128;
//     coding_parameters.cp_tx0 = 0;
//     coding_parameters.cp_ty0 = 0;
//     coding_parameters.tile_size_on = true;

    opj_image_t *jpeg2k_image = init_jpeg2k_image(coding_parameters,
						  cols(), rows(),
						  num_components,
						  color_space, dst);

    int *values = jpeg2k_image->comps[0].data;
    int minval = values[0];
    int maxval = minval;
    cout << "DiskImageResourceJP2::write(): values[0] = " << values[0] << endl;
    cout << "DiskImageResourceJP2::write(): width = " << jpeg2k_image->comps[0].w << endl;
    cout << "DiskImageResourceJP2::write(): height = " << jpeg2k_image->comps[0].h
	 << endl;
    for (int index = 0;
	 index < jpeg2k_image->comps[0].w * jpeg2k_image->comps[0].h;
	 index++)
    {
      if (minval > values[index])
	minval = values[index];
      if (maxval < values[index])
	maxval = values[index];;
    }
    cout << "DiskImageResourceJP2::write(): minval = " << minval << endl;
    cout << "DiskImageResourceJP2::write(): maxval = " << maxval << endl;

    cout << "DiskImageResourceJP2::write(): "
	 << "creating a JP2 encoder" << endl;

    // get a JP2 compression context pointer
    opj_cinfo_t* cinfo = opj_create_compress(CODEC_JP2);

    // Set up event manager (not necessary)
    setup_event_handling(opj_common_ptr(cinfo));

    cout << "DiskImageResourceJP2::write(): "
	 << "setting up JP2 encoder" << endl;


    cout << "num resolutions = " << coding_parameters.numresolution << endl;
    // setup the encoder parameters using the current image
    opj_setup_encoder(cinfo, &coding_parameters, jpeg2k_image);

    cout << "DiskImageResourceJP2::write(): "
	 << "opening codec i/o stream" << endl;

    // open a codec input/output byte stream for writing
    opj_cio_t *cio = opj_cio_open(opj_common_ptr(cinfo), 0, 0);
    if (cio == 0)
    {
      cout << "DiskImageResourceJP2::write(): "
	   << "failed to open a codec input output stream." << endl;

      vw_throw(IOErr() << "DiskImageResourceJP2::write(): "
	       "failed to open a codec input output stream.");
    }

    cout << "DiskImageResourceJP2::write(): "
	 << "encoding image" << endl;

    cout << "jpeg2k_image x0 = " << jpeg2k_image->x0 << endl;
    cout << "jpeg2k_image y0 = " << jpeg2k_image->y0 << endl;
    cout << "jpeg2k_image width = " << jpeg2k_image->x1 << endl;
    cout << "jpeg2k_image height = " << jpeg2k_image->y1 << endl;
    cout << "jpeg2k_image->comps = " << jpeg2k_image->comps << endl;
    cout << "jpeg2k_image->numcomps = " << jpeg2k_image->numcomps << endl;

    bool success = opj_encode(cinfo, cio, jpeg2k_image,
			      coding_parameters.index);
    if (!success)
    {
      cout << "DiskImageResourceJP2::write(): "
	   << "failed to encode image." << endl;

      opj_cio_close(cio);
      vw_throw(IOErr() << "DiskImageResourceJP2::write(): "
	       "failed to encode image.");
    }

    int codestream_length = cio_tell(cio);
    int bytes_to_end =  cio_numbytesleft(cio);
    
    bool convert_to_jpx = true;
    uint8* output_buffer = 0;
    if(convert_to_jpx)
    { // Convert to jpx, and add GML if desired
      JP2File jp2f(cio->buffer, codestream_length);
      //jp2f.print();
      if(jp2f.convert_to_jpx(true /*will contain GML*/) != 0)
        vw_throw(IOErr() << "DiskImageResourceJP2::write(): "
	       "failed to convert jp2 to backward-compatible jpx");
      //gml_print_code("/home/ttemplet/gmljp2/gmlsamples/minimalgml.xml");
      //cartography::GeoReference georef; //FIXME: libvwFileIO does not include GeoReference.cc, so this causes a linker error
      char* gml = make_gml(/*georef*/);
      if(jp2f.add_gml((uint8*)gml, strlen(gml), false) != 0)
        vw_throw(IOErr() << "DiskImageResourceJP2::write(): "
	       "failed to add GML to jpx");
      //jp2f.print();
      codestream_length = jp2f.bytes();
      output_buffer = new uint8[codestream_length];
      jp2f.serialize(output_buffer);
    }
    else
      output_buffer = cio->buffer;

    cout << "DiskImageResourceJP2::write(): "
	 << "writing encoded image: " << m_filename << endl;

    cout << "DiskImageResourceJP2::write(): "
	 << "encoded buffer length = " << codestream_length << endl;
    cout << "DiskImageResourceJP2::write(): "
	 << "bytes to end of encode stream = " << bytes_to_end << endl;

    // write the buffer to disk
    FILE *output_fp = fopen(m_filename.c_str(), "wb");
    if (output_fp == 0)
      vw_throw(IOErr() << "DiskImageResourceJP2::write(): "
	       "failed to open " << m_filename);

    fwrite(output_buffer, 1, codestream_length, output_fp);

    cout << "DiskImageResourceJP2::write(): "
	 << "closing streams and files" << endl;

    fclose(output_fp);

    // close and free the byte stream
    opj_cio_close(cio);

    // free remaining compression structures
    opj_destroy_compress(cinfo);

    cout << "DiskImageResourceJP2::write(): "
	 << "cleaning up" << endl;

    // if allocated output_buffer ourselves, free it here
    if(convert_to_jpx)
      delete[] output_buffer;
    // delete image structures
//     delete [] jpeg2k_image->comps;
//     delete jpeg2k_image;
  }

  // A FileIO hook to open a file for reading
  DiskImageResource*
  DiskImageResourceJP2::construct_open(std::string const& filename)
  {
    return new DiskImageResourceJP2(filename);
  }

  // A FileIO hook to open a file for writing
  DiskImageResource*
  DiskImageResourceJP2::construct_create(std::string const& filename,
					 ImageFormat const& format)
  {
    return new DiskImageResourceJP2(filename, format);
  }
}
