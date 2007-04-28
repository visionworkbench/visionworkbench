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
#include <vw/Math/Matrix.h>
#include <vw/Core/Exception.h>
#include <vw/Cartography/DiskImageResourceGeoReferenceHelperJP2.h>
#include <vw/FileIO/DiskImageResourceJP2.h>
#include <vw/FileIO/JP2.h>
#include <vw/Cartography/GeoReference.h>


// Reference: GML in jpx: http://www.opengeospatial.org/standards/gmljp2


namespace {

  vw::JP2SuperBox* make_gmljp2_boxes(vw::uint8* d, vw::uint64 nbytes, bool allocated)
  {
    vw::JP2DataBox* b;
    vw::JP2SuperBox* a;
    vw::JP2SuperBox* a2;
    //vw::JP2SuperBox::JP2BoxIterator pos;
    char* tmpstr;
    
    // find Contiguous Codestream box
    //pos.reset();
    //b = (JP2DataBox*)find_box(0x6A703263, &pos); // "jp2c"
    //if(!b)
    //  return -1;

    // add outer Association box before Contiguous Codestream box
    a = new vw::JP2SuperBox(0x61736F63 /*"asoc"*/);
    //insert_box_before(pos, a);
    //NOTE: a will be deleted from sub_boxes list

    // add Label box with "gml.data" to outer Association box
    tmpstr = new char[9];
    strcpy(tmpstr, "gml.data");
    b = new vw::JP2DataBox(0x6C626C20 /*"lbl\040"*/, (vw::uint8*)tmpstr, 8, true);
    //NOTE: tmpstr will be deleted by JP2DataBox
    a->insert_box_last(b);
    //NOTE: b will be deleted by JP2SuperBox

    // add inner Association box to outer Association box
    a2 = new vw::JP2SuperBox(0x61736F63 /*"asoc"*/);
    a->insert_box_last(a2);
    //NOTE: a2 will be deleted by JP2SuperBox

    // add Label box with "gml.root-instance" to inner Association box
    tmpstr = new char[18];
    strcpy(tmpstr, "gml.root-instance");
    b = new vw::JP2DataBox(0x6C626C20 /*"lbl\040"*/, (vw::uint8*)tmpstr, 17, true);
    //NOTE: tmpstr will be deleted by JP2DataBox
    a2->insert_box_last(b);
    //NOTE: b will be deleted by JP2SuperBox

    // add XML box with GML root-instance data to inner Association box
    b = new vw::JP2DataBox(0x786D6C20 /*"xml\040"*/, d, nbytes, allocated);
    //NOTE: d will be deleted by JP2DataBox iff allocated
    a2->insert_box_last(b);
    //NOTE: b will be deleted by JP2SuperBox

    return a;
  }

}

namespace vw {
namespace cartography {

  /*static*/ void DiskImageResourceGeoReferenceHelperJP2::read_georeference( FileMetadata* georef_, DiskImageResource* r_ ) {
    DiskImageResourceJP2* r = (DiskImageResourceJP2*)r_;
    GeoReference* georef = (GeoReference*)georef_;
    JP2BoxList* additional_boxes = (JP2BoxList*)r->additional_boxes();
    JP2ReaderRequirementsList* additional_requirements = (JP2ReaderRequirementsList*)r->additional_requirements();
    //FIXME
  }
  
  /*static*/ void DiskImageResourceGeoReferenceHelperJP2::write_georeference( DiskImageResource* r_, FileMetadata const* georef_ ) {
    DiskImageResourceJP2* r = (DiskImageResourceJP2*)r_;
    GeoReference const* georef = (GeoReference const*)georef_;
    JP2BoxList* additional_boxes = (JP2BoxList*)r->additional_boxes();
    JP2ReaderRequirementsList* additional_requirements = (JP2ReaderRequirementsList*)r->additional_requirements();
    
    JP2SuperBox* a = make_gmljp2_boxes( (uint8*)(georef->gml_str().c_str()), georef->gml_str().length(), false );
    additional_boxes->push_back( a );
    // 67 is the (standard) requirement number for GMLJP2
    additional_requirements->push_back( std::make_pair( 67, false ) );
  }

}} // namespace vw::cartography
