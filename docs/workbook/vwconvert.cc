// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <iostream>
#include <vw/Image.h>
#include <vw/FileIO.h>

int main( int argc, char *argv[] ) {
  try {
    VW_ASSERT( argc==3, vw::ArgumentErr() << "Invalid command-line args." );
    vw::ImageView<float> image;
    read_image( image, argv[1] );
    write_image( argv[2], image );
  }
  catch( vw::Exception& e ) {
    std::cerr << "Error: " << e.what() << std::endl;
    std::cerr << "Usage: vwconvert <source> <destination>" << std::endl;
    return 1;
  }
  return 0;
}
