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


#include <iostream>
#include <vw/Image.h>
#include <vw/FileIO.h>

int main( int argc, char *argv[] ) {
  try {
    VW_ASSERT( argc==3, vw::ArgumentErr() << "Invalid command-line args." );
    vw::ImageView<vw::PixelRGBA<float> > image;
    read_image( image, argv[1] );
    write_image( argv[2], image );
  }
  catch( const vw::Exception& e ) {
    std::cerr << "Error: " << e.what() << std::endl;
    std::cerr << "Usage: vwconvert <source> <destination>" << std::endl;
    return 1;
  }
  return 0;
}
