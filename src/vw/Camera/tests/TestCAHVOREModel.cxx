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


#include <gtest/gtest_VW.h>

#include <vw/Math/Vector.h>
#include <vw/Camera/CAHVOREModel.h>
#include <vw/Camera/CAHVModel.h>
#include <test/Helpers.h>

using namespace vw;
using namespace vw::camera;
using namespace vw::test;

TEST( CAHVOREModel, ForwardReverse ) {
  CAHVOREModel cahvore(Vector3(0.606185,-0.043367,-0.234891),
                       Vector3(0.712013,0.037316,0.701174),
                       Vector3(353.341,474.873,350.82),
                       Vector3(44.0102,16.904,683.916),
                       Vector3(0.712953,0.038186,0.700171),
                       Vector3(3e-06,-0.013032,-0.00754),
                       Vector3(0.000942,0.00228,0.001613),
                       3, 0.37 );

  for ( uint32 i = 100; i < 901; i += 100 ) {
    for ( uint32 j = 100; j < 901; j+= 100 ) {
      Vector2 pixel( i,j );
      Vector3 unit = cahvore.pixel_to_vector( pixel );
      Vector3 point = cahvore.C + 30*unit;
      Vector2 result = cahvore.point_to_pixel( point );
      EXPECT_VECTOR_NEAR( result, pixel, 3e-2 );
    }
  }
}

TEST( CAHVOREModel, StringWriteRead ) {
  CAHVOREModel cahvore(Vector3(0.606185,-0.043367,-0.234891),
                       Vector3(0.712013,0.037316,0.701174),
                       Vector3(353.341,474.873,350.82),
                       Vector3(44.0102,16.904,683.916),
                       Vector3(0.712953,0.038186,0.700171),
                       Vector3(3e-06,-0.013032,-0.00754),
                       Vector3(0.000942,0.00228,0.001613),
                       3, 0.37 );
  UnlinkName file("CAHVORETemp.txt");
  cahvore.write( file );

  CAHVOREModel cahvore2( file );

  EXPECT_VECTOR_NEAR( cahvore2.C, cahvore.C, 1e-5 );
  EXPECT_VECTOR_NEAR( cahvore2.A, cahvore.A, 1e-5 );
  EXPECT_VECTOR_NEAR( cahvore2.H, cahvore.H, 1e-5 );
  EXPECT_VECTOR_NEAR( cahvore2.V, cahvore.V, 1e-5 );
  EXPECT_VECTOR_NEAR( cahvore2.O, cahvore.O, 1e-5 );
  EXPECT_VECTOR_NEAR( cahvore2.R, cahvore.R, 1e-5 );
  EXPECT_VECTOR_NEAR( cahvore2.E, cahvore.E, 1e-5 );
  EXPECT_NEAR( cahvore2.P, cahvore.P, 1e-5 );
}

TEST( CAHVOREModel, Linearize ) {
  CAHVOREModel cahvore(Vector3(0.606185,-0.043367,-0.234891),
                       Vector3(0.712013,0.037316,0.701174),
                       Vector3(353.341,474.873,350.82),
                       Vector3(44.0102,16.904,683.916),
                       Vector3(0.712953,0.038186,0.700171),
                       Vector3(3e-06,-0.013032,-0.00754),
                       Vector3(0.000942,0.00228,0.001613),
                       3, 0.37 );
  CAHVModel cahv =
    linearize_camera(cahvore, Vector2i(1024,1024),
                     Vector2i(1024,1024));

  // This is not the best test, but it insures consistent performance
  // while I compress the code.
  EXPECT_VECTOR_NEAR(Vector3(0.606185,-0.043367,-0.234891),
                     cahv.C, 1e-5);
  EXPECT_VECTOR_NEAR(Vector3(0.712254,0.0290096,0.701322),
                     cahv.A, 1e-5);
  EXPECT_VECTOR_NEAR(Vector3(358.009,279.762,354.175),
                     cahv.H, 1e-2);
  EXPECT_VECTOR_NEAR(Vector3(177.463,13.6499,548.543),
                     cahv.V, 1e-2);
}
