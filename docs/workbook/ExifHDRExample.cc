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


#include <vw/vw.h>
#include <stdio.h>
#include <vector>
using namespace vw;
using namespace vw::HDR;

int main(int argc, char** argv) {
  std::vector<Vector<double> > curves;
  std::vector<string> files;
  for(int i = 1; i < argc; i++)
    files.push_back(argv[i]);
  // Process HDR stack using Exif tags
  ImageView<PixelRGB<double> > hdr_exif = process_ldr_images_exif(files,
                                          curves);
  write_image("hdr.exr", hdr_exif);

  // Apply Drago tone-mapping operator.
  ImageView<PixelRGB<double> > tone_mapped = drago_tone_map(hdr_exif);
  write_image("tm.jpg", tone_mapped);

  // Apply gamma correction and save.
  ImageView<PixelRGB<double> > gamma = pow(tone_mapped, 1.0/2.2);
  write_image("tm_gamma.jpg", gamma);

  // Re-apply camera response curves and save.
  // First must invert curves calculated earlier.
  std::vector<Vector<double> > inverse_curves(curves.size());
  for (int i = 0; i < curves.size(); i++) {
    invert_curve(curves[i], inverse_curves[i],
                 VW_HDR_RESPONSE_POLYNOMIAL_ORDER);
  }
  psi(tone_mapped, inverse_curves);
  write_image("tm_curved.jpg", tone_mapped);

  // Apply gamma correction after response curves.
  // Usually gives best results.
  ImageView<PixelRGB<double> > tm_c_g = pow(tone_mapped, 1.0/2.2);
  write_image("tm_c_g.jpg", tm_c_g);  

  return 0;
}
