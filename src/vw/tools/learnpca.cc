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
#include <vw/InterestPoint/LearnPCA.h>

int main(int argc, char *argv[])
{
  if (argc < 2) {
    std::cout << "learnpca <training image files...>"
              << std::endl;
    return 0;
  }

  LearnPCA lpca("pca_basis.exr", "pca_avg.exr");

  for (int i = 1; i < argc; i++) {
    vw::DiskImageView<vw::PixelRGB<vw::uint8> > disk_image(argv[i]);
    lpca.processImage(disk_image);
  }
  lpca.runPCA();

  return 0;
}

