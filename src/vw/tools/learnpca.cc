// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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

