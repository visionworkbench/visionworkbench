#include <iostream>
#include "learnpca.h"

int main(int argc, char *argv[])
{

  DiskImageView<PixelRGB<uint8> > disk_image(argv[1]);
  LearnPCA lpca("pca_basis.exr", "pca_avg.exr");

  lpca.processImage(disk_image);
  lpca.runPCA();

  return 0;
}

