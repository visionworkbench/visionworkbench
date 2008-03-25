#include <iostream>
#include "learnpca.h"

int main(int argc, char *argv[])
{
  LearnPCA lpca("pca_basis.exr", "pca_avg.exr");

  for (int i = 1; i < argc; i++) {
    DiskImageView<PixelRGB<uint8> > disk_image(argv[i]);
    lpca.processImage(disk_image);
  }
  lpca.runPCA();

  return 0;
}

