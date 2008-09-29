#include <iostream>
#include <vw/InterestPoint/LearnPCA.h>

using namespace std;

int main(int argc, char *argv[])
{
  if (argc < 2) {
    cout << "learnpca <training image files...>" << endl;
    return 0;
  }

  LearnPCA lpca("pca_basis.exr", "pca_avg.exr");

  for (int i = 1; i < argc; i++) {
    DiskImageView<PixelRGB<uint8> > disk_image(argv[i]);
    lpca.processImage(disk_image);
  }
  lpca.runPCA();

  return 0;
}

