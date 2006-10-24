#include <vw/Image.h>
#include <vw/FileIO.h>

typedef vw::ImageView<vw::PixelRGB<float> > Image;

int main() {
  Image mural;
  read_image( mural, "mural.jpg" );
  write_image( "mural_rotate_180.jpg", rotate_180(mural) );
  write_image( "mural_flip_vertical.jpg", flip_vertical(mural) );
  write_image( "mural_flip_horizontal.jpg", flip_horizontal(mural) );
  write_image( "mural_rotate_90_cw.jpg", rotate_90_cw(mural) );
  write_image( "mural_rotate_90_ccw.jpg", rotate_90_ccw(mural) );
  write_image( "mural_transpose.jpg", transpose(mural) );
  write_image( "mural_crop.jpg", crop(mural,80,60,160,120) );
  write_image( "mural_subsample.jpg", subsample(mural,2) );
  write_image( "mural_threshold.jpg", threshold(mural,0.5) );
  write_image( "mural_clamp.jpg" , clamp(mural,0.25,0.75) );
  return 0;
}
