#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Mosaic.h>

int main() {

  // Generate the variations on the mural image
  vw::ImageView<vw::PixelRGB<vw::float32> > mural;
  read_image( mural, "images/mural.jpg" );
  write_image( "images/mural_rotate_180.jpg", rotate_180(mural) );
  write_image( "images/mural_flip_vertical.jpg", flip_vertical(mural) );
  write_image( "images/mural_flip_horizontal.jpg", flip_horizontal(mural) );
  write_image( "images/mural_rotate_90_cw.jpg", rotate_90_cw(mural) );
  write_image( "images/mural_rotate_90_ccw.jpg", rotate_90_ccw(mural) );
  write_image( "images/mural_transpose.jpg", transpose(mural) );
  write_image( "images/mural_crop.jpg", crop(mural,80,60,160,120) );
  write_image( "images/mural_subsample.jpg", subsample(mural,2) );
  write_image( "images/mural_threshold.jpg", threshold(mural,0.5) );
  write_image( "images/mural_clamp.jpg" , clamp(mural,0.25,0.75) );

  // Generate the hand-lips composite
  vw::mosaic::ImageComposite<vw::PixelRGBA<vw::float32> > composite;
  vw::ImageView<vw::PixelRGBA<vw::float32> > hand, lips;
  read_image( hand, "images/hand.png" );
  read_image( lips, "images/lips.png" );
  composite.insert( hand, 0, 0 );
  composite.insert( lips, 131, 302 );
  composite.prepare();
  write_image( "images/hand-lips-blend.jpg", composite );
  // Composite source images onto white for printing
  hand += 1.0 - select_channel(hand,3);
  write_image( "images/hand.jpg", hand );
  lips += 1.0 - select_channel(lips,3);
  write_image( "images/lips.jpg", lips );
    
  // Generate the walker quadtree
  vw::ImageView<vw::PixelRGBA<vw::float32> > walker;
  read_image( walker, "images/Walker.jpg" );
  vw::mosaic::ImageQuadTreeGenerator<vw::PixelRGBA<vw::float32> > qtree( "images/Walker", walker );
  qtree.generate();
  // Composite three patches onto white for printing
  vw::ImageView<vw::PixelRGBA<vw::float32> > patch;
  read_image( patch, "images/Walker.qtree/Walker.png" );
  patch += 1.0 - select_channel(patch,3);
  write_image( "images/Walker.qtree/Walker.jpg", patch );
  read_image( patch, "images/Walker.qtree/Walker/1.png" );
  patch += 1.0 - select_channel(patch,3);
  write_image( "images/Walker.qtree/Walker/1.jpg", patch );
  read_image( patch, "images/Walker.qtree/Walker/1/2.png" );
  patch += 1.0 - select_channel(patch,3);
  write_image( "images/Walker.qtree/Walker/1/2.jpg", patch );

  return 0;
}
