/// deplate.cc
///
/// Converts a plate file to a nested set of directories and tiles on
/// disk.

#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Plate/PlateFile.h>

using namespace vw;
using namespace vw::platefile;

int main( int argc, char *argv[] ) {
  
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <platefile name> <directory name>";
    exit(0);
  }

  // Open the plate file
  PlateFile platefile(argv[1]);

  // Create the output directory
  std::string directory_name = argv[2];
  if ( !fs::exists(directory_name) )
    fs::create_directory(directory_name);

  int nlevels = 4;
  for (int n = 0; n < nlevels; ++n) {
    int block_cols = pow(2,n);
    int block_rows = pow(2,n);
    
    for (int j = 0; j < block_rows; ++j) {
      for (int i = 0; i < block_cols; ++i) {
        
        ImageView<PixelRGB<uint8> > tile;
        try {
          IndexRecord rec = platefile.read(tile, i, j, n);

          if (rec.valid()) {
          
            std::cout << "\t--> Extracting: [ " << i << " " << j << " " <<  n << "]\n";
            
            std::ostringstream ostr;
            ostr << directory_name << "/" << n;
            if ( !fs::exists(ostr.str()) )
              fs::create_directory(ostr.str());
            
            ostr << "/" << i;
            if ( !fs::exists(ostr.str()) )
              fs::create_directory(ostr.str());
            
            ostr << "/" << j << "." << rec.block_filetype();
            std::cout << "Writing " << ostr.str() << "\n";
            write_image(ostr.str(), tile);

          }
          
        } catch (TileNotFoundErr &e) {
          // Do nothing for tiles that aren't found.
        }
      }
    }    

    platefile.save();

  }
}
