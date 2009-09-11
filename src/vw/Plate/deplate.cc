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
  platefile.print();

  // Create the output directory
  std::string directory_name = argv[2];
  if ( !fs::exists(directory_name) )
    fs::create_directory(directory_name);

  // The number of levels is hard-coded for now.
  int nlevels = 4;
  for (int n = 0; n < nlevels; ++n) {
    int block_cols = pow(2,n);
    int block_rows = pow(2,n);
    
    for (int j = 0; j < block_rows; ++j) {
      for (int i = 0; i < block_cols; ++i) {
        
        ImageView<PixelRGB<uint8> > tile;
        try {
          std::cout << "\t--> Extracting: [ " << j << " " << i << " @ level " <<  n << "]\n";
          IndexRecord rec = platefile.read(tile, i, j, n);
          if (!rec.valid())
            vw_throw(TileNotFoundErr() << "\tTile was found, but was marked invalid.");
          std::cout << "done.\n\n\n";

          
          
          // Create the level directory (if it doesn't exist)
          std::ostringstream ostr;
          ostr << directory_name << "/" << n;
          if ( !fs::exists(ostr.str()) )
            fs::create_directory(ostr.str());
          
          // Create the column directory (if it doesn't exist)
          ostr << "/" << i;
          if ( !fs::exists(ostr.str()) )
            fs::create_directory(ostr.str());
          
          // Create the file (with the row as the filename)
          ostr << "/" << j << "." << rec.block_filetype();
          std::cout << "Writing " << ostr.str() << "\n";
          write_image(ostr.str(), tile);

        } catch (TileNotFoundErr &e) {
          std::cout << e.what() << "\n";
        }
      }
    }    

    platefile.save();

  }
}
