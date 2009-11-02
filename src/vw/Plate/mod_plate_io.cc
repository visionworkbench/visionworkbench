
// Vision Workbench
#include <vw/Plate/RemoteIndex.h>
#include <vw/Plate/PlateFile.h>

// Apache and the Apache Runtime
#include "httpd.h"
#include "http_config.h"
#include "http_protocol.h"
#include "ap_config.h"
#include "http_core.h"

// STL
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// Boost
#include <boost/algorithm/string.hpp>

using namespace vw;
using namespace vw::platefile;

#define MODPLATE_DEBUG 1


class PlateModule {


  // ---------------------------------------------------
  //                 Content Handlers
  // ---------------------------------------------------

  int error_handler(request_rec *r, std::string error_message) {
    r->content_type = "text/html";   
    if (MODPLATE_DEBUG) {
      ap_rprintf(r, "%s", error_message.c_str());
      return OK;
    } else {
      return HTTP_NOT_FOUND; 
    }
  }
  
  int image_handler(request_rec *r, int platefile_id,
                    int col, int row, int depth, std::string file_suffix) {

    // --------------  Access Plate Index -----------------

    IndexRecord idx_record;
    std::string plate_filename;
    try {
      //platefile::Index idx(plate_filename);
      vw_out(0) << "------------ CREATING REMOTE INDEX ----------------\n";
      platefile::RemoteIndex idx(platefile_id, "foobar");
      // if (!m_index)
      //   return error_handler(r, "Error: connection is not open to remote index server.");
      vw_out(0) << "------------ ISSUING READ REQUEST ----------------\n";
      idx_record = idx.read_request(col,row,depth);

      vw_out(0) << "------------ ISSUING INFO REQUEST ----------------\n";
      plate_filename = idx.platefile_name();
      
      //      For debugging:
      // r->content_type = "text/html";      
      // ap_rprintf(r, "<p>Platefile name: %s.</p>\n", plate_filename.c_str());
      // // ap_rprintf(r, "<p>Index Max depth: %d.</p>\n", idx.version());
      // ap_rprintf(r, "<p>Index Record (%d):</p>", idx_record.status()); 
      // ap_rprintf(r, "<p>Index Record (%d):</p>", idx_record.status() != INDEX_RECORD_VALID); 
      // ap_rprintf(r, "<p>  Blob ID: %d</p>", idx_record.blob_id()); 
      // ap_rprintf(r, "<p>  Blob Offset: %d</p>", int(idx_record.blob_offset())); 
      // ap_rprintf(r, "<p>  Tile Size: %d</p>", idx_record.tile_size()); 
      // //       ap_rprintf(r, "<p>  Blob File Type: %s</p>", idx_record.tile_filetype().c_str());
      // return OK;
    
    } catch(vw::Exception &e) {
      std::ostringstream ostr;
      ostr << "An error occured reading the plate index: " << e.what() << "\n";
      return error_handler(r, ostr.str());
    }  

    // ---------------- Return the image ------------------

    if (idx_record.status() != INDEX_RECORD_VALID)
      return error_handler(r, "Error: the index record was not valid.");
    
    try {
      // Compute the blob filename and open the plate index.  TODO: This
      // should really be opened in some sort of plugin initialization
      // function rather than on each query!!
      std::ostringstream ostr;
      ostr << plate_filename << "/plate_" << idx_record.blob_id() << ".blob";
      std::string blob_filename = ostr.str();

      // Open the blob file and read the data.
      Blob blob(blob_filename, true);
      boost::shared_array<uint8> data = blob.read_data(idx_record.blob_offset());

      // Read in the data and send it over the web connection.  TODO:
      // This is probably not as fast as memory mapping the blob file,
      // so this could probably be sped up.
      r->content_type = "image/png";      
      ap_rwrite(data.get(), idx_record.tile_size(), r);
      return OK;
    } catch (Exception &e) {
      return error_handler(r, "An unknown error occured.");
    }

  }

  int wtml_handler(request_rec *r, std::string filename) {
    r->content_type = "text/plain";      

    ap_rprintf(r, "<?xml version='1.0' encoding='UTF-8'?>\n");
    ap_rprintf(r, "<Folder Name=\"ARC Test Data\" Group=\"View\">\n");

    ap_rprintf(r, "  <ImageSet Generic=\"False\" DataSetType=\"Earth\" BandPass=\"Visible\" Name=\"Bluemarble\" Url=\"http://localhost:31337/plate_example/{1}/{2}/{3}.png\" BaseTileLevel=\"0\" TileLevels=\"2\" BaseDegreesPerTile=\"360\" FileType=\".png\" BottomsUp=\"False\" Projection=\"Toast\" QuadTreeMap=\"0123\" CenterX=\"0\" CenterY=\"0\" OffsetX=\"0\" OffsetY=\"0\" Rotation=\"0\" Sparse=\"False\" ElevationModel=\"False\" StockSet=\"False\">\n");
    ap_rprintf(r, "    <ThumbnailUrl>http://localhost:31337/plate_example/0/0/0.png</ThumbnailUrl>\n");
    ap_rprintf(r, "    <Credits>NASA</Credits>\n");
    ap_rprintf(r, "  </ImageSet>\n");
    ap_rprintf(r, "</Folder>\n");
    return OK;
  }

  // Private variables
  boost::shared_ptr<std::ofstream> m_output_stream;
  //  boost::shared_ptr<platefile::RemoteIndex> m_index;
  

public: 
  PlateModule() {

    // Redirect the output stream to a file.
    m_output_stream.reset( new std::ofstream("/tmp/mod_plate.log") );
    vw_log().set_console_stream(*(m_output_stream.get()));
    vw_out(0) << "\n\nInitializing mod_plate module.\n";

  }

  ~PlateModule() {}

  void apache_module_init() {
  }

  void apache_module_destroy() {
  }

  int apache_module_callback(request_rec *r) {
    if (r->header_only) return OK;

    // First, we do a little bit of trimming to get rid of any
    // front-slashes that are leading or trailing the path_info.
    std::string path_info = r->path_info;
    if (path_info.size() == 0) 
      return DECLINED;

    int start = 0, end = path_info.size();
    if (path_info[start] == '/')
      ++start;
    std::string trimmed_path_info = path_info.substr(start, end-start);

    // --------------  Parse the URL String -----------------

    // Split the path into tokens
    std::vector<std::string> path_tokens;
    boost::split( path_tokens, trimmed_path_info, boost::is_any_of("/") );

    // r->content_type = "text/plain";      
    // ap_rprintf(r, "path info: %d\n", path_tokens.size());
    // return OK;

    // Image Query [platefile_id, level, col, row]
    if (path_tokens.size() == 4 && path_tokens[path_tokens.size()-1].size() != 0) {        
      vw_out(0) << "Handling IMAGE request: " << r->path_info << "\n";
    
      // Split the final element into [ row, file_suffix ]
      std::vector<std::string> final_tokens;
      boost::split( final_tokens, 
                    path_tokens[path_tokens.size()-1], boost::is_any_of(".") );
      if (final_tokens.size() != 2) 
        return error_handler(r, "Error parsing file type.  Expected an image filename.");

      // Parse the tokens
      int platefile_id = atoi(path_tokens[0].c_str());
      int depth = atoi(path_tokens[1].c_str());
      int col = atoi(path_tokens[2].c_str());
      int row = atoi(final_tokens[0].c_str());
      std::string file_suffix = final_tokens[1];

      return image_handler(r, platefile_id, col, row, depth, file_suffix);
    
      // WTML query
    } else if (path_tokens.size() == 1 && path_tokens[0].size() != 0) {
      vw_out(0) << "Handling WTML request: " << r->path_info << "\n";

      // Split the element into [ name, file_suffix ]
      std::vector<std::string> final_tokens;
      boost::split( final_tokens, 
                    path_tokens[0], boost::is_any_of(".") );
      if (final_tokens.size() != 2 || final_tokens[1] != "wtml")
        return error_handler(r, "Error parsing file type. Expected a query for WTML file.\n");

      return wtml_handler(r, path_tokens[0]);

    } else {                              // Invalid query
      return error_handler(r, "Path must have 3 or 1 tokens.\n");
    }

  }

};

// ---------------------------------------------------
// Create a single instance of the SystemLog
// ---------------------------------------------------
namespace {
  vw::RunOnce mod_plate_once = VW_RUNONCE_INIT;
  boost::shared_ptr<PlateModule> mod_plate_ptr;
  void init_mod_plate() {
    mod_plate_ptr = boost::shared_ptr<PlateModule>(new PlateModule());
  }
}

/// Static method to access the singleton instance of the plate module object.
PlateModule& mod_plate() {
  mod_plate_once.run( init_mod_plate );
  return *mod_plate_ptr;
}

// --------------------- Apache C++ Entry Points ------------------------

extern "C" void mod_plate_init() {
  mod_plate().apache_module_init();
}

extern "C" void mod_plate_destroy() {
  mod_plate().apache_module_destroy();
}

extern "C" int mod_plate_callback(request_rec *r) {
  mod_plate().apache_module_callback(r);
}



