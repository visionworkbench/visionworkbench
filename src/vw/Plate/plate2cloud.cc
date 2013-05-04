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


#include <limits>
#include <string>
#include <vw/Plate/PlateFile.h>
#include <vw/Plate/TileManipulation.h>
using namespace vw;
using namespace vw::platefile;

#include <boost/tokenizer.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/program_options.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
namespace po = boost::program_options;
namespace fs = boost::filesystem;
#include <kml/dom.h>
#include <kml/engine/kml_file.h>
#include <kml/base/file.h>


// plate2cloud: Query the given specified plate, create the map images
// and kml meta-data, and upload them to Google Cloud Storage using the
// Python program gsutil.

// The user needs to have gsutil installed as well as a Google Cloud
// Storage account. The account settings and other preferences for
// gsutil must be set in $HOME/.boto before calling this program.

struct Options {
  std::string output_name, planet, url_name, mod_plate_base_url;
  ssize_t last_level;
  int num_machines, machine_index;
  bool resume_job;
  std::string gs_prog; // path to gsutil
  
  Options():planet("Earth"), last_level(std::numeric_limits<ssize_t>::max()),
            num_machines(1), machine_index(0),
            resume_job(false),
            gs_prog("gsutil"){}
  //gs_prog = "$HOME/Python-2.7.2/build/bin/python $HOME/bin/gsutil-mp/gsutil";
  
};

void handle_arguments(int argc, char* argv[], Options& opt) {
  
  po::options_description general_options("Extract plate into KML tiles");
  general_options.add_options()
    ("output_name,o", po::value(&opt.output_name), "Output name for the KML result.")
    ("planet,p", po::value(&opt.planet), "Planet name (Earth/Moon/Mars).")
    ("mod_plate", po::value(&opt.mod_plate_base_url),
     "Do not create copy of images. Use mod plate images with user-defined base url (default: off).")
    ("last_level,l", po::value(&opt.last_level), "Last level to generate (default: do all levels).")
    ("num_machines,n", po::value(&opt.num_machines), "Number of machines to use (default: 1).")
    ("machine_index,i", po::value(&opt.machine_index), "Machine index, between 0 and num_machines-1 (default: 0).")
    ("resume,r", "Resume previously interrupted job (default: false).")
    ("gs_prog", "Path to gsutil.")
    ("help,h", "Display this help message.");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("url", po::value(&opt.url_name), "");

  po::options_description options("");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("url",-1);

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
    po::notify( vm );
  } catch (const po::error& e) {
    vw_throw( ArgumentErr() << "Error parsing input:\n\t"
              << e.what() << options );
  }

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " <plate url> [options]\n";

  if ( vm.count("help") || opt.url_name.empty() )
    vw_throw( ArgumentErr() << usage.str() << general_options );
  if ( vm.count("resume") ) opt.resume_job = true;
  if ( *(opt.url_name.rbegin()) == '/' )
    opt.url_name = opt.url_name.substr(0,opt.url_name.size()-1);
  if ( opt.output_name.empty() ) {
    size_t folder_divide, ext_divide;
    folder_divide = opt.url_name.rfind("/");
    ext_divide = opt.url_name.rfind(".");
    opt.output_name =
      opt.url_name.substr(folder_divide+1,ext_divide-folder_divide-1);
  }

  // Rm trailing slashes from output_name
  while (!opt.output_name.empty() && opt.output_name[opt.output_name.size() - 1] == '/'){
    opt.output_name.resize(opt.output_name.size() - 1);
  }

  // Validate the planet name
  boost::algorithm::to_lower(opt.planet);
  if ( opt.planet != "earth" && opt.planet != "moon" && opt.planet != "mars"){
    vw_throw( ArgumentErr() << "Unknown planet.\n\t" << options );
  }

  if (opt.num_machines <= 0 || opt.machine_index < 0 || opt.machine_index >= opt.num_machines){
      vw_throw( ArgumentErr() << "Invalid number of machines or machine index.\n\t" << options );
  }
  
  return;
}

std::string get_dir3(std::string pathName){
  
  // Out of the directory path: dir1/dir2/dir3/dir4/dir5/.../ return the string dir1/dir2/dir3
  // If the input path is smaller than 3 directories, such as dir1/dir2, return "".

  typedef boost::tokenizer<boost::char_separator<char> >  tokenizer;
  boost::char_separator<char> slash("/");
  tokenizer tokens(pathName, slash);
  
  // Get the last non-empty token
  std::string outDir = "";
  int count = 0, three = 3;
  for (tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter){
    if (*tok_iter == "") continue;
    outDir += *tok_iter + "/";
    count++;
    if (count >= three) break;
  }

  if (count < three) return "";

  // Strip the trailing slash
  int len = outDir.size();
  if (len >0 && outDir[len - 1] == '/'){
    outDir.resize(len - 1);
  }

  return outDir;
}

void get_resume_job_data(// Inputs
                         bool resume_job,
                         const std::string& base_folder, 
                         // Outputs
                         long long int& resume_pos,
                         std::string& resume_file
                         ){

  // Find out the position from which to resume a previously interrupted job.
  
  resume_pos = -1; // Must be negative
  
  resume_file = base_folder + "_job.txt";
  if (!resume_job) return;

  std::ifstream jfile(resume_file.c_str());
  if (! (jfile >> resume_pos) ){
    std::cout << "Warning: Could not read the resuming position from " << resume_file << std::endl;
  }
  
  std::cout << "Resuming from position " << resume_pos << std::endl;
  return;
}

kmldom::LatLonAltBoxPtr
create_latlonaltbox( double south, double west, double degree_size,
                     kmldom::KmlFactory* factory ) {
  kmldom::LatLonAltBoxPtr latlonaltbox =
    factory->CreateLatLonAltBox();
  latlonaltbox->set_north( south + degree_size );
  latlonaltbox->set_south( south );
  latlonaltbox->set_east( west + degree_size );
  latlonaltbox->set_west( west );
  return latlonaltbox;
}

kmldom::LatLonBoxPtr
create_latlonbox( double south, double west, double degree_size,
                  kmldom::KmlFactory* factory ) {
  kmldom::LatLonBoxPtr latlonbox =
    factory->CreateLatLonBox();
  latlonbox->set_north( south + degree_size );
  latlonbox->set_south( south );
  latlonbox->set_east( west + degree_size );
  latlonbox->set_west( west );
  return latlonbox;
}

kmldom::LodPtr
create_lod( ssize_t min, ssize_t max, kmldom::KmlFactory* factory ) {
  kmldom::LodPtr lod = factory->CreateLod();
  lod->set_minlodpixels( min );
  lod->set_maxlodpixels( max );
  return lod;
}

std::string prefix_from_location( ssize_t col, ssize_t row, ssize_t level ) {
  std::string result = "";
  size_t count = 0;
  for ( ssize_t i = level; i >= 0; i-- ) {
    if ( (count % 3 == 0) && ( count != 0 ) )
      result += "/";
    size_t mask = 0x01 << i;
    result += ( ( ( col & mask  ) >> i ) |
                ( ( ( row & mask ) >> i ) << 1 ) ) + 48;
    count++;
  }
  return result;
}

void create_network_links(std::string const& base_folder,
                          fs::path const& path,
                          const std::list<TileHeader> & links,
                          kmldom::KmlFactory* factory,
                          kmldom::FolderPtr & folder){
  
  // Add network links to the given kml stored in factory
  
  BOOST_FOREACH( TileHeader const& lower, links ) {
    double ldeg_delta = 360.0 / double( 1 << lower.level() );
    double llon_west = -180.0 + lower.col()*ldeg_delta;
    double llat_south = 180 - ldeg_delta*(lower.row()+1);

    kmldom::NetworkLinkPtr netlink = factory->CreateNetworkLink();
    kmldom::RegionPtr region = factory->CreateRegion();
    region->set_lod( create_lod( 128, -1, factory ) );
    region->set_latlonaltbox( create_latlonaltbox( llat_south, llon_west,
                                                   ldeg_delta, factory ) );
    netlink->set_region(region);
    fs::path link_path( base_folder + "/" +
                        prefix_from_location( lower.col(), lower.row(),
                                              lower.level() ) );
    kmldom::LinkPtr linkptr = factory->CreateLink();
    linkptr->set_viewrefreshmode(kmldom::VIEWREFRESHMODE_ONREGION);
    if ( ( lower.level() % 3 == 0 ) && lower.level() != 0 ) {
      // Pulling out the directory change
      std::string link_string = path.filename()+"/"+link_path.filename()+".kml";
      linkptr->set_href( link_string );
    } else {
      linkptr->set_href( fs::path(link_path).replace_extension(".kml").filename() );
    }
    netlink->set_link(linkptr);
    folder->add_feature( netlink );
  }

  return;
}

void get_files_in_dir(std::string gs_prog, std::string dir_name, // inputs
                      std::set<std::string> & dir_list           // outputs
                      ){

  // Query the Google Storage server and extract the list of files in
  // the requested directory as a set.  This function is used to
  // verify if all files were uploaded. This flow is rather immature.
  
  dir_list.clear();
  if (dir_name == "") return;

  string tmpFile = "gs_output.txt";
  string cmd = gs_prog + " ls -l gs://" + dir_name + "* > " + tmpFile;

  std::cout << "Executing: " << cmd << std::endl;
  for (int s = 0; s < 2; s++){
    int ans = system(cmd.c_str());
    std::cout << "success is " << ans << std::endl;
    if (ans == 0) break;
  }
  system("date");
  
  string line;
  std::ifstream fh (tmpFile.c_str());
  if (fh.is_open()){
    while ( fh.good() ){
      getline (fh,line);

      // Parse the line, extact the files
      typedef boost::tokenizer<boost::char_separator<char> >  tokenizer;
      boost::char_separator<char> space(" ");
      tokenizer tokens(line, space);
      for (tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter){
        std::string tok = *tok_iter;
        if (tok.size() >= 5 && tok.substr(0, 5) == "gs://"){
          std::string file = tok.substr(5, tok.size() - 5);
          dir_list.insert(file);
        }
      }
      
    }
    fh.close();
  }

  // Wipe the temporary file
  cmd = string("rm -f ") + tmpFile;
  system(cmd.c_str());

}

void do_upload(std::string base_folder, long long int curr_pos, std::string resume_file, std::string gs_prog){

  // Upload to google storage the current directory, and then wipe it.

  string cmd = string("sleep 1; ")
    + gs_prog + " mb gs://" + base_folder + " 2>/dev/null; " // make the bucket
    + "cd " + base_folder + "; "                             // go to the corresp. local directory
    + gs_prog + " -m cp -a public-read -R * gs://"
    + base_folder;// + ">/dev/null 2>&1";                    // copy the dir to the bucket
  std::cout << "Executing " << cmd << std::endl;
  system(cmd.c_str());

  cmd = "rm -rf " + base_folder;
  std::cout << "Executing " << cmd << std::endl;
  system(cmd.c_str());

  // Write the position at which we are now, to know where
  // to resume if this job got interrupted.
  std::ofstream jfile(resume_file.c_str());
  jfile << curr_pos << std::endl;
  jfile.close();
  
  return;
}

std::string path2basename (std::string pathName){

  // Out of path/to/mydir extract the string "mydir". Can handle strings having 
  // multiple slashes such as /path/to/dir/// and also dots, such as
  // /path/to/d.i.r

  typedef boost::tokenizer<boost::char_separator<char> >  tokenizer;
  boost::char_separator<char> slash("/");
  tokenizer tokens(pathName, slash);

  if (tokens.begin() == tokens.end()){
    vw_throw( ArgumentErr() << "Invalid directory: " << pathName );
  }

  // Get the last non-empty token
  std::string dirName = "";
  for (tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter){
    if (*tok_iter != "")  dirName = *tok_iter;
  }

  if (dirName == ""){
    vw_throw( ArgumentErr() << "Invalid directory: " << pathName );
  }

  return dirName;
  
}

class GroundOverlayEngine {
  size_t m_draw_order;

protected:
  kmldom::GroundOverlayPtr
  create_linkless_overlay( TileHeader const& tile,
                           kmldom::KmlFactory* factory,
                           bool lowest_overlay ) const {
    double deg_delta = 360.0 / double( 1 << tile.level() );
    double lon_west = -180.0 + tile.col()*deg_delta;
    double lat_south = 180 - deg_delta*(tile.row()+1);

    kmldom::GroundOverlayPtr goverlay = factory->CreateGroundOverlay();
    kmldom::RegionPtr region = factory->CreateRegion();
    if ( tile.level() == 1 ) {
      // Top layerish. Layer 0 is never drawn because it's
      // boundaries are illegal (Lat range too large). So this step
      // is to make sure layer 1 can always be seen when zoomed out.
      region->set_lod( create_lod( 1, 513, factory ) );
    } else if ( lowest_overlay ) {
      // End of branch
      region->set_lod( create_lod( 128, -1, factory ) );
    } else {
      // Original code went to 1024
      region->set_lod( create_lod( 128, 513, factory ) );
    }
    region->set_latlonaltbox( create_latlonaltbox( lat_south, lon_west,
                                                   deg_delta, factory ) );
    goverlay->set_region( region );
    goverlay->set_latlonbox( create_latlonbox( lat_south, lon_west,
                                               deg_delta, factory ) );
    goverlay->set_draworder( m_draw_order );
    return goverlay;
  }
public:
  GroundOverlayEngine( size_t draw_order = 50 ) : m_draw_order(draw_order){}
};

class LocalEquiEngine : public GroundOverlayEngine {
  boost::shared_ptr<PlateFile> m_platefile;
public:
  LocalEquiEngine( boost::shared_ptr<PlateFile> platefile,
                   size_t draw_order = 50 ) : GroundOverlayEngine(draw_order),
                                              m_platefile(platefile) {}

  kmldom::GroundOverlayPtr
  operator()( fs::path const& kml_location, TileHeader const& tile,
              kmldom::KmlFactory* factory, bool lowest_overlay, bool do_write ) const {

    string imgType;
    if (do_write){
      // Write to disk and return the extension of the image file name written
      std::pair<std::string, TileHeader> returnVal
        = m_platefile->read_to_file(kml_location.string(), tile.col(), tile.row(),
                                    tile.level(), tile.transaction_id(), true);
      imgType = returnVal.second.filetype();
    }else{
      // Return the extension only, without writing to disk
      boost::tuple<std::string, TileHeader, TileData>
        ret = m_platefile->img_file_name(kml_location.string(), tile.col(), tile.row(),
                                         tile.level(), tile.transaction_id(), true);
      imgType = ret.get<1>().filetype(); 
    }

    // The image has exactly the same name as the kml file, but with a .jpg or .png extension
    // instead of .kml.
    string imgName  = fs::path(kml_location).replace_extension(imgType).filename().string();

    kmldom::GroundOverlayPtr goverlay = this->create_linkless_overlay( tile, factory, lowest_overlay );
    kmldom::IconPtr icon = factory->CreateIcon();
    icon->set_href(imgName);

    goverlay->set_icon(icon);
    return goverlay;
  }
};

class ModPlateEquiEngine : public GroundOverlayEngine {
  std::string m_base_url, m_extension;
public:
  ModPlateEquiEngine( std::string const& base_url,
                      std::string const& extension = "jpg",
                      size_t draw_order = 50 ) :  GroundOverlayEngine(draw_order),
                                                  m_base_url(base_url),
                                                  m_extension(extension) {}

  kmldom::GroundOverlayPtr
  operator()( fs::path /*kml_location*/, TileHeader const& tile,
              kmldom::KmlFactory* factory, bool lowest_overlay, bool do_write ) const {
    kmldom::GroundOverlayPtr goverlay =
      this->create_linkless_overlay( tile, factory, lowest_overlay );
    kmldom::IconPtr icon = factory->CreateIcon();
    std::ostringstream ostr;
    ostr << m_base_url << "/" << tile.level() << "/" << tile.row()
         << "/" << tile.col() << "." << m_extension;
    icon->set_href(ostr.str());
    goverlay->set_icon( icon );
    return goverlay;
  }
};

template <class EngineT>
void draw_kml_level ( std::string const& base_folder,
                      std::string const& planet,
                      ssize_t level, ssize_t last_level,
                      int num_machines, int machine_index,
                      long long int resume_pos, std::string resume_file,
                      BBox2i const& region, TransactionOrNeg id,
                      boost::shared_ptr<PlateFile> plate,
                      EngineT const& goverlay_engine,
                      std::string & dir3,
                      std::set<std::string> & dir3_list,
                      std::string gs_prog
                      ) {

  if (level > last_level){
    return;
  }

  // Upload to the cloud the current level of kml, and call itself recursively for 
  // higher levels.

  // The flag do_write is used for debugging. It can be set to false to simulate just
  // quering the plate file system, without writing any of the output fetched from it.
  bool do_write = true;
      
  // Turning on the flag do_missing will make the code work in the mode where
  // it does verification of whether all files were submitted to the cloud,
  // and resubmits the ones which were not. This flow is rather immature.
  bool do_missing = false;
  if (do_missing) do_write = false;    
  
  int start = 1; // Don't modify this lightly.
  static long long int kml_count = start - 1; // Will equal 'start' at 1st iter
  static long long int dir_count = 0;
  
  // Do not delete this pointer. LibKML is managing this itself somewhere.
  kmldom::KmlFactory* factory = kmldom::KmlFactory::GetFactory();

  // Searching for transactions
  std::list<TileHeader> tile_records;
  tile_records = plate->search_by_region(level,region,TransactionRange(id));

  // Drawing Level
  BOOST_FOREACH( TileHeader const& t, tile_records ) {

    std::list<TileHeader> links =
      plate->search_by_region(level+1,region*2,TransactionRange(id));
    
    if (level == last_level){
      // Put no links to the next level if we intend the current
      // level to be the last
      links.clear(); 
    }

    fs::path path( base_folder + "/" +
                   prefix_from_location( t.col(), t.row(), t.level() ) );

    fs::path curr_dir = path.parent_path();

    // If dir3 of the current directory changed, we must
    // recompute the list of files in dir3.
    std::string curr_dir3 = "";
    if (do_missing){
      curr_dir3 = get_dir3(curr_dir.string());
      if (curr_dir3 != dir3){
        dir3 = curr_dir3;
        if (dir3 != ""){
          dir_count++;
          if (dir_count % num_machines == machine_index){
            get_files_in_dir(gs_prog, dir3, // inputs
                             dir3_list      // outputs
                             );
          }
        }
      }
    }
    
    if (do_write || do_missing) fs::create_directories(curr_dir);

    // Producing this layer's kml
    kmldom::FolderPtr folder = factory->CreateFolder();
    folder->set_name( path.string() );
    create_network_links(base_folder, path, links, factory, folder); // child links

    kml_count++;
    bool willDo = false;
    if (do_missing){
      willDo =  ( (dir3 != "") && (dir_count % num_machines == machine_index) );
    }else{
      willDo = ( (kml_count > resume_pos) && (kml_count % num_machines == machine_index) );
    }
      
    if (willDo){ // curr machine

      // Start kml writing

      // Writing the image to disk. Creating the ground overlay.
      kmldom::GroundOverlayPtr goverlay = goverlay_engine( path, t, factory, !links.size(), do_write );
      folder->add_feature(goverlay);
      string img_file = curr_dir.string() + "/" + goverlay->get_icon()->get_href();

      if ( do_missing && dir3 != "" ){
        if ( dir3_list.find(img_file) != dir3_list.end() ){
          //std::cout << "Found " << img_file << " in " << dir3 << std::endl;
        }else{
          std::cout << "Missing file: " << img_file << " in " << dir3 << std::endl;
          // Write the missing file
          bool l_do_write = true;
          kmldom::GroundOverlayPtr goverlay = goverlay_engine( path, t, factory, !links.size(), l_do_write );
        }
      }
        
      // Writing the kml to disk.
      kmldom::KmlPtr kml = factory->CreateKml();
      kml->set_feature( folder );
      kml->set_hint("target=" + planet);
      kmlengine::KmlFilePtr kml_ptr = kmlengine::KmlFile::CreateFromImport(kml);
      std::string xml_data;
      kml_ptr->SerializeToString(&xml_data);
      std::string kml_file = fs::path(path).replace_extension(".kml").string();
      if  (kml_count == start){
        // Rename 0.kml to be the opening kml
        kml_file = (fs::path(base_folder) / (path2basename(base_folder) + ".kml")).string();
      }

      bool is_file_missing = false;
      if ( do_missing && dir3 != "" ){
        if ( dir3_list.find(kml_file) != dir3_list.end() ){
          //std::cout << "Found " << kml_file << " in " << dir3 << std::endl;
        }else{
          is_file_missing = true;
          std::cout << "Missing file: " << kml_file << " in " << dir3 << std::endl;
        }
      }
        
      if (do_write || is_file_missing){
        //std::cout << "Writing " << kml_file << std::endl;
        kmlbase::File::WriteStringToFile(xml_data, kml_file);
      }
      // end kml writing

      int period = 20000; // upload after this many kml files have been created
      if ( (kml_count/num_machines) % period == 0 ){
        do_upload(base_folder, kml_count, resume_file, gs_prog);
      }
    }

    // Recursive call
    BOOST_FOREACH( TileHeader const& lower, links ) {
      draw_kml_level( base_folder, planet, lower.level(), last_level,
                      num_machines, machine_index, resume_pos, resume_file,
                      BBox2i( lower.col(), lower.row(), 1, 1 ),
                      id, plate, goverlay_engine,
                      dir3, dir3_list,
                      gs_prog
                      );
    }
  }
  return;
}

int main( int argc, char *argv[] ) {

  Options opt;
  try {
    handle_arguments( argc, argv, opt );

    long long int resume_pos;
    string resume_file;
    get_resume_job_data(opt.resume_job, opt.output_name, // inputs
                        resume_pos, resume_file          // outputs
                        );

    boost::shared_ptr<PlateFile> platefile( new PlateFile(opt.url_name) );

    // Create out output directory and chdir
    fs::create_directory( opt.output_name );

    std::string dir3 = "";
    std::set<std::string> dir3_list;
    
    if ( opt.mod_plate_base_url.empty() ) {
      // Standard local platefile
      LocalEquiEngine engine( platefile );
      draw_kml_level( opt.output_name, opt.planet, 0, opt.last_level,
                      opt.num_machines, opt.machine_index, 
                      resume_pos, resume_file,
                      BBox2i(0,0,1,1),
                      TransactionOrNeg(-1), platefile, engine,
                      dir3, dir3_list, opt.gs_prog
                      );
    } else {
      // Recursive call that goes depth first in generation.
      ModPlateEquiEngine engine( opt.mod_plate_base_url );
      draw_kml_level( opt.output_name, opt.planet, 0, opt.last_level,
                      opt.num_machines, opt.machine_index, 
                      resume_pos, resume_file,
                      BBox2i(0,0,1,1),
                      TransactionOrNeg(-1), platefile, engine,
                      dir3, dir3_list, opt.gs_prog
                      );
    }
 
    // Upload the last batch of files
    do_upload(opt.output_name, std::numeric_limits<long long int>::max(), resume_file, opt.gs_prog); 

  } catch ( const ArgumentErr& e ) {
    vw_out() << e.what() << std::endl;
    return 1;
  } catch ( const Exception& e ) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
