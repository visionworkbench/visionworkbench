// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Image.h>
#include <vw/Plate/PlateFile.h>
#include <vw/Plate/TileManipulation.h>
using namespace vw;
using namespace vw::platefile;

#include <boost/foreach.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/math/distributions/fisher_f.hpp>

using namespace std;

// --- Functions to Apply --------------------------

// Reduce Base class
//   defines interface
template <typename ImplT>
struct ReduceBase {
  inline ImplT& impl() { return static_cast<ImplT&>(*this); }
  inline ImplT const& impl() const { return static_cast<ImplT const&>(*this); }

  template <class PixelT>
  inline void operator()( list<ImageView<PixelT> > const& input,
                          list<TileHeader> const& input_header,
                          ImageView<PixelT> & output) {
    impl()(input,input_header,output);
  }
};

// Weighted Average Implementation
//   mimick this with own functions
struct WeightedAverage : public ReduceBase<WeightedAverage> {

  template <class PixelT>
  inline void operator()( list<ImageView<PixelT> > const& input,
                          list<TileHeader> const& /*input_header*/,
                          ImageView<PixelT> & output) {
    // Input Images will always have an alpha channel. That is a requirement of PlateFiles.
    int num_channels = CompoundNumChannels<PixelT>::value;
    std::vector<ImageView<float32> > sum_weighted_data;
    for ( uint8 i = 0; i < num_channels-1; i++ )
      sum_weighted_data.push_back(ImageView<float32>(input.front().cols(),
                                                     input.front().rows()));
    ImageView<float32> summed_weights(input.front().cols(),
                                      input.front().rows());

    // Summing multiple images
    typedef typename list< ImageView<PixelT > >::const_iterator image_iter;
    for ( image_iter image = input.begin();
          image != input.end(); image++ ) {
      summed_weights += channel_cast<float32>(select_channel(*image,num_channels-1));
      // Iterating over non alpha channels
      for ( uint8 i = 0; i < num_channels-1; i++ )
        sum_weighted_data[i] += channel_cast<float32>(select_channel(*image,num_channels-1))*channel_cast<float32>(select_channel(*image,i));
    }

    // Normalizing
    for ( uint8 i = 0; i < num_channels-1; i++ )
      sum_weighted_data[i] /= summed_weights;
    output.set_size(input.front().cols(),
                    input.front().rows());
    for ( uint8 i = 0; i < num_channels-1; i++ )
      select_channel(output,i) = sum_weighted_data[i];
    // Setting output alpha
    select_channel(output,num_channels-1) = threshold(summed_weights,0,
                                                      ChannelRange<typename PixelChannelType<PixelT>::type>::min(),
                                                      ChannelRange<typename PixelChannelType<PixelT>::type>::max());
  }
};

// Weighted Variance Implementation
struct WeightedVar2 : public ReduceBase<WeightedVar2> {

  template <class PixelT>
  inline void operator()( list<ImageView<PixelT> > const& input,
                          list<TileHeader> const& /*input_header*/,
                          ImageView<PixelT> & var2) {
    // Input Images will always have an alpha channel. That is a requirement of PlateFiles.
    int num_channels = CompoundNumChannels<PixelT>::value;
    std::vector<ImageView<float32> > sum_weighted_data; //store the mean
    std::vector<ImageView<float32> > sum_weighted_data2; //store the variance

    for ( uint8 i = 0; i < num_channels-1; i++ ){
      sum_weighted_data.push_back(ImageView<float32>(input.front().cols(),
                                                     input.front().rows()));
      sum_weighted_data2.push_back(ImageView<float32>(input.front().cols(),
                                                     input.front().rows()));
    }

    ImageView<float32> summed_weights(input.front().cols(),
                                      input.front().rows());

    // Summing multiple images
    typedef typename list< ImageView<PixelT > >::const_iterator image_iter;
    for ( image_iter image = input.begin();image != input.end(); image++ ) {
      summed_weights += channel_cast<float32>(select_channel(*image,num_channels-1));

      // Iterating over non alpha channels
      for ( uint8 i = 0; i < num_channels-1; i++ ){
        sum_weighted_data[i] += channel_cast<float32>(select_channel(*image,num_channels-1))*channel_cast<float32>(select_channel(*image,i));

        sum_weighted_data2[i] += channel_cast<float32>(select_channel(*image,num_channels-1))
                                 *channel_cast<float32>(select_channel(*image,i))*channel_cast<float32>(select_channel(*image,i));
      }
    }

    // Normalizing
    for ( uint8 i = 0; i < num_channels-1; i++ ){
      sum_weighted_data[i] /= summed_weights;
      sum_weighted_data2[i] = sum_weighted_data2[i]/summed_weights - sum_weighted_data[i];
    }

    var2.set_size(input.front().cols(),
                  input.front().rows());
    for ( uint8 i = 0; i < num_channels-1; i++ )
      select_channel(var2,i) = sum_weighted_data2[i];
    // Setting output alpha
    select_channel(var2,num_channels-1) =
      threshold(summed_weights,0,
                ChannelRange<typename PixelChannelType<PixelT>::type>::min(),
                ChannelRange<typename PixelChannelType<PixelT>::type>::max());
  }
};

// Robust Mean Implementation
struct RobustMean : public ReduceBase<RobustMean> {
private:
  float smart_weighted_mean( vector<float> & weights,
                             vector<float> const& samples,
                             float const sign_level=0.3,
                             float const learn_rate=0.1,
                             float const error_tol=1e-5,
                             int32 const max_iter=1000 ) {
    namespace bm = boost::math;
    switch ( samples.size() ) {
    case 1:
      weights[0]=1;
      return samples[0];
    case 2:
      weights[0]=1;
      weights[1]=1;
      return (samples[0]+samples[1])/2;
    default:

      float weighted_mean = 0;
      std::vector<float> prev_wt = weights; // the previous weights
      for (int i = 0; i < max_iter; ++i) {

        // accumulation of all sums
        float SW1 = 0;    // sum of weights
        float SW2 = 0;    // sum of squared weights
        float SWX = 0;    // weighted sum of samples
        float WX2 = 0;    // weighted sum of squared data
        for (size_t j = 0; j < samples.size(); ++j) {
          SW1 += weights[j];
          SW2 += weights[j]*weights[j];
          SWX += weights[j]*samples[j];
          WX2 += weights[j]*samples[j]*samples[j];
        }

        // weighted mean
        weighted_mean = SWX/SW1;

        float DN = SW1*WX2-SWX*SWX; // denominator
        float v2 = SW1-SW2/SW1;     // the second degree of freedom
        float sse_wt = 0;           // squared sum of weight differences
        for (size_t j = 0; j < samples.size(); ++j) {

          // F statistic
          float DF = samples[j]-weighted_mean; // difference from the mean
          float FS = (SW1*SW1-SW2)*DF*DF/DN;   // F statistic

          // degree of freedom
          float v1 = 1-weights[j]/SW1;         // the first degree of freedom

          // p-value calculation
          float p_value = 1;                   // p-value
          if (SW1 > 1 && DN > 0) {             // basic assumptions
            bm::fisher_f dist(v1,v2);
            p_value = 1-cdf(dist, FS);
          }

          // gradient decent update
          weights[j] += learn_rate*(p_value-sign_level);
          if ( weights[j] < 0 ) weights[j] = 0; // 0.0 is lower bound of weight
          if ( weights[j] > 1 ) weights[j] = 1; // 1.0 is upper bound of weight

          // squared sum of weight differences
          sse_wt += (weights[j]-prev_wt[j])*(weights[j]-prev_wt[j]);
          prev_wt[j] = weights[j];              // update previous weight

        }

        // terminal condition: mean squared difference of weights is
        // smaller than the tolerance
        if ( sse_wt/weights.size() < error_tol ) break;
      }

      return  weighted_mean;
    }
  }

public:
  template <class PixelT>
  inline void operator()( list<ImageView<PixelT> > const& input,
                          list<TileHeader> const& /*input_header*/,
                          ImageView<PixelT> & output ) {
    uint8 num_v_channels = CompoundNumChannels<PixelT>::value-1;
    output.set_size( input.front().cols(),
                     input.front().rows() );

    // Iterating through every pixel within a tile
    for ( int32 ix = 0; ix < input.front().cols(); ix++ ) {
      for ( int32 iy = 0; iy < input.front().rows(); iy++ ) {
        vector<std::vector<float> > t_samples;
        t_samples.resize(num_v_channels);
        vector<float> t_weight;

        // Seperating data out into vectors, only loading up data that
        // is not empty.
        typedef typename list<ImageView<PixelT> >::const_iterator iter_type;
        for ( iter_type iter = input.begin();
              iter != input.end(); iter++ ) {
          if ( (*iter)(ix,iy)[num_v_channels] == 0 )
            continue;
          t_weight.push_back( float((*iter)(ix,iy)[num_v_channels])/float(ChannelRange<typename PixelChannelType<PixelT>::type>::max()) );
          for ( uint8 iz = 0; iz < num_v_channels; iz++ )
            t_samples[iz].push_back( (*iter)(ix,iy)[iz] );
        }

        if ( t_weight.empty() ) {
          output(ix,iy) = PixelT();
          continue;
        }

        typedef typename PixelChannelType<PixelT>::type ChannelT;
        output(ix,iy)[num_v_channels] = ChannelRange<ChannelT>::max();
        for ( uint8 iz = 0; iz < num_v_channels; iz++ ) {
          vector<float> copy = t_weight;
          output(ix,iy)[iz] =
            boost::numeric_cast<ChannelT>(smart_weighted_mean( copy,
                                                               t_samples[iz]));
        }
      } // iy - end loop
    }   // ix - first loop
  }     // end operator()
};

// --- Standard Terminal Argument ------------------

// Standard Arguments
struct Options {
  // Input
  string url;
  int32 level;
  int32 start_trans_id;
  int32 end_trans_id;

  // Output
  string function;
  int32 transaction_id;

  // For spawning multiple jobs
  int32 job_id;
  int32 num_jobs;
};

void handle_arguments(int argc, char *argv[], Options& opt) {
  po::options_description general_options("Perform weighted averages of all layers within a tile inside a plate file");
  general_options.add_options()
    ("job_id,j", po::value(&opt.job_id)->default_value(0), "")
    ("num_jobs,n", po::value(&opt.num_jobs)->default_value(1), "")
    ("start_t", po::value(&opt.start_trans_id)->default_value(0), "Input starting transaction ID range.")
    ("end_t", po::value(&opt.end_trans_id), "Input ending transaction ID range.")
    ("level,l", po::value(&opt.level)->default_value(-1), "Level inside the plate in which to process. -1 will error out and show the number of levels available.")
    ("function,f", po::value(&opt.function)->default_value("WeightedAvg"), "Functions that are available are [WeightedAvg RobustMean WeightedVar]")
    ("transaction-id,t",po::value(&opt.transaction_id)->default_value(2000), "Transaction id to write to")
    ("help", "Display this help message");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("input-file", po::value<std::string>(&opt.url), "");

  po::options_description options("");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("input-file",-1);

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
    po::notify( vm );
  } catch (po::error &e) {
    vw_throw( ArgumentErr() << "Error parsing input:\n\t"
              << e.what() << options );
  }

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " <plate_filename> [options]\n";

  if ( vm.count("help") || opt.url.empty() )
    vw_throw( ArgumentErr() << usage.str() << options );
}

// --- Meta Application of Above Functions ----------

// apply_reduce
template <typename ReduceT, class PixelT>
void apply_reduce( boost::shared_ptr<PlateFile> platefile,
                   std::list<BBox2i> const& workunits,
                   Options& opt, ReduceBase<ReduceT>& reduce) {

  TerminalProgressCallback tpc("plate.platereduce", "Processing");
  double inc_tpc = 1.0/float(workunits.size());
  BOOST_FOREACH( const BBox2i& workunit, workunits) {
    tpc.report_incremental_progress(inc_tpc);
    for ( int ix = 0; ix < workunit.width(); ix++ ) {
      for ( int iy = 0; iy < workunit.height(); iy++ ) {
        Vector2i location(ix,iy);
        location += workunit.min();

        // Polling for Tiles
        std::list<TileHeader> tile_records;
        tile_records = platefile->search_by_location(location[0],
                                                     location[1],
                                                     opt.level,
                                                     opt.start_trans_id,
                                                     opt.end_trans_id, true);

        // No Tiles? No Problem!
        if (tile_records.empty())
          continue;

        // Loading images
        std::list<ImageView<PixelT> > tiles;
        BOOST_FOREACH( const TileHeader& tile, tile_records ) {
          ImageView<PixelT> new_tile;
          platefile->read( new_tile, location[0],
                           location[1], opt.level,
                           tile.transaction_id(), true );
          tiles.push_back(new_tile);
        }

        // Calling function
        ImageView<PixelT> result;
        reduce(tiles, tile_records, result);

        platefile->write_request();
        platefile->write_update(result,
                                location[0], location[1],
                                opt.level, opt.transaction_id);
        platefile->write_complete();
      }
    }
  }
  tpc.report_finished();
}

// Function that runs the apply_reduce over the plate file
template <typename ReduceT>
void do_run( Options& opt, ReduceBase<ReduceT>& reduce ) {
  boost::shared_ptr<PlateFile> platefile =
    boost::shared_ptr<PlateFile>( new PlateFile(opt.url) );

  if ( opt.level < 0 ||
       opt.level >= platefile->num_levels() ) {
    vw_throw( ArgumentErr() << "In correct level selection, "
              << opt.level << ".\n\nPlatefile " << opt.url << " has "
              << platefile->num_levels() << " levels internally.\n" );
  }

  // This is arbitrary, just needed to divide up jobs
  int32 region_size = 1 << opt.level;
  BBox2i full_region(0,0,region_size,region_size);
  std::list<BBox2i> workunits = bbox_tiles(full_region,4,4);
  std::list<BBox2i> mworkunits;
  int32 count = 0;
  BOOST_FOREACH(const BBox2i& c, workunits) {
    if (count==opt.num_jobs)
      count=0;
    if (count==opt.job_id)
      mworkunits.push_back(c);
    count++;
  }
  vw_out() << "Job " << opt.job_id << "/" << opt.num_jobs << " has "
           << mworkunits.size() << " work units.\n";

  switch(platefile->pixel_format()) {
  case VW_PIXEL_GRAYA:
    switch(platefile->channel_type()) {
    case VW_CHANNEL_UINT8:
      apply_reduce<ReduceT, PixelGrayA<uint8> >(platefile, mworkunits,
                                                opt, reduce);
      break;
    case VW_CHANNEL_INT16:
      apply_reduce<ReduceT, PixelGrayA<int16> >(platefile, mworkunits,
                                                opt, reduce);
      break;
    case VW_CHANNEL_FLOAT32:
      apply_reduce<ReduceT, PixelGrayA<float32> >(platefile, mworkunits,
                                                  opt, reduce);
      break;
    default:
      vw_throw(InputErr() << "Platefile contains unsupported channel type.\n" );
    }
    break;
  case VW_PIXEL_RGBA:
    switch(platefile->channel_type()) {
    case VW_CHANNEL_UINT8:
      apply_reduce<ReduceT, PixelRGBA<uint8> >(platefile, mworkunits,
                                               opt, reduce);
      break;
    default:
      vw_throw(InputErr() << "Platefile contains unsupported channel type.\n" );
    }
    break;
  default:
    vw_throw(InputErr() << "Platefile contains a pixel type thats unsupported.\n" );
  }
}

int main( int argc, char *argv[] ) {

  Options opt;
  try {
    handle_arguments( argc, argv, opt );

    // Handing out jobs now
    boost::to_lower(opt.function);
    if ( opt.function == "weightedavg" ) {
      WeightedAverage f;
      do_run<WeightedAverage>( opt, f );
    } else if ( opt.function == "robustmean" ) {
      RobustMean f;
      do_run<RobustMean>( opt, f );
    }
    else if ( opt.function == "weightedvar" ) {
      WeightedVar2 f;
      do_run<WeightedVar2>( opt, f );
    } else {
      vw_throw( ArgumentErr() << "Unknown function, " << opt.function << "\n" );
    }

  } catch ( const ArgumentErr& e ) {
    vw_out() << e.what() << std::endl;
    return 1;
  } catch ( const Exception& e ) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
