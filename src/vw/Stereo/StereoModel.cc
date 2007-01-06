#include <vw/Stereo/StereoModel.h>
#include <vw/Math/Vector.h>

using namespace std;

namespace vw {
namespace stereo {

  ImageView<Vector3> StereoModel::operator()(ImageView<PixelDisparity<double> > const& disparity_map,
                                             ImageView<double> &error) {

    // Error analysis
    double mean_error = 0.0;
    double max_error = 0.0;
    int point_count = 0;
    int divergent = 0;
    
    // Allocate xyz image and get pointer to buffer
    ImageView<Vector3> xyz(disparity_map.cols(), disparity_map.rows());
    error.set_size(disparity_map.cols(), disparity_map.rows());
        
    // Compute 3D position for each pixel in the disparity map
    cout << "StereoModel: Applying camera models\n";
    for (unsigned int y = 0; y < disparity_map.rows(); y++) {
      if (y % 100 == 0) {
        printf("\tStereoModel computing points: %0.2f%% complete.\r", 100.0f*float(y)/disparity_map.rows());
        fflush(stdout);      
      }
      for (unsigned int x = 0; x < disparity_map.cols(); x++) {
        if ( !disparity_map(x,y).missing() ) {
          xyz(x,y) = (*this)(Vector2( x, y ),
                             Vector2( x + disparity_map(x,y).h(),
                                      y + disparity_map(x,y).v() ), 
                             error(x,y) );          

          if (error(x,y) >= 0) {
            // Keep track of error statistics
            if (error(x,y) > max_error)
              max_error = error(x,y);
            mean_error += error(x,y);
            ++point_count;
          } else {	   
            // rays diverge or are parallel
            xyz(x,y) = Vector3();
            divergent++;
          }
        } else {
          xyz(x,y) = Vector3();
          error(x,y) = 0;
        }
      }
    }  
    
    if (divergent != 0) 
      cout << "WARNING in StereoModel: " << divergent << " rays diverged or were parallel!\n";
     
    cout << "\tStereoModel computing points: Done.                  \n";
    cout << "\tMean error = " << mean_error/double(point_count)
         << ",  Max error = " << max_error << endl;
    return xyz;
  }
    
}} // namespace vw::stereo
