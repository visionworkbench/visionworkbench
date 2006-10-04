#include <vw/Stereo/StereoModel.h>
#include <vw/Math/Vector.h>

using namespace std;

namespace vw {
namespace stereo {

  Vector3 StereoModel::ray_ray_closest_intersection(Vector3 const& pointA,
                                                    Vector3 const& vecFromA,
                                                    Vector3 const& pointB,
                                                    Vector3 const& vecFromB,
                                                    double& error) const {

    double det, a, b, c, d, e, s, t, sNumerator, tNumerator;
    
    Vector3 vecBToA = pointA - pointB;
    // This routine is based on pseudo code from Schneider's "Geometric
    // Tools for Computer Graphics" book
    a = 1.0;				   // assumes dir vecs are unit vecs
    b = vw::math::dot_prod(vecFromA, vecFromB);
    c = 1.0;				   // assumes dir vecs are unit vecs
    d = vw::math::dot_prod(vecFromA, vecBToA);
    e = vw::math::dot_prod(vecFromB, vecBToA);
    
    det = a*c - b*b;
    
    if (det == 0.0)	{		   // the rays are parallel
      error = -1.0;
      return Vector3();
    }
    
    sNumerator = b*e - c*d;
    tNumerator = a*e - b*d;
    
    if ((sNumerator < 0) || (tNumerator < 0)) { // the rays diverge
      error = -1.0;
      return Vector3();
    }
    
    // Then the values of line parameters of the closest points are:
    s = sNumerator / det;
    t = tNumerator / det;
    
    // Take the midpoint of the closest points as the intersection
    // point of the rays
    Vector3 closestPointA = pointA + (s * vecFromA);
    Vector3 closestPointB = pointB + (t * vecFromB);
    Vector3 errorVec = closestPointA - closestPointB;
    error = norm_2(errorVec);
//     printf("ERROR: %e\n", error);
//     std::cout << closestPointA << "\n";
//     std::cout << closestPointB << "\n";
//     std::cout << errorVec << "\n";
    return 0.5 * (closestPointA + closestPointB);
  }

  Vector3 StereoModel::triangulate_point_test(Vector3 const& pointA,
                                              Vector3 const& vecFromA,
                                              Vector3 const& pointB,
                                              Vector3 const& vecFromB,
                                              double& error) const {

    Vector3 vecAToB = pointA -pointB;
    double a, b, c, d, e, f, s, t, sNumerator, tNumerator;
    double den;
    
    // This routine is based on pseudo code from Schneider's "Geometric
    // Tools for Computer Graphics" book... however, the pseudo code there
    // is incorrect based on their derivation (note the sign of b). 
    // It has been corrected below.
    
    // Basically, you have pointA + s*vecFromA and pointB + t*vecFromB
    // and you want to solve for s and t such that you find the minimum
    // distance (ideally, the intersection) between A and B
    a = dot_prod(vecFromA, vecFromA);		
    b = dot_prod(vecFromA, vecFromB);
    c = dot_prod(vecFromB, vecFromB);	
    d = dot_prod(vecFromA, vecAToB);
    e = dot_prod(vecFromB, vecAToB);
    f = dot_prod(vecAToB, vecAToB);

    den = b*b - a*c;
    
    if (den == 0.0)	{		   // the rays are parallel
      printf("StereoModel: triangulate_point: Err1: det == 0 \n");
      error = -1.0;
      return Vector3();
    }

    sNumerator = d*c - b*e;
    tNumerator = b*d - a*e;
        
    // In the future, it would be good to add a check to see if the
    // rays diverge. This has to be done by looking at the signs of 
    // sNumerator and tNumerator for each quadrant in the image 
    // coordinate frame. 

    //  For example: (check correctness and syntax)
    //  if vecFromA(X,Y) > 0 and vecFromB (X,Y) > 0 
    //  if ((sNumerator < 0) || (tNumerator < 0)) // the rays diverge
    //  printf("StereoModel_CAHV: TriangulatePoint: Err2: rays diverge \n");
    //  return;

    // The scaling parameters at the location of the closest point between the
    // two vectors are:
    s = sNumerator / den;
    t = tNumerator / den;
    
    // Use the location between the closest points as the intersection point
    // of the rays  
    
    Vector3 closestPointA = pointA + (s * vecFromA);
    Vector3 closestPointB = pointB + (t * vecFromB);
    Vector3 errorVec = closestPointA - closestPointB;
    error = norm_2(errorVec);
    
    return 0.5 * (closestPointA + closestPointB);
  }

  ImageView<Vector3> StereoModel::operator()(ImageView<PixelDisparity<double> > const& disparity_map,
                                             ImageView<double> &error) {

    error.set_size(disparity_map.cols(), disparity_map.rows());
    
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
      printf("\tStereoModel computing points: %0.2f%% complete.\r", 100.0f*float(y)/disparity_map.rows());
      fflush(stdout);      
      for (unsigned int x = 0; x < disparity_map.cols(); x++) {
        if ( !disparity_map(x,y).missing() ) {
          Vector3 closest_point = (*this)(Vector2( x, y ),
                                          Vector2( x + disparity_map(x,y).h(),
                                                   y + disparity_map(x,y).v() ), 
                                          error(x,y) );
          
          xyz(x, y) = closest_point;
          if (error(x,y) > 0) {
            // This can serve as a new form of outlier rejection
            if (error(x,y) > max_error)
              max_error = error(x,y);
            mean_error += error(x,y);
            ++point_count;
          } else {	   
            // rays diverge or are parallel
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
