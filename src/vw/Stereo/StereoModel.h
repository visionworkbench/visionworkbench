#ifndef __VW_STEREO_STEREOMODEL_H__
#define __VW_STEREO_STEREOMODEL_H__

#include <vw/Camera/CameraModel.h>
#include <vw/Stereo/DisparityMap.h>

namespace vw { 
namespace stereo {

  class StereoModel {

  public:
    
    //------------------------------------------------------------------
    // Constructors / Destructors
    //------------------------------------------------------------------	
    StereoModel(vw::camera::CameraModel const& camera_model1, 
                vw::camera::CameraModel const& camera_model2) : 
      m_camera1(&camera_model1), m_camera2(&camera_model2) {}

    //------------------------------------------------------------------
    // Public Methods
    //------------------------------------------------------------------

    /// Apply a stereo model to a disparity map to produce an image of
    /// XYZ points.  Missing pixels in the disparity map will result
    /// in zero vector pixels in the point image.
    ImageView<Vector3> operator()(ImageView<PixelDisparity<double> > const& disparity_map,
                                  ImageView<double> &error );

    /// Apply a stereo model to a single pair of image coordinates.
    /// Returns an xyz point.  The error is set to -1 if the rays were
    /// parallel or divergent, otherwise it returns the 2-norm of the
    /// distance between the rays at their nearest point of
    /// intersection.
    Vector3 operator()(Vector2 const& pix1, Vector2 const& pix2, double& error ) const;

    /// Returns the dot product of the two rays emanating from camera
    /// 1 and camera 2 through pix1 and pix2 respectively.  This can
    /// effectively be interpreted as the angle (in radians) between
    /// the two rays, and it is a useful test for checking when the
    /// two rays are close to parallel.
    double convergence_angle(Vector2 const& pix1, Vector2 const& pix2) const;
  
  protected:
    //------------------------------------------------------------------
    // Protected Methods
    //------------------------------------------------------------------

    /// Return the 2-norm of the error vector ( the vector from the
    /// closest point of intersectio of A to the closest point of
    /// intersection of B ), or -1 if the rays are parallel or
    /// divergent.
    Vector3 triangulate_point(Vector3 const& pointA,
                              Vector3 const& vecFromA,
                              Vector3 const& pointB,
                              Vector3 const& vecFromB,
                              double& error) const {

      Vector3 v12 = cross_prod(vecFromA, vecFromB);
      Vector3 v1 = cross_prod(v12, vecFromA);
      Vector3 v2 = cross_prod(v12, vecFromB);

      Vector3 closestPointA = pointA + dot_prod(v2, pointB-pointA)/dot_prod(v2, vecFromA)*vecFromA;
      Vector3 closestPointB = pointB + dot_prod(v1, pointA-pointB)/dot_prod(v1, vecFromB)*vecFromB;

      Vector3 errorVec = closestPointA - closestPointB;
      error = norm_2(errorVec);
    
      return 0.5 * (closestPointA + closestPointB);
    }

  private:

    //------------------------------------------------------------------
    // Internal Variables
    //------------------------------------------------------------------
    const vw::camera::CameraModel*		m_camera1;
    const vw::camera::CameraModel*		m_camera2;
  };

}}	// namespace vw::stereo

#endif	// __VW_STEREO_STEREOMODEL_H__

