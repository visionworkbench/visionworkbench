#ifndef __VW_STEREO_STEREOMODEL_H__
#define __VW_STEREO_STEREOMODEL_H__

#include <vw/Camera/CameraModel.h>
#include <vw/Stereo/DisparityMap.h>

#include <vector>

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
    //    ImageView<Vector3> operator()(ImageView<PixelDisparity<double> > const& disparity_map);

    /// Apply a stereo model to a single pair of image coordinates.
    /// Returns an xyz point.  The error is set to -1 if the rays were
    /// parallel or divergent, otherwise it returns the 2-norm of the
    /// distance between the rays at their nearest point of
    /// intersection.
    inline Vector3 operator()(Vector2 const& pix1, Vector2 const& pix2, double& error ) {
      
      try {
        // determine range by triangulation
        Vector3 originA = m_camera1->camera_center(pix1);
        Vector3 vecFromA = m_camera1->pixel_to_vector(pix1);
        
        Vector3 originB = m_camera2->camera_center(pix2);
        Vector3 vecFromB = m_camera2->pixel_to_vector(pix2);

        Vector3 result =  triangulate_point(originA, vecFromA,
                                            originB, vecFromB, 
                                            error);

        if (pix1.y() >= 5643 && pix1.y() <= 5650) {
          std::cout << "Pix 1: " << pix1 << "    pix2: " << pix2 << "\n";
          std::cout << "\tA:  " << originA << "     " << vecFromA << "\n";
          std::cout << "\tB:  " << originB << "     " << vecFromB << "\n";
          std::cout << "\tResult:  " << result << "    " << error << "\n\n";
        }

        return result;
      } catch (vw::camera::PixelToRayErr &e) {
        error = 0;
        return Vector3();
      }
    }
  
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

#endif	// __VW_ORBITINGPUSHBROOM_STEREOMODEL_H__

