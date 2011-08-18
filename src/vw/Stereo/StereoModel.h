// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_STEREO_STEREOMODEL_H__
#define __VW_STEREO_STEREOMODEL_H__

#include <vw/Stereo/DisparityMap.h>

namespace vw {

// forward declaration
namespace camera {
  class CameraModel;
}

namespace stereo {

  class StereoModel {
    const camera::CameraModel *m_camera1, *m_camera2;
    bool m_least_squares;

  public:

    //------------------------------------------------------------------
    // Constructors / Destructors
    //------------------------------------------------------------------
    StereoModel(camera::CameraModel const* camera_model1,
                camera::CameraModel const* camera_model2,
                bool least_squares_refine = false) :
      m_camera1(camera_model1), m_camera2(camera_model2),
      m_least_squares(least_squares_refine) {}

    //------------------------------------------------------------------
    // Public Methods
    //------------------------------------------------------------------

    /// Apply a stereo model to a disparity map to produce an image of
    /// XYZ points.  Missing pixels in the disparity map will result
    /// in zero vector pixels in the point image.
    ///
    /// Users really shouldn't use this method, the ideal method is
    /// the 'stereo_triangulate' in StereoView.h.
    ImageView<Vector3> operator()(ImageView<PixelMask<Vector2f> > const& disparity_map,
                                  ImageView<double> &error ) const;

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
                              double& error) const;

    void refine_point( Vector2 const& pix1,
                       Vector2 const& pix2,
                       Vector3& point ) const;
  };

}}      // namespace vw::stereo

#endif  // __VW_STEREO_STEREOMODEL_H__

