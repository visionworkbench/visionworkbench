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


#ifndef __VW_STEREO_STEREOMODEL_H__
#define __VW_STEREO_STEREOMODEL_H__

#include <vw/Math/Vector.h>

namespace vw {

  template <class PixelT> class ImageView;
  template <class PixelT> struct PixelMask;

  // forward declaration
  namespace camera {
    class CameraModel;
  }

namespace stereo {

  // Two-ray triangulation. Triangulate the point by finding the
  // midpoint of the segment joining the closest points on the two
  // rays emanating from the camera.
  vw::Vector3 triangulate_pair(vw::Vector3 const& dir0, vw::Vector3 const& ctr0, 
                               vw::Vector3 const& dir1, vw::Vector3 const& ctr1, 
                               vw::Vector3& errorVec);
  
  /// Class for triangulating rays emanating at given pixels in given cameras 
  class StereoModel {

  public:

    //------------------------------------------------------------------
    // Constructors / Destructors
    //------------------------------------------------------------------
    // Constructor with n cameras
    StereoModel(std::vector<const camera::CameraModel *> const& cameras,
                double angle_tol = 0.0);
    // Constructor with two cameras
    StereoModel(camera::CameraModel const* camera_model1,
                camera::CameraModel const* camera_model2,
                double angle_tol = 0.0);

    // This class does not modify the input models.
    virtual ~StereoModel() {}

    //------------------------------------------------------------------
    // Public methods
    //------------------------------------------------------------------

    /// Apply a stereo model to a disparity map to produce an image of
    /// XYZ points.  Missing pixels in the disparity map will result
    /// in zero vector pixels in the point image.
    ///
    /// Users really shouldn't use this method, the ideal method is
    /// the 'stereo_triangulate' in StereoView.h.
    ImageView<Vector3> operator()(ImageView<PixelMask<Vector2f>> const& disparity_map,
                                  ImageView<double> &error) const;

    /// Apply a stereo model to multiple or just two image coordinates.
    /// Returns an xyz point. The error is set to 0 if triangulation
    /// did not succeed, otherwise it is the vector between the closest points on the rays.
    virtual Vector3 operator()(std::vector<Vector2> const& pixVec,
                               Vector3& errorVec) const;
    virtual Vector3 operator()(std::vector<Vector2> const& pixVec,
                               double & error) const;
    virtual Vector3 operator()(Vector2 const& pix1, Vector2 const& pix2,
                               Vector3& errorVec) const;
    virtual Vector3 operator()(Vector2 const& pix1, Vector2 const& pix2,
                               double & error) const;

    /// Returns the dot product of the two rays emanating from camera
    /// 1 and camera 2 through pix1 and pix2 respectively.  This can
    /// effectively be interpreted as the angle (in radians) between
    /// the two rays, and it is a useful test for checking when the
    /// two rays are close to parallel.
    double convergence_angle(Vector2 const& pix1, Vector2 const& pix2) const;

    // The default 1-cos(x) function does badly when x close to 0,
    // which then leads to incorrect angle tolerance.
    static double robust_1_minus_cos(double x);

  protected:

    std::vector<const camera::CameraModel *> m_cameras;
    double m_angle_tol;

    //------------------------------------------------------------------
    // Protected methods
    //------------------------------------------------------------------

    /// Return the 2-norm of the error vector ( the vector from the
    /// closest point of intersection of A to the closest point of
    /// intersection of B ), or -1 if the rays are parallel or divergent.
    static Vector3 triangulate_point(std::vector<Vector3> const& camDirs,
                                     std::vector<Vector3> const& camCtrs,
                                     Vector3& errorVec);
    
    static bool are_nearly_parallel(double angle_tol,
                                    std::vector<Vector3> const& camDirs);
  };

}}      // namespace vw::stereo

#endif  // __VW_STEREO_STEREOMODEL_H__

