// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


//
// Copyright 2006 Carnegie Mellon University. All rights reserved.
//
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
//
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
//
// __END_LICENSE__

// TestPinholeModelCalibrate.h
#include <cxxtest/TestSuite.h>

#include <vw/Camera/PinholeModelCalibrate.h>
#include <vw/Math.h>

#include <cstdlib>
#include <ctime>

using namespace vw;
using namespace vw::camera;

class TestPinholeModelCalibrate : public CxxTest::TestSuite
{
  double mean_error(const PinholeModel& m, const std::vector<vw::Vector2>& pixels, const std::vector<vw::Vector3>& points) {
    double mean = 0;
    for (int i = 0; i < pixels.size(); i++)
      mean += vw::math::norm_2(pixels[i] - m.point_to_pixel(points[i]));
    return mean / pixels.size();
  }

  double mean_error(const PinholeModel& m, const std::vector<vw::Vector2>& pixels, const std::vector<vw::Vector3>& points, const std::vector<int>& indices) {
    double mean = 0;
    for (int j = 0; j < indices.size(); j++) 
      mean += vw::math::norm_2(pixels[indices[j]] - m.point_to_pixel(points[indices[j]]));
    return mean / indices.size();
  }


  inline double rand(double max = 100) {
    return max*std::rand()/RAND_MAX;
  }

  unsigned long sequence;
  double repeatable_srand(unsigned long init = 10) {
    sequence = init;
  }

  double repeatable_rand(double max = 100) {
    sequence = sequence*1103515245 + 12345;
    return max * static_cast<double>(sequence) / std::numeric_limits<unsigned long>::max();
  }

public:
  void test_pinholemodel_serialize_four() {
    std::srand(std::time(0));
    
    for (int i = 0; i < 10; i++) {
      Vector3 cc(rand(), rand(), rand());
      Vector3 rv(rand(), rand(), rand());
      Matrix3x3 rm(vw::math::axis_angle_to_matrix(rv));
      double fu = rand(), fv = rand(), cu = rand(), cv = rand();
      Vector3 u(rand(), rand(), rand());
      Vector3 v(rand(), rand(), rand());
      Vector3 w(rand(), rand(), rand());
      Vector4 tsaiv(vw::Vector4(rand(), rand(), rand(), rand()));
      TsaiLensDistortion tsai(tsaiv);
      PinholeModel m(cc, rm, fu, fv, cu, cv, u, v, w, tsai);

      vw::Vector<double> serial(serialize_pinholemodel<PinholeModelSerializeIntrinsic, PinholeModelSerializeRotation, PinholeModelSerializeTranslation, PinholeModelSerializeTSAI>(m));

      // check if correctly serialized
      TS_ASSERT_EQUALS(serial(0), fu);
      TS_ASSERT_EQUALS(serial(1), fv);
      TS_ASSERT_EQUALS(serial(2), cu);
      TS_ASSERT_EQUALS(serial(3), cv);
      Vector3 model_rv(m.camera_pose().axis_angle());
      TS_ASSERT_EQUALS(serial(4), model_rv.x());
      TS_ASSERT_EQUALS(serial(5), model_rv.y());
      TS_ASSERT_EQUALS(serial(6), model_rv.z());
      Vector3 model_cc(m.camera_center());
      TS_ASSERT_EQUALS(serial(7), model_cc.x());
      TS_ASSERT_EQUALS(serial(8), model_cc.y());
      TS_ASSERT_EQUALS(serial(9), model_cc.z());
      TS_ASSERT_EQUALS(serial(10), tsaiv(0));
      TS_ASSERT_EQUALS(serial(11), tsaiv(1));
      TS_ASSERT_EQUALS(serial(12), tsaiv(2));
      TS_ASSERT_EQUALS(serial(13), tsaiv(3));

      // test deserialization
      PinholeModel d;
      d.set_coordinate_frame(u, v, w);
      deserialize_pinholemodel<PinholeModelSerializeIntrinsic, PinholeModelSerializeRotation, PinholeModelSerializeTranslation, PinholeModelSerializeTSAI>(d, serial);
      double d_fu, d_fv, d_cu, d_cv;
      d.intrinsic_parameters(d_fu, d_fv, d_cu, d_cv);
      TS_ASSERT_EQUALS(fu, d_fu);
      TS_ASSERT_EQUALS(fv, d_fv);
      TS_ASSERT_EQUALS(cu, d_cu);
      TS_ASSERT_EQUALS(cv, d_cv);
      TS_ASSERT_DELTA(d.camera_pose().axis_angle().x(), m.camera_pose().axis_angle().x(), 1e-14);
      TS_ASSERT_DELTA(d.camera_pose().axis_angle().y(), m.camera_pose().axis_angle().y(), 1e-14);
      TS_ASSERT_DELTA(d.camera_pose().axis_angle().z(), m.camera_pose().axis_angle().z(), 1e-14);
      TS_ASSERT_EQUALS(d.camera_center().x(), m.camera_center().x());
      TS_ASSERT_EQUALS(d.camera_center().y(), m.camera_center().y());
      TS_ASSERT_EQUALS(d.camera_center().z(), m.camera_center().z());   
      Vector4 tsai_d( dynamic_cast<TsaiLensDistortion&>(*(d.lens_distortion())).distortion_parameters() );
      TS_ASSERT_EQUALS(tsai_d(0), tsaiv(0));
      TS_ASSERT_EQUALS(tsai_d(1), tsaiv(1));
      TS_ASSERT_EQUALS(tsai_d(2), tsaiv(2));
      TS_ASSERT_EQUALS(tsai_d(3), tsaiv(3));
    }
  }

  void test_pinholemodel_serialize_three() {
    std::srand(std::time(0));
    
    for (int i = 0; i < 10; i++) {
      Vector3 cc(rand(), rand(), rand());
      Vector3 rv(rand(), rand(), rand());
      Matrix3x3 rm(vw::math::axis_angle_to_matrix(rv));
      double fu = rand(), fv = rand(), cu = rand(), cv = rand();
      Vector3 u(rand(), rand(), rand());
      Vector3 v(rand(), rand(), rand());
      Vector3 w(rand(), rand(), rand());
      Vector4 tsaiv(vw::Vector4(rand(), rand(), rand(), rand()));
      TsaiLensDistortion tsai(tsaiv);
      PinholeModel m(cc, rm, fu, fv, cu, cv, u, v, w, tsai);

      vw::Vector<double> serial(serialize_pinholemodel<PinholeModelSerializeIntrinsic, PinholeModelSerializeRotation, PinholeModelSerializeTranslation>(m));

      // check if correctly serialized
      TS_ASSERT_EQUALS(serial(0), fu);
      TS_ASSERT_EQUALS(serial(1), fv);
      TS_ASSERT_EQUALS(serial(2), cu);
      TS_ASSERT_EQUALS(serial(3), cv);
      Vector3 model_rv(m.camera_pose().axis_angle());
      TS_ASSERT_EQUALS(serial(4), model_rv.x());
      TS_ASSERT_EQUALS(serial(5), model_rv.y());
      TS_ASSERT_EQUALS(serial(6), model_rv.z());
      Vector3 model_cc(m.camera_center());
      TS_ASSERT_EQUALS(serial(7), model_cc.x());
      TS_ASSERT_EQUALS(serial(8), model_cc.y());
      TS_ASSERT_EQUALS(serial(9), model_cc.z());

      // test deserialization
      PinholeModel d;
      d.set_coordinate_frame(u, v, w);
      deserialize_pinholemodel<PinholeModelSerializeIntrinsic, PinholeModelSerializeRotation, PinholeModelSerializeTranslation>(d, serial);
      double d_fu, d_fv, d_cu, d_cv;
      d.intrinsic_parameters(d_fu, d_fv, d_cu, d_cv);
      TS_ASSERT_EQUALS(fu, d_fu);
      TS_ASSERT_EQUALS(fv, d_fv);
      TS_ASSERT_EQUALS(cu, d_cu);
      TS_ASSERT_EQUALS(cv, d_cv);
      TS_ASSERT_DELTA(d.camera_pose().axis_angle().x(), m.camera_pose().axis_angle().x(), 1e-14);
      TS_ASSERT_DELTA(d.camera_pose().axis_angle().y(), m.camera_pose().axis_angle().y(), 1e-14);
      TS_ASSERT_DELTA(d.camera_pose().axis_angle().z(), m.camera_pose().axis_angle().z(), 1e-14);
      TS_ASSERT_EQUALS(d.camera_center().x(), m.camera_center().x());
      TS_ASSERT_EQUALS(d.camera_center().y(), m.camera_center().y());
      TS_ASSERT_EQUALS(d.camera_center().z(), m.camera_center().z());    
      }
    }

  void test_pinholemodel_serialize_two() {
    std::srand(std::time(0));
    
    for (int i = 0; i < 10; i++) {
      Vector3 cc(rand(), rand(), rand());
      Vector3 rv(rand(), rand(), rand());
      Matrix3x3 rm(vw::math::axis_angle_to_matrix(rv));
      double fu = rand(), fv = rand(), cu = rand(), cv = rand();
      Vector3 u(rand(), rand(), rand());
      Vector3 v(rand(), rand(), rand());
      Vector3 w(rand(), rand(), rand());
      Vector4 tsaiv(vw::Vector4(rand(), rand(), rand(), rand()));
      TsaiLensDistortion tsai(tsaiv);
      PinholeModel m(cc, rm, fu, fv, cu, cv, u, v, w, tsai);

      vw::Vector<double> serial(serialize_pinholemodel<PinholeModelSerializeIntrinsic, PinholeModelSerializeRotation>(m));

      // check if correctly serialized
      TS_ASSERT_EQUALS(serial(0), fu);
      TS_ASSERT_EQUALS(serial(1), fv);
      TS_ASSERT_EQUALS(serial(2), cu);
      TS_ASSERT_EQUALS(serial(3), cv);
      Vector3 model_rv(m.camera_pose().axis_angle());
      TS_ASSERT_EQUALS(serial(4), model_rv.x());
      TS_ASSERT_EQUALS(serial(5), model_rv.y());
      TS_ASSERT_EQUALS(serial(6), model_rv.z());

      // test deserialization
      PinholeModel d;
      d.set_coordinate_frame(u, v, w);
      deserialize_pinholemodel<PinholeModelSerializeIntrinsic, PinholeModelSerializeRotation>(d, serial);
      double d_fu, d_fv, d_cu, d_cv;
      d.intrinsic_parameters(d_fu, d_fv, d_cu, d_cv);
      TS_ASSERT_EQUALS(fu, d_fu);
      TS_ASSERT_EQUALS(fv, d_fv);
      TS_ASSERT_EQUALS(cu, d_cu);
      TS_ASSERT_EQUALS(cv, d_cv);
      TS_ASSERT_DELTA(d.camera_pose().axis_angle().x(), m.camera_pose().axis_angle().x(), 1e-14);
      TS_ASSERT_DELTA(d.camera_pose().axis_angle().y(), m.camera_pose().axis_angle().y(), 1e-14);
      TS_ASSERT_DELTA(d.camera_pose().axis_angle().z(), m.camera_pose().axis_angle().z(), 1e-14);
      }
    }

  void test_pinholemodel_serialize_one() {
    std::srand(std::time(0));
    
    for (int i = 0; i < 10; i++) {
      Vector3 cc(rand(), rand(), rand());
      Vector3 rv(rand(), rand(), rand());
      Matrix3x3 rm(vw::math::axis_angle_to_matrix(rv));
      double fu = rand(), fv = rand(), cu = rand(), cv = rand();
      Vector3 u(rand(), rand(), rand());
      Vector3 v(rand(), rand(), rand());
      Vector3 w(rand(), rand(), rand());
      Vector4 tsaiv(vw::Vector4(rand(), rand(), rand(), rand()));
      TsaiLensDistortion tsai(tsaiv);
      PinholeModel m(cc, rm, fu, fv, cu, cv, u, v, w, tsai);

      vw::Vector<double> serial(serialize_pinholemodel<PinholeModelSerializeIntrinsic>(m));

      // check if correctly serialized
      TS_ASSERT_EQUALS(serial(0), fu);
      TS_ASSERT_EQUALS(serial(1), fv);
      TS_ASSERT_EQUALS(serial(2), cu);
      TS_ASSERT_EQUALS(serial(3), cv);

      // test deserialization
      PinholeModel d;
      d.set_coordinate_frame(u, v, w);
      deserialize_pinholemodel<PinholeModelSerializeIntrinsic>(d, serial);
      double d_fu, d_fv, d_cu, d_cv;
      d.intrinsic_parameters(d_fu, d_fv, d_cu, d_cv);
      TS_ASSERT_EQUALS(fu, d_fu);
      TS_ASSERT_EQUALS(fv, d_fv);
      TS_ASSERT_EQUALS(cu, d_cu);
      TS_ASSERT_EQUALS(cv, d_cv);
    }
  }

  void test_pinholemodel_calibrate() {
    repeatable_srand(); // initialize reproducible srand; local implementation
    
    Vector3 cc(1, 1, 1);
    Vector3 rv(1, 1, 1);
    Matrix3x3 rm(vw::math::axis_angle_to_matrix(rv));
    double fu = 1000, fv = 500, cu = 200, cv = 100;
    Vector3 u(1, 0, 0);
    Vector3 v(0, 1, 0);
    Vector3 w(0, 0, 1);
    Vector4 tsaiv(1e-2, 1e-2, 1e-2, 1e-2);
    TsaiLensDistortion tsai(tsaiv);
    PinholeModel m(cc, rm, fu, fv, cu, cv, u, v, w, tsai);

    for (int ni = 0; ni < 3; ni++) {
      std::vector<vw::Vector3> points; // in 3 space
      std::vector<vw::Vector2> pixels; // from projection through m 
      int n = 30;
      for (int i = 0; i < n; i++) {
	Vector3 p(repeatable_rand(), repeatable_rand(), repeatable_rand());
	points.push_back(p);
	vw::Vector2 noise(repeatable_rand(10), repeatable_rand(10)); // add some noise!
	pixels.push_back(m.point_to_pixel(p) + noise); 
      }
    
      double mean = mean_error(m, pixels, points);

      // see if the optimizer improves the results; the mean error should decrease as 
      // the number of variables the optimizer gets to play with increases    
      // These tests can fail on real data depending on your srand implementation without anything being actually broken
      // - which is why a "fake" srand/rand function pair local to this class was implemented to generate 
      // pseudo pseudo random values (i.e. sequence that is always the same) for point positions and noise - 
      // thus the test should be completely repeatable
      { 
	PinholeModel c(m);
	pinholemodel_calibrate<PinholeModelSerializeIntrinsic>(c, pixels, points, 1000);
	double new_mean = mean_error(c, pixels, points);
	TS_ASSERT_LESS_THAN(new_mean, mean);
      }

      { 
	PinholeModel c(m);
	pinholemodel_calibrate<PinholeModelSerializeTSAI>(c, pixels, points, 1000);
	double new_mean = mean_error(c, pixels, points);
	TS_ASSERT_LESS_THAN(new_mean, mean);
      }

      { 
	PinholeModel c(m);
	pinholemodel_calibrate<PinholeModelSerializeRotation>(c, pixels, points, 1000);
	double new_mean = mean_error(c, pixels, points);
	TS_ASSERT_LESS_THAN(new_mean, mean);
      }

      { 
	PinholeModel c(m);
	pinholemodel_calibrate<PinholeModelSerializeTranslation>(c, pixels, points, 1000);
	double new_mean = mean_error(c, pixels, points);
	TS_ASSERT_LESS_THAN(new_mean, mean);
	mean = new_mean;
      }

      { 
	PinholeModel c(m);
	pinholemodel_calibrate<PinholeModelSerializeTranslation, PinholeModelSerializeRotation>(c, pixels, points, 1000);
	double new_mean = mean_error(c, pixels, points);
	TS_ASSERT_LESS_THAN(new_mean, mean);
	mean = new_mean;
      }

      { 
	PinholeModel c(m);
	pinholemodel_calibrate<PinholeModelSerializeTranslation, PinholeModelSerializeRotation, PinholeModelSerializeIntrinsic>(c, pixels, points, 1000);
	double new_mean = mean_error(c, pixels, points);
	TS_ASSERT_LESS_THAN(new_mean, mean);
	mean = new_mean;
      }

      { 
	PinholeModel c(m);
	pinholemodel_calibrate<PinholeModelSerializeTranslation, PinholeModelSerializeRotation, PinholeModelSerializeIntrinsic, PinholeModelSerializeTSAI>(c, pixels, points, 1000);
	double new_mean = mean_error(c, pixels, points);
	TS_ASSERT_LESS_THAN(new_mean, mean);
	mean = new_mean;
      }
    }
  }

  void test_pinholemodel_calibrate_ransac() {
    repeatable_srand(); // initialize reproducible srand; local implementation
    
    Vector3 cc(1, 1, 1);
    Vector3 rv(1, 1, 1);
    Matrix3x3 rm(vw::math::axis_angle_to_matrix(rv));
    double fu = 1000, fv = 500, cu = 200, cv = 100;
    Vector3 u(1, 0, 0);
    Vector3 v(0, 1, 0);
    Vector3 w(0, 0, 1);
    Vector4 tsaiv(1e-2, 1e-2, 1e-2, 1e-2);
    TsaiLensDistortion tsai(tsaiv);
    PinholeModel m(cc, rm, fu, fv, cu, cv, u, v, w, tsai);

    std::vector<vw::Vector3> points; // in 3 space
    std::vector<vw::Vector2> pixels; // from projection through m 
    int n = 30;
    for (int i = 0; i < n; i++) {
      Vector3 p(repeatable_rand(), repeatable_rand(), repeatable_rand());
      points.push_back(p);
      vw::Vector2 noise(repeatable_rand(20), repeatable_rand(20)); // add some noise!
      pixels.push_back(m.point_to_pixel(p) + noise); 
    }
    
    double mean = mean_error(m, pixels, points);

    // these tests verify mainly that RANSAC actually does respect the inlier_threshold value passed in
    // (The maximum error in the image plane for a resulting camera model can be equal to inlier_threshold)
    // and that the number of inliers is not "too small" 

    // play with these parameters to make the test more or less stringent and/or fast
    double inlier_threshold = 15;
    const int ransac_inlier_threshold = 10; // how many inliers to we require
    const int ransac_iter = 20; // number of ransac iterations
    const int lm_iter = 5; // number of levenberg marquardt iterations at every ransac iteration
    { 
      PinholeModel c(m);
      std::vector<int> inliers( pinholemodel_calibrate_ransac<PinholeModelSerializeIntrinsic>(c, pixels, points, inlier_threshold, ransac_iter, lm_iter) );
      TS_ASSERT_LESS_THAN(mean_error(c, pixels, points, inliers), inlier_threshold);
      TS_ASSERT_LESS_THAN(ransac_inlier_threshold, inliers.size()); // only critical if this fails repeatably, as ransac classification depends on random numbers
    }

    { 
      PinholeModel c(m);
      std::vector<int> inliers( pinholemodel_calibrate_ransac<PinholeModelSerializeTSAI>(c, pixels, points, inlier_threshold, ransac_iter, lm_iter) );
      TS_ASSERT_LESS_THAN(mean_error(c, pixels, points, inliers), inlier_threshold);
      TS_ASSERT_LESS_THAN(ransac_inlier_threshold, inliers.size()); // only critical if this fails repeatably, as ransac classification depends on random numbers
    }

    { 
      PinholeModel c(m);
      std::vector<int> inliers( pinholemodel_calibrate_ransac<PinholeModelSerializeRotation>(c, pixels, points, inlier_threshold, ransac_iter, lm_iter) );
      TS_ASSERT_LESS_THAN(mean_error(c, pixels, points, inliers), inlier_threshold);
      TS_ASSERT_LESS_THAN(ransac_inlier_threshold, inliers.size()); // only critical if this fails repeatably, as ransac classification depends on random numbers
    }

    { 
      PinholeModel c(m);
      std::vector<int> inliers( pinholemodel_calibrate_ransac<PinholeModelSerializeTranslation>(c, pixels, points, inlier_threshold, ransac_iter, lm_iter) );
      TS_ASSERT_LESS_THAN(mean_error(c, pixels, points, inliers), inlier_threshold);
      TS_ASSERT_LESS_THAN(ransac_inlier_threshold, inliers.size()); // only critical if this fails repeatably, as ransac classification depends on random numbers
    }

    { 
      PinholeModel c(m);
      std::vector<int> inliers( pinholemodel_calibrate_ransac<PinholeModelSerializeTranslation, PinholeModelSerializeRotation>(c, pixels, points, inlier_threshold, ransac_iter, lm_iter) );
      TS_ASSERT_LESS_THAN(mean_error(c, pixels, points, inliers), inlier_threshold);
      TS_ASSERT_LESS_THAN(ransac_inlier_threshold, inliers.size()); // only critical if this fails repeatably, as ransac classification depends on random numbers
    }

    { 
      PinholeModel c(m);
      std::vector<int> inliers( pinholemodel_calibrate_ransac<PinholeModelSerializeTranslation, PinholeModelSerializeRotation, PinholeModelSerializeIntrinsic>(c, pixels, points, inlier_threshold, ransac_iter, lm_iter) );
      TS_ASSERT_LESS_THAN(mean_error(c, pixels, points, inliers), inlier_threshold);
      TS_ASSERT_LESS_THAN(ransac_inlier_threshold, inliers.size()); // only critical if this fails repeatably, as ransac classification depends on random numbers
    }

    { 
      PinholeModel c(m);
      std::vector<int> inliers( pinholemodel_calibrate_ransac<PinholeModelSerializeTranslation, PinholeModelSerializeRotation, PinholeModelSerializeIntrinsic, PinholeModelSerializeTSAI>(c, pixels, points, inlier_threshold, ransac_iter, lm_iter) );
      TS_ASSERT_LESS_THAN(mean_error(c, pixels, points, inliers), inlier_threshold);
      TS_ASSERT_LESS_THAN(ransac_inlier_threshold, inliers.size()); // only critical if this fails repeatably, as ransac classification depends on random numbers
    }
  }
};
