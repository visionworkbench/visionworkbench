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


/// \file Math.h
///
/// A convenience header that includes the header files in vw/Math.
///
#ifndef __VW_MATH_H__
#define __VW_MATH_H__

#include <vw/config.h>

#include <vw/Math/Functions.h>
#include <vw/Math/Functors.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/BBox.h>
#include <vw/Math/Quaternion.h>
#include <vw/Math/ConjugateGradient.h>
#include <vw/Math/NelderMead.h>
#include <vw/Math/ParticleSwarmOptimization.h>
#include <vw/Math/EulerAngles.h>
#include <vw/Math/Statistics.h>
#include <vw/Math/DisjointSet.h>
#include <vw/Math/MinimumSpanningTree.h>
#include <vw/Math/RANSAC.h>

#if defined(VW_HAVE_PKG_LAPACK) && VW_HAVE_PKG_LAPACK==1
#include <vw/Math/LinearAlgebra.h>
#include <vw/Math/LevenbergMarquardt.h>
#include <vw/Math/Geometry.h>
#include <vw/Math/PoseEstimation.h>
#endif

#endif // __VW_MATH_H__
