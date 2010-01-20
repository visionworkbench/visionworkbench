// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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
#include <vw/Math/FundamentalMatrix.h>
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
#include <vw/Math/PoseEstimation.h>
#include <vw/Math/RANSAC.h>

#if defined(VW_HAVE_PKG_LAPACK) && VW_HAVE_PKG_LAPACK==1
#include <vw/Math/LinearAlgebra.h>
#include <vw/Math/LevenbergMarquardt.h>
#include <vw/Math/Geometry.h>
#endif

#endif // __VW_MATH_H__
