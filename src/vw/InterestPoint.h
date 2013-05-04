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


/// \file InterestPoint.h
///
/// A convenience header that includes the header files in vw/InterestPoint.
///
#ifndef __VW_INTERESTPOINT_H__
#define __VW_INTERESTPOINT_H__

// Data Types
#include <vw/InterestPoint/InterestData.h>
#include <vw/InterestPoint/InterestTraits.h>
#include <vw/InterestPoint/ImageOctave.h>
#include <vw/InterestPoint/ImageOctaveHistory.h>
#include <vw/InterestPoint/IntegralImage.h>

// Operators
#include <vw/InterestPoint/InterestOperator.h>
#include <vw/InterestPoint/BoxFilter.h>
#include <vw/InterestPoint/IntegralInterestOperator.h>

// Detection
#include <vw/InterestPoint/Extrema.h>
#include <vw/InterestPoint/Localize.h>
#include <vw/InterestPoint/Detector.h>
#include <vw/InterestPoint/IntegralDetector.h>

// Description
#include <vw/InterestPoint/LearnPCA.h>
#include <vw/InterestPoint/WeightedHistogram.h>
#include <vw/InterestPoint/Descriptor.h>
#include <vw/InterestPoint/IntegralDescriptor.h>

// Matching
#include <vw/InterestPoint/Matcher.h>

#endif // __VW_INTERESTPOINT_H__

