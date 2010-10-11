// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Image.h
///
/// A convenience header that includes the header files in vw/Image.
///
#ifndef __VW_IMAGE_H__
#define __VW_IMAGE_H__

#include <vw/Image/PixelIterator.h>
#include <vw/Image/ImageViewBase.h>
#include <vw/Image/PixelAccessors.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageIO.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/PixelTypeInfo.h>
#include <vw/Image/PixelMath.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/MaskViews.h>
#include <vw/Image/PerPixelViews.h>
#include <vw/Image/PerPixelAccessorViews.h>
#include <vw/Image/UtilityViews.h>
#include <vw/Image/ImageMath.h>
#include <vw/Image/ImageResource.h>
#include <vw/Image/ImageResourceImpl.h>
#include <vw/Image/ImageResourceView.h>
#include <vw/Image/ViewImageResource.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/EdgeExtension.h>
#include <vw/Image/Interpolation.h>
#include <vw/Image/Convolution.h>
#include <vw/Image/Filter.h>
#include <vw/Image/Transform.h>
#include <vw/Image/BlockProcessor.h>
#include <vw/Image/BlockRasterize.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/Palette.h>
#include <vw/Image/SparseImageCheck.h>

#endif // __VW_IMAGE_H__
