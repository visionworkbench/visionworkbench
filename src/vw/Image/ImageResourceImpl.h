// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// Header to avoid ifdefs in user code. Include this instead of the specific
// impls.
//
#ifndef __VW_IMAGE_IMAGERESOURCEIMPL_H__
#define __VW_IMAGE_IMAGERESOURCEIMPL_H__

#include <vw/config.h>

#define __INSIDE_VW_IMAGE_IMAGERESOURCEIMPL_H__ 1

#if defined(VW_HAVE_PKG_OPENCV) && VW_HAVE_PKG_OPENCV == 1
# include <vw/Image/ImageResourceOpenCV.h>
#endif

#include <vw/Image/ImageResourceStream.h>

#undef __INSIDE_VW_IMAGE_IMAGERESOURCEIMPL_H__

#endif
