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


// Header to avoid ifdefs in user code. Include this instead of the specific
// impls.
//
#ifndef __VW_IMAGE_IMAGERESOURCEIMPL_H__
#define __VW_IMAGE_IMAGERESOURCEIMPL_H__

#include <vw/vw_config.h>

#define __INSIDE_VW_IMAGE_IMAGERESOURCEIMPL_H__ 1

#if defined(VW_HAVE_PKG_OPENCV) && VW_HAVE_PKG_OPENCV == 1
# include <vw/Image/ImageResourceOpenCV.h>
#endif

#include <vw/Image/ImageResourceStream.h>

#undef __INSIDE_VW_IMAGE_IMAGERESOURCEIMPL_H__

#endif
