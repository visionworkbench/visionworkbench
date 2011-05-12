// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_STEREO_MIXTURE_COMPONENT_H__
#define __VW_STEREO_MIXTURE_COMPONENT_H__

#include <vw/Image.h>

template<class ImplT>
struct MixtureComponentBase : private boost::noncopyable {
    inline ImplT& impl() { return static_cast<ImplT&>(*this); }
    inline ImplT const& impl() const { return static_cast<ImplT const&>(*this); }

  protected:
    MixtureComponentBase() { }
};

#endif//__VW_STEREO_MIXTURE_COMPONENT_H__
