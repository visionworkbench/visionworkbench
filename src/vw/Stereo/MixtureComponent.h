// __BEGIN_LICENSE__
// __END_LICENSE__

#ifndef __VW_STEREO_MIXTURE_COMPONENT_H__
#define __VW_STEREO_MIXTURE_COMPONENT_H__

#include <vw/Image.h>

template<class ImplT>
struct MixtureComponentBase {
  inline ImplT& impl() { return static_cast<ImplT&>(*this); }
  inline ImplT const& impl() const { return static_cast<ImplT const&>(*this); }

protected:
  MixtureComponentBase() { }
  MixtureComponentBase(MixtureComponentBase const&) { }
  MixtureComponentBase& operator=(MixtureComponentBase const&) { return *this; }
};

#endif//__VW_STEREO_MIXTURE_COMPONENT_H__
