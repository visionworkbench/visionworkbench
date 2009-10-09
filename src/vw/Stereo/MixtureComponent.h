#pragma once

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
  
      
