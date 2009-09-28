#pragma once
#include "MixtureComponent.h"

#include <vw/Image.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Math.h>

namespace vw {
  template<class PrecisionT>
    class UniformMixtureComponent : public MixtureComponentBase<UniformMixtureComponent<PrecisionT> > {
  public:
  UniformMixtureComponent(Vector2i window_size) : w(window_size(0), window_size(1)), p(window_size(0), window_size(1))
    {
      fill(w, 1.);
      fill(p, 1.);
    }
    
    void reset(BBox2i left_window) {
      fill(w, 1.);
      fill(p, 1.);
    }
    
    void fit_parameters() { }

    void update_posterior() { }
    
    PrecisionT log_likelihood() { return 0; }

    ImageView<PrecisionT> & weights() { return w; }
    
    ImageView<PrecisionT> const& weights() const { return w; }
    
    ImageView<PrecisionT> const& prob() const { return p; }
    
  private:
    ImageView<PrecisionT> w;  // weights
    ImageView<PrecisionT> p; // posterior probabilities of data in window    
  };

}
