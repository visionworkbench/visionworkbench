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


#ifndef __VW_STEREO_UNIFORM_MIXTURE_COMPONENT_H__
#define __VW_STEREO_UNIFORM_MIXTURE_COMPONENT_H__

#include <vw/Stereo/MixtureComponent.h>

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

#endif//__VW_STEREO_UNIFORM_MIXTURE_COMPONENT__
