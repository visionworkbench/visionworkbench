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


#ifndef __VW_PLATE_TOAST_PLATEMANAGER_H__
#define __VW_PLATE_TOAST_PLATEMANAGER_H__

#include <vw/Plate/PlateManager.h>

namespace vw {
namespace platefile {

  template <class PixelT>
  class ToastPlateManager : public PlateManager<PixelT> {
  protected:

    // Transforms an input image with input georef to toast projection
    // at the most ideal matching pyramid level.
    virtual void transform_image( cartography::GeoReference const& georef,
                                  ImageViewRef<PixelT>& image,
                                  TransformRef& txref, int& level ) const;

 public:
    ToastPlateManager(boost::shared_ptr<PlateFile> platefile) :
      PlateManager<PixelT>(platefile) {}

    // Non-working .. just matching base type
    cartography::GeoReference georeference( int level ) const;

    void generate_mipmap_tile(ImageView<PixelT>& dest,
        int col, int row, int level, TransactionOrNeg transaction_id, bool preblur) const;

    ImageView<PixelT> fetch_child_tile(
        int x, int y, int level, TransactionOrNeg transaction_id) const;

    virtual void mipmap(uint32 starting_level, BBox2i const& bbox, TransactionOrNeg input_transaction_id, bool preblur, const ProgressCallback &progress_callback = ProgressCallback::dummy_instance()) const;
  };

}} // namespace vw::plate

#endif // __VW_PLATE_TOAST_PLATEMANAGER_H__
