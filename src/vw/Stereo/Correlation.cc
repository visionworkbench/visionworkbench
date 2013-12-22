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

#include <vw/Math/BBox.h>
#include <vw/Image/Statistics.h>
#include <vw/Stereo/Correlation.h>

namespace vw {
namespace stereo {

  bool
  SearchParamLessThan::operator()(SearchParam const& A,
                                  SearchParam const& B) const {
    // The amount of computation for correlation for a given
    // SearchParam instance is proportional to the product of the
    // dimensions of its search boxes.
    return search_volume(A) < search_volume(B);
  }

  inline int32 area( BBox2i const& a ) {
    int32 width = a.width();
    int32 heigh = a.height();
    if ( width < 0 || heigh < 0 )
      return 0;
    return width * heigh;
  }

  inline void expand_bbox( BBox2i& a, BBox2i const& b ) {

    // Treat the cases of empty boxes. Those come in many
    // shapes, some rather pathological.
    if (b.empty()) return;
    if (a.empty()) { a = b; return; }

    a.min().x() = std::min(a.min().x(),b.min().x());
    a.min().y() = std::min(a.min().y(),b.min().y());
    a.max().x() = std::max(a.max().x(),b.max().x());
    a.max().y() = std::max(a.max().y(),b.max().y());
  }

  bool subdivide_regions( ImageView<PixelMask<Vector2i> > const& disparity,
                          BBox2i const& current_bbox,
                          std::vector<SearchParam>& list,
                          Vector2i const& kernel_size,
                          int32 fail_count ) {

    // 1.) Is this region too small? Must we stop?
    if ( prod(current_bbox.size()) <= 200 ||
         current_bbox.width() < 16 || current_bbox.height() < 16 ){
      BBox2i expanded = current_bbox;
      expanded.expand(1);
      expanded.crop( bounding_box( disparity ) );
      PixelAccumulator<EWMinMaxAccumulator<Vector2i> > accumulator;
      for_each_pixel( crop(disparity, expanded), accumulator );
      if ( !accumulator.is_valid() ) return true;

      list.push_back( SearchParam( current_bbox,
                                   BBox2i(accumulator.minimum(),
                                          accumulator.maximum() + Vector2i(1,1) ) ) );
      return true;
    }

    // 2) Divide the section into 4 quadrants, does it reduce total search?
    Vector2i split_pt = current_bbox.size()/2;
    BBox2i q1( current_bbox.min(), current_bbox.min()+split_pt );
    BBox2i q4( current_bbox.min()+split_pt, current_bbox.max() );
    BBox2i q2( current_bbox.min() + Vector2i(split_pt[0],0),
               Vector2i(current_bbox.max()[0],current_bbox.min()[1]+split_pt[1]) );
    BBox2i q3( current_bbox.min() + Vector2i(0,split_pt[1]),
               Vector2i(current_bbox.min()[0]+split_pt[0],current_bbox.max()[1]) );
    BBox2i q1_search, q2_search, q3_search, q4_search;

    int32 split_search = 0;
    { // Q1
      PixelAccumulator<EWMinMaxAccumulator<Vector2i> > accumulator;
      for_each_pixel( crop(disparity,q1), accumulator );
      if ( accumulator.is_valid() ) {
        q1_search = BBox2i(accumulator.minimum(),
                           accumulator.maximum()+Vector2i(1,1));
        split_search += area(q1_search) * prod(q1.size()+kernel_size);
      }
    }
    { // Q2
      PixelAccumulator<EWMinMaxAccumulator<Vector2i> > accumulator;
      for_each_pixel( crop(disparity,q2), accumulator );
      if ( accumulator.is_valid() ) {
        q2_search = BBox2i(accumulator.minimum(),
                           accumulator.maximum()+Vector2i(1,1));
        split_search += area(q2_search) * prod(q2.size()+kernel_size);
      }
    }
    { // Q3
      PixelAccumulator<EWMinMaxAccumulator<Vector2i> > accumulator;
      for_each_pixel( crop(disparity,q3), accumulator );
      if ( accumulator.is_valid() ) {
        q3_search = BBox2i(accumulator.minimum(),
                           accumulator.maximum()+Vector2i(1,1));
        split_search += area(q3_search) * prod(q3.size()+kernel_size);
      }
    }
    { // Q4
      PixelAccumulator<EWMinMaxAccumulator<Vector2i> > accumulator;
      for_each_pixel( crop(disparity,q4), accumulator );
      if ( accumulator.is_valid() ) {
        q4_search = BBox2i(accumulator.minimum(),
                           accumulator.maximum()+Vector2i(1,1));
        split_search += area(q4_search) * prod(q4.size()+kernel_size);
      }
    }

    // 3) Find current search v2
    BBox2i current_search_region;
    if ( q1_search != BBox2i() )
      current_search_region = q1_search;
    if ( q2_search != BBox2i() && current_search_region == BBox2i() )
      current_search_region = q2_search;
    else
      expand_bbox( current_search_region, q2_search );
    if ( q3_search != BBox2i() && current_search_region == BBox2i() )
      current_search_region = q3_search;
    else
      expand_bbox( current_search_region, q3_search );
    if ( q4_search != BBox2i() && current_search_region == BBox2i() )
      current_search_region = q4_search;
    else
      expand_bbox( current_search_region, q4_search );
    int32 current_search = area(current_search_region) * prod(current_bbox.size()+kernel_size);

    if ( split_search > current_search*0.9 && fail_count == 0 ) {
      // Did bad .. maybe next level will have better luck?
      std::vector<SearchParam> failed;
      if (!subdivide_regions( disparity, q1, list, kernel_size, fail_count + 1 ) )
        failed.push_back(SearchParam(q1,q1_search));
      if (!subdivide_regions( disparity, q2, list, kernel_size, fail_count + 1 ) )
        failed.push_back(SearchParam(q2,q2_search));
      if (!subdivide_regions( disparity, q3, list, kernel_size, fail_count + 1 ) )
        failed.push_back(SearchParam(q3,q3_search));
      if (!subdivide_regions( disparity, q4, list, kernel_size, fail_count + 1 ) )
        failed.push_back(SearchParam(q4,q4_search));
      if ( failed.size() == 4 ) {
        // all failed, push back this region as a whole
        list.push_back( SearchParam( current_bbox,
                                     current_search_region ) );
        return true;
      } else if ( failed.size() == 3 ) {
        // 3 failed to split can I merge ?
        std::vector<SearchParam>::const_iterator it1 = failed.begin(), it2 = failed.begin();
        ++it2;
        if ( ( it1->first.min().x() == it2->first.min().x() ||
               it1->first.min().y() == it2->first.min().y() ) &&
             it1->second == it2->second ) {
          BBox2i merge = it1->first;
          expand_bbox( merge, it2->first );
          list.push_back( SearchParam( merge, it1->second ) );
          list.push_back( *++it2 );
          return true;
        }
        ++it1; ++it2;
        if ( ( it1->first.min().x() == it2->first.min().x() ||
               it1->first.min().y() == it2->first.min().y() ) &&
             it1->second == it2->second ) {
          BBox2i merge = it1->first;
          expand_bbox( merge, it2->first );
          list.push_back( SearchParam( merge, it1->second ) );
          list.push_back( failed.front() );
          return true;
        }
        it1 = failed.begin();
        if ( ( it1->first.min().x() == it2->first.min().x() ||
               it1->first.min().y() == it2->first.min().y() ) &&
             it1->second == it2->second ) {
          BBox2i merge = it1->first;
          expand_bbox( merge, it2->first );
          list.push_back( SearchParam( merge, it1->second ) );
          list.push_back( *++it1 );
          return true;
        }
        // Push only the bombed regions, possibly a merge step could go here
        list.insert( list.end(), failed.begin(), failed.end() );
      } else if ( failed.size() == 2 ) {
        // 2 failed to split .. can I merge?
        if ( ( failed.front().first.min().x() == failed.back().first.min().x() ||
               failed.front().first.min().y() == failed.back().first.min().y() ) &&
             failed.front().second == failed.back().second ) {
          BBox2i merge = failed.front().first;
          expand_bbox( merge, failed.back().first );
          list.push_back( SearchParam( merge, failed.front().second ) );
          return true;
        }
        list.insert( list.end(), failed.begin(), failed.end() );
      } else if ( failed.size() == 1 ) {
        // Only 1 failed to split .. push it back
        list.push_back( failed.front() );
      }
      return true;
    } else if ( split_search > current_search*0.9 && fail_count > 0 ) {
      // Bad again .. back up
      return false;
    } else {
      // Good split
      subdivide_regions( disparity, q1, list, kernel_size );
      subdivide_regions( disparity, q2, list, kernel_size );
      subdivide_regions( disparity, q3, list, kernel_size );
      subdivide_regions( disparity, q4, list, kernel_size );
    }
    return true;
  }

}} // end namespace vw::stereo
