// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_CARTOGRAPHY_BRESENHAMLINE_H__
#define __VW_CARTOGRAPHY_BRESENHAMLINE_H__

#include <vw/Core/FundamentalTypes.h>
#include <vw/Math/Vector.h>

namespace vw {
namespace cartography {

  class BresenhamLine {
    vw::int32 x0, y0, x1, y1;
    vw::int32 x, y;
    bool steep;
    vw::int32 deltax, deltay, error, ystep;

    void setup() {
      steep = abs(y1-y0) > abs(x1-x0);
      if (steep) {
        std::swap(x0,y0);
        std::swap(x1,y1);
      }
      if ( x0 > x1 ) {
        std::swap(x0,x1);
        std::swap(y0,y1);
      }
      deltax = x1 - x0;
      deltay = abs(y1-y0);
      error = deltax / 2;
      ystep = y0 < y1 ? 1 : -1;
      x = x0; y = y0;
    }

  public:
    BresenhamLine( vw::Vector2i const& start, vw::Vector2i const& stop ) :
      x0(start[0]), y0(start[1]), x1(stop[0]), y1(stop[1]) {
      setup();
    }

    BresenhamLine( vw::int32 mx0, vw::int32 my0, vw::int32 mx1, vw::int32 my1 ) :
      x0(mx0), y0(my0), x1(mx1), y1(my1) {
      setup();
    }

    vw::Vector2i operator*() const {
      if (steep)
        return vw::Vector2i(y,x);
      else
        return vw::Vector2i(x,y);
    }

    void operator++() {
      x++;
      error -= deltay;
      if ( error < 0 ) {
        y += ystep;
        error += deltax;
      }
    }

    bool is_good() const { return x < x1; }
  };

}}

#endif//__VW_CARTOGRAPHY_BRESENHAMLINE_H__
