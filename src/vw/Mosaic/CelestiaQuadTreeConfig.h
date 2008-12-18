// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
//
// Copyright 2006 Carnegie Mellon University. All rights reserved.
//
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
//
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
//
// __END_LICENSE__

/// \file CelestiaQuadTreeConfig.h
///
/// A configuration class that provides callbacks for
/// QuadTreeGenerator that generate Celestia "virtual textures".
///
/// Helpful doc on the subject: "Creating Textures for Celestia"
/// <http://www.lns.cornell.edu/~seb/celestia/textures.html>
///
#ifndef __VW_MOSAIC_CELESTIAQUADTREECONFIG_H__
#define __VW_MOSAIC_CELESTIAQUADTREECONFIG_H__

#include <vw/Mosaic/QuadTreeGenerator.h>

namespace vw {
namespace mosaic {

  // This class is overkill, but it exists by analogy to others
  // like it for consistency.
  class CelestiaQuadTreeConfig {
  public:
    void configure( QuadTreeGenerator& qtree ) const;

    // Makes paths of the form "path/name/level1/tx_2_1.jpg"
    // 2_1 is the tile at (x,y) location (2,1), (0,0) is upper-left
    static std::string image_path( QuadTreeGenerator const& qtree, std::string const& name );

  };

} // namespace mosaic
} // namespace vw

#endif // __VW_MOSAIC_CELESTIAQUADTREECONFIG_H__
