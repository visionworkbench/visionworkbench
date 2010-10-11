#include <vw/Mosaic/QuadTreeConfig.h>
#include <vw/Mosaic/CelestiaQuadTreeConfig.h>
#include <vw/Mosaic/GigapanQuadTreeConfig.h>
#include <vw/Mosaic/GMapQuadTreeConfig.h>
#include <vw/Mosaic/KMLQuadTreeConfig.h>
#include <vw/Mosaic/TMSQuadTreeConfig.h>
#include <vw/Mosaic/UniviewQuadTreeConfig.h>

using namespace vw;
using namespace vw::mosaic;

boost::shared_ptr<QuadTreeConfig> QuadTreeConfig::make(const std::string& type) {
  typedef boost::shared_ptr<QuadTreeConfig> ptr_t;
  std::string utype = boost::to_upper_copy(type);

  if (utype == "CELESTIA")
    return ptr_t(new CelestiaQuadTreeConfig());
  else if (utype == "GIGAPAN" || utype == "GIGAPAN_NOPROJ")
    return ptr_t(new GigapanQuadTreeConfig());
  else if (utype == "GMAP")
    return ptr_t(new GMapQuadTreeConfig());
  else if (utype == "KML")
    return ptr_t(new KMLQuadTreeConfig());
  else if (utype == "TMS")
    return ptr_t(new TMSQuadTreeConfig());
  else if (utype == "UNIVIEW")
    return ptr_t(new UniviewQuadTreeConfig());
  else
    vw_throw(NoImplErr() << "Unknown quad tree type: " << utype);
}
