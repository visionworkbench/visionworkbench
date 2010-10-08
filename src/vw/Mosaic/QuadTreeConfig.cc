#include <vw/Mosaic/QuadTreeConfig.h>
#include <vw/Mosaic/CelestiaQuadTreeConfig.h>
#include <vw/Mosaic/GigapanQuadTreeConfig.h>
#include <vw/Mosaic/GMapQuadTreeConfig.h>
#include <vw/Mosaic/KMLQuadTreeConfig.h>
#include <vw/Mosaic/TMSQuadTreeConfig.h>
#include <vw/Mosaic/ToastQuadTreeConfig.h>
#include <vw/Mosaic/UniviewQuadTreeConfig.h>

using namespace vw;
using namespace vw::mosaic;

boost::shared_ptr<QuadTreeConfig> QuadTreeConfig::make(const std::string& type) {
  typedef boost::shared_ptr<QuadTreeConfig> ptr_t;

  if (type == "celestia")
    return ptr_t(new CelestiaQuadTreeConfig());
  else if (type == "gigapan" || type == "gigapan-noproj")
    return ptr_t(new GigapanQuadTreeConfig());
  else if (type == "gmap")
    return ptr_t(new GMapQuadTreeConfig());
  else if (type == "kml")
    return ptr_t(new KMLQuadTreeConfig());
  else if (type == "tms")
    return ptr_t(new TMSQuadTreeConfig());
  else if (type == "uniview")
    return ptr_t(new UniviewQuadTreeConfig());
  else
    vw_throw(NoImplErr() << "Unknown quad tree type: " << type);
}
