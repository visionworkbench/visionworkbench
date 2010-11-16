#include <vw/tools/Common.h>
#include <vw/FileIO/DiskImageResource.h>
#include <boost/scoped_ptr.hpp>

vw::ImageFormat vw::tools::taste_image(const std::string& filename) {
  boost::scoped_ptr<vw::SrcImageResource> src(vw::DiskImageResource::open(filename));
  return src->format();
}
