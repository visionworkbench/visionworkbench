#include <vw/Plate/Datastore.h>
#include <vw/Plate/HTTPUtils.h>
#include <vw/Plate/detail/Blobstore.h>

namespace vw { namespace platefile {

Datastore* Datastore::open(const Url& url) {
  return new detail::Blobstore(url);
  //vw_throw(NoImplErr() << "Unsupported datastore scheme: " << url.scheme());
}

Datastore* Datastore::open(const Url& url, const IndexHeader& d) {
  return new detail::Blobstore(url, d);
  //vw_throw(NoImplErr() << "Unsupported datastore scheme: " << url.scheme());
}

Datastore::TileSearch& Datastore::get(TileSearch& buf, uint32 level, uint32 row, uint32 col, TransactionRange range, uint32 limit) {
  this->head(buf, level, row, col, range, limit);
  return this->populate(buf);
}

Datastore::TileSearch& Datastore::get(TileSearch& buf, uint32 level,   const BBox2u& region, TransactionRange range, uint32 limit) {
  this->head(buf, level, region, range, limit);
  return this->populate(buf);
}

uint32 Datastore::id() const                    { return index_header().platefile_id(); }
uint32 Datastore::num_levels() const            { return index_header().num_levels(); }
uint32 Datastore::tile_size() const             { return index_header().tile_size(); }
std::string Datastore::tile_filetype() const    { return index_header().tile_filetype(); }
PixelFormatEnum Datastore::pixel_format() const { return static_cast<PixelFormatEnum>(index_header().pixel_format()); }
ChannelTypeEnum Datastore::channel_type() const { return static_cast<ChannelTypeEnum>(index_header().channel_type()); }

}}
