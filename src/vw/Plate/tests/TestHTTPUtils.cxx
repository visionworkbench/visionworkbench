#include <gtest/gtest.h>
#include <vw/Plate/HTTPUtils.h>

using namespace std;
using namespace vw;
using namespace vw::platefile;

TEST(HTTPUtils, UrlUnescape) {
  EXPECT_EQ("rawr pants", url_unescape("rawr%20pants"));
  EXPECT_EQ("rawr pants", url_unescape("rawr+pants"));
  EXPECT_EQ("?",  url_unescape("%3f"));
  EXPECT_EQ("?",  url_unescape("%3F"));
  EXPECT_EQ("/woo/",      url_unescape("%2fwoo%2f"));
}

TEST(HTTPUtils, UrlEscape) {
  EXPECT_EQ("%2fpants%2fcheese%2f", url_escape("/pants/cheese/", ""));
  EXPECT_EQ("%2fpants%2fcheese%2f", url_escape("/pants/cheese/"));
  EXPECT_EQ("/pants/cheese/", url_escape("/pants/cheese/", "/"));
  EXPECT_EQ("wee+moo", url_escape("wee moo"));
}

TEST(HTTPUtils, Empty) {
  EXPECT_EQ("", url_escape(""));
  EXPECT_EQ("", url_unescape(""));
}

TEST(HTTPUtils, QueryMap) {
  QueryMap a;
  EXPECT_EQ("",     a.serialize());
  EXPECT_EQ("rawr", a.get("key", string("rawr")));

  QueryMap b("pants=cheese&rawr=moo");
  EXPECT_EQ("?pants=cheese&rawr=moo", b.serialize());
  EXPECT_EQ("moo", b.get("rawr", string("pants")));
  EXPECT_EQ("pants", b.get("rawr2", string("pants")));
}
