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

  Url u1("");
  EXPECT_EQ("file", u1.scheme());
  EXPECT_EQ(""    , u1.hostname());
  EXPECT_EQ(0     , u1.port());
  EXPECT_EQ(""    , u1.path());
  EXPECT_EQ(""    , u1.fragment());
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

TEST(HTTPUtils, Http) {
  Url u("http://rawr/pants?cheese=moo&pants=foo#wee=blarg");

  EXPECT_EQ("http",   u.scheme());
  EXPECT_EQ("rawr",   u.netloc());
  EXPECT_EQ("/pants", u.path());

  EXPECT_EQ("moo",  u.query().get("cheese", string("not provided")));
  EXPECT_EQ("foo",  u.query().get("pants",  string("not provided")));
  EXPECT_EQ("nope", u.query().get("rawr",   string("nope")));

  EXPECT_EQ("wee=blarg", u.fragment());
}

TEST(HTTPUtils, File) {
  Url a("/pants/cheese");
  EXPECT_EQ("file",          a.scheme());
  EXPECT_EQ("",              a.netloc());
  EXPECT_EQ("/pants/cheese", a.path());

  Url b("file:///pants/cheese");
  EXPECT_EQ("file",          a.scheme());
  EXPECT_EQ("",              a.netloc());
  EXPECT_EQ("/pants/cheese", a.path());

  EXPECT_EQ("file://file/", Url("file").url());
  EXPECT_EQ("file:///",     Url("file:").url());
  EXPECT_EQ("file:///",     Url("file:/").url());
  EXPECT_EQ("file:///",     Url("file://").url());
  //EXPECT_EQ("file:///moo",  Url("file:moo").url()); // XXX: This one is known-broken
  EXPECT_EQ("file:///moo",  Url("file:/moo").url());
  EXPECT_EQ("file:///moo",  Url("file:///moo").url());
  EXPECT_EQ("file://host/pants", Url("file://host/pants").url());
}

TEST(HTTPUtils, Tricky) {
  Url u("/pants:moo/");
  EXPECT_EQ("file", u.scheme());
  EXPECT_EQ("/pants:moo/", u.path());
}

TEST(HTTPUtils, PF) {
  Url u1("pf:///exchange/platefilename.plate");
  EXPECT_EQ("pf", u1.scheme());
  EXPECT_TRUE(u1.hostname().empty());
  EXPECT_EQ(0, u1.port());
  EXPECT_EQ("/exchange/platefilename.plate", u1.path());

  Url u3("pf://12.123.12.123:25/exchange/platefilename.plate");
  EXPECT_EQ("pf", u3.scheme());
  EXPECT_EQ("12.123.12.123", u3.hostname());
  EXPECT_EQ(25, u3.port());
  EXPECT_EQ("/exchange/platefilename.plate", u3.path());

  Url u4("pf://123.23.23.123/exchange/platefilename.plate");
  EXPECT_EQ("pf", u4.scheme());
  EXPECT_EQ("123.23.23.123", u4.hostname());
  EXPECT_EQ(0, u4.port());
  EXPECT_EQ("/exchange/platefilename.plate", u4.path());
}

TEST(HTTPUtils, Change) {
  Url u("pf://meow/pants/cheese?rawr=foo");
  u.netloc() = "woof";
  u.path() = u.path() + "/append";
  EXPECT_EQ("pf://woof/pants/cheese/append?rawr=foo", u.url());
}
