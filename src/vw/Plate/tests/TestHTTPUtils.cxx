// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


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

  Url u1;
  EXPECT_EQ("", u1.scheme());
  EXPECT_EQ("", u1.hostname());
  EXPECT_EQ(0 , u1.port());
  EXPECT_EQ("", u1.netloc());
  EXPECT_EQ("/", u1.path());
  EXPECT_EQ("", u1.fragment());
}

TEST(HTTPUtils, QueryMap) {
  QueryMap a;
  EXPECT_EQ("",     a.serialize());
  EXPECT_EQ("rawr", a.get("key", string("rawr")));

  QueryMap b("pants=cheese&rawr=moo");
  // This technically isn't guaranteed... there's no map ordering.
  EXPECT_EQ("?pants=cheese&rawr=moo", b.serialize());
  EXPECT_EQ("moo", b.get("rawr", string("pants")));
  EXPECT_EQ("pants", b.get("rawr2", string("pants")));

  EXPECT_TRUE(b.has("rawr"));
  EXPECT_NO_THROW(b.get<string>("rawr"));
  EXPECT_THROW(b.get<int>("rawr2"), LogicErr);
  EXPECT_THROW(b.get<int>("rawr"), ArgumentErr);

  EXPECT_EQ("none", b.get("moon", string("none")));
  b.set("moon", "stars");
  EXPECT_EQ("stars", b.get("moon", string("none")));
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

  EXPECT_EQ("file://file",       Url("file").string());

  EXPECT_EQ("file://file%3a",    Url("file:").string());
  EXPECT_EQ("file:",             Url("file:").path());

  EXPECT_EQ("file://file%3a/",   Url("file:/").string());
  EXPECT_EQ("file:/",            Url("file:/").path());

  EXPECT_EQ("file:///",          Url("file://").string());

  EXPECT_EQ("file://file%3a/moo", Url("file:/moo").string());
  EXPECT_EQ("file:/moo",          Url("file:/moo").path());

  EXPECT_EQ("file:///moo",       Url("file:///moo").string());
  EXPECT_EQ("file://host/pants", Url("file://host/pants").string());

  Url c("rawr");
  Url d("file://rawr");

  EXPECT_EQ("rawr",        c.path());
  EXPECT_EQ("rawr",        d.path());
  EXPECT_EQ("file://rawr", c.string());
  EXPECT_EQ("file://rawr", d.string());

  Url e("pants/cheese");
  Url f("file://pants/cheese");

  EXPECT_EQ("pants/cheese",        e.path());
  EXPECT_EQ("pants/cheese",        f.path());
  EXPECT_EQ("file://pants/cheese", e.string());
  EXPECT_EQ("file://pants/cheese", f.string());
}

TEST(HTTPUtils, Tricky) {
  Url u("/pants:moo/");
  EXPECT_EQ("file", u.scheme());
  EXPECT_EQ("/pants:moo/", u.path());
}

TEST(HTTPUtils, PF) {
  Url u1("pf:///exchange/platefilename.plate");
  EXPECT_EQ("pf", u1.scheme());
  EXPECT_TRUE(u1.netloc().empty());
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
  EXPECT_NO_THROW(u.netloc("woof"));
  EXPECT_NO_THROW(u.path(u.path() + "/append"));
  EXPECT_EQ("pf://woof/pants/cheese/append?rawr=foo", u.string());
}

TEST(HTTPUtils, Errors) {
  Url u;
  EXPECT_THROW(u.string(),     LogicErr);
  EXPECT_THROW(u.path(""),     ArgumentErr);
  EXPECT_THROW(u.path("rawr"), ArgumentErr);

  EXPECT_THROW(Url(""),     ArgumentErr);
}

TEST(HTTPUtils, Operator) {
  const char data[] = "http://rawr/pants?cheese=moo&pants=foo#wee=blarg";
  std::istringstream in(data);
  std::ostringstream out;

  Url u1;

  in  >> u1;
  out << u1;

  EXPECT_EQ(data, u1.string());
  EXPECT_EQ(data, out.str());
}

TEST(HTTPUtils, PathSplit) {
  // Absolute
  //   "/"             -> [""]
  //   "/pants"        -> ["", "pants"]
  //   "/pants/cheese" -> ["", "pants", "cheese"]
  // Relative
  //   ""              -> [] (ASSERTION FAILED)
  //   "pants"         -> ["pants"]
  //   "pants/cheese"  -> ["pants", "cheese"]

  Url::split_t split;

  split = Url("/").path_split();
  EXPECT_EQ(1, split.size());
  EXPECT_EQ("", split[0]);

  split = Url("/pants").path_split();
  EXPECT_EQ(2, split.size());
  EXPECT_EQ("", split[0]);
  EXPECT_EQ("pants", split[1]);

  split = Url("/pants/cheese").path_split();
  EXPECT_EQ(3, split.size());
  EXPECT_EQ("", split[0]);
  EXPECT_EQ("pants", split[1]);
  EXPECT_EQ("cheese", split[2]);

  EXPECT_THROW(split = Url("").path_split(), ArgumentErr);

  split = Url("pants").path_split();
  EXPECT_EQ(1, split.size());
  EXPECT_EQ("pants", split[0]);

  split = Url("pants/cheese").path_split();
  EXPECT_EQ(2, split.size());
  EXPECT_EQ("pants", split[0]);
  EXPECT_EQ("cheese", split[1]);
}

struct range_t {
  typedef std::string* T;
  typedef boost::iterator_range<T> r_t;
  const T data;
  range_t(const T data_) : data(data_) {}
  r_t operator()(size_t b, size_t e) {return r_t(data + b, data + e);}
};

TEST(HTTPUtils, PathJoin) {
  std::string A[] = {"", "pants", "cheese"};
  range_t r(A);

  Url j;
  j.scheme("file");

  j.path_join(r(0,1));
  EXPECT_EQ("/", j.path());
  EXPECT_EQ("file:///", j.string());

  j.path_join(r(0,2));
  EXPECT_EQ("/pants", j.path());
  EXPECT_EQ("file:///pants", j.string());

  j.path_join(r(0,3));
  EXPECT_EQ("/pants/cheese", j.path());
  EXPECT_EQ("file:///pants/cheese", j.string());

  EXPECT_THROW(j.path_join(r(1,1)), ArgumentErr);

  j.path_join(r(1,2));
  EXPECT_EQ("pants", j.path());
  EXPECT_EQ("file://pants", j.string());

  j.path_join(r(1,3));
  EXPECT_EQ("pants/cheese", j.path());
  EXPECT_EQ("file://pants/cheese", j.string());
}

TEST(HTTPUtils, PlatefileUrl) {
  EXPECT_THROW(PlatefileUrl u("http://rawr", "pants"), ArgumentErr);

  PlatefileUrl u("http://rawr", "pants.plate");
  EXPECT_THROW(u.name("pants2"), ArgumentErr);

  EXPECT_EQ("http://rawr/pants.plate", u.string());
  EXPECT_EQ("pants.plate", u.name());
  EXPECT_EQ("http://rawr/", u.base().string());

  EXPECT_NO_THROW(u.name("pants2.plate"));
  EXPECT_EQ("http://rawr/pants2.plate", u.string());
  EXPECT_EQ("pants2.plate", u.name());
  EXPECT_EQ("http://rawr/", u.base().string());

  EXPECT_NO_THROW(u.base("moo://moo?bob=waffles"));
  EXPECT_EQ("moo://moo/pants2.plate?bob=waffles", u.string());
  EXPECT_EQ("pants2.plate", u.name());
  EXPECT_EQ("moo://moo/?bob=waffles", u.base().string());
}
