#include <iostream>

#include <vw/config.h>

#include <vw/vw.h>

using namespace vw;

int main() {
  std::cout << "Hello, world!" << std::endl;

  PixelRGB<float> p1(0.2,0.3,0.4);
  PixelRGB<int> p2(1,2,3);
  std::cout << p1 << std::endl;
  std::cout << p1+p2 << std::endl;
  std::cout << p1-p2 << std::endl;
  std::cout << 2*p1 << std::endl;
  std::cout << p1*2 << std::endl;
  std::cout << p2/2.0 << std::endl;
  std::cout << -p1 << std::endl;
  std::cout << -p2 << std::endl;
  float x = 10.0;

  std::cout << _numeric(ChannelRange<int8>::max()) << std::endl;
  std::cout << _numeric(ChannelRange<int16>::max()) << std::endl;
  std::cout << _numeric(ChannelRange<int32>::max()) << std::endl;
  std::cout << _numeric(ChannelRange<int64>::max()) << std::endl;
  std::cout << _numeric(ChannelRange<uint8>::max()) << std::endl;
  std::cout << _numeric(ChannelRange<uint16>::max()) << std::endl;
  std::cout << _numeric(ChannelRange<uint32>::max()) << std::endl;
  std::cout << _numeric(ChannelRange<uint64>::max()) << std::endl;
  std::cout << _numeric(ChannelRange<float32>::max()) << std::endl;
  std::cout << _numeric(ChannelRange<float64>::max()) << std::endl;
  std::cout << _numeric(ChannelRange<PixelRGB<std::complex<uint8> > >::max()) << std::endl;
  std::cout << p1 << std::endl;
  std::cout << channel_cast<uint8>(10*p1) << std::endl;
  std::cout << _numeric(channel_cast<uint8>(x)) << std::endl;
  std::cout << channel_cast_rescale<uint8>(p1) << std::endl;
  return 0;
}
