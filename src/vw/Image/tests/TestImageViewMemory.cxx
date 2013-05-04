// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


#include <test/Helpers.h>
#include <vw/Image/ImageView.h>

using namespace vw;
using namespace vw::test;
using namespace std;

// spike operator new... every nothrow new fails [this lets us ask for
// unlimited memory without the system actually trying to give it]
void* operator new[](size_t, const nothrow_t&) throw() {
  return 0;
}

void operator delete[](void*, const nothrow_t&) throw() {
  /* noop */
}

TEST(ImageView, TooMuchMemory) {
  int32 cols, rows, planes;
  cols = rows = 1<<15;
  planes = 1 << 9;

  ImageView<uint8> img;

  // We can probably only trip this error on 32-bit platforms
  if (sizeof(size_t) <= sizeof(int32))
    EXPECT_THROW_MSG(img.set_size(cols, rows, planes), ArgumentErr, "too many pixels");

  rows = planes = 1;
  EXPECT_THROW_MSG(img.set_size(cols, rows, planes), ArgumentErr, "too many bytes");
}
