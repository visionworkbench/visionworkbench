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


uniform sampler2DRect i1;
uniform float f1;

void main() {
   float pi = 3.141592654;
   float half_pi = 3.141592654 / 2.0;

   vec4 y = texture2DRect(i1, gl_TexCoord[0].st);
   vec4 x = vec4(f1);
   vec4 result;
         vec4 atan_value = atan(y/x);
/* channel-r */
   if(x.r > 0.0) {
      result.r = atan_value.r;
   }
   else if(x.r < 0.0) {
      result.r = atan_value.r + sign(y.r) * pi;
   }
   else if(y.r == 0.0) {
      result.r = atan_value.r;
   }
   else {
      result.r = half_pi * sign(y.r);
   }
/* channel-g */
   if(x.g > 0.0) {
      result.g = atan_value.g;
   }
   else if(x.g < 0.0) {
      result.g = atan_value.g + sign(y.g) * pi;
   }
   else if(y.g == 0.0) {
      result.g = atan_value.g;
   }
   else {
      result.g = half_pi * sign(y.g);
   }
/* channel-b */
   if(x.b > 0.0) {
      result.b = atan_value.b;
   }
   else if(x.b < 0.0) {
      result.b = atan_value.b + sign(y.b) * pi;
   }
   else if(y.b == 0.0) {
      result.b = atan_value.b;
   }
   else {
      result.b = half_pi * sign(y.b);
   }
/* channel-a */
   if(x.a > 0.0) {
      result.a = atan_value.a;
   }
   else if(x.a < 0.0) {
      result.a = atan_value.a + sign(y.a) * pi;
   }
   else if(y.a == 0.0) {
      result.a = atan_value.a;
   }
   else {
      result.a = half_pi * sign(y.a);
   }

   gl_FragColor.rgba = result;
}


float sign(float n) {
   if(n < 0.0) {
      return -1.0;
   }
   else {
      return 1.0;
   }
}
