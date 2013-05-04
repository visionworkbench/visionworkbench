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
uniform sampler2DRect i2;

void main() {
   float pi = 3.141592654;
   float half_pi = 3.141592654 / 2.0;

   vec4 y = texture2DRect(i1, gl_TexCoord[0].st);
   vec4 x = texture2DRect(i2, gl_TexCoord[0].st);
   vec4 result = vec4(0.0);
   vec4 atan_value = atan(y/x);
/* Per-Channel Variables */
   float ch_x;
   float ch_y;
   float ch_atan;
   float ch_output;
/* channel-r */
   ch_x = x.r;
   ch_y = y.r;
   ch_atan = atan_value.r;
   ch_output = ch_atan;
   if(ch_x == 0.0) {
      ch_output = half_pi * sign(ch_y);
   }
   else if(ch_x < 0.0) {
      ch_output = ch_atan + sign(ch_y) * pi;
   }
        result.r = ch_output;
/* channel-g */
   ch_x = x.g;
   ch_y = y.g;
   ch_atan = atan_value.g;
   ch_output = ch_atan;
   if(ch_x == 0.0) {
      ch_output = half_pi * sign(ch_y);
   }
   else if(ch_x < 0.0) {
      ch_output = ch_atan + sign(ch_y) * pi;
   }
        result.g = ch_output;
/* channel-b */
   ch_x = x.b;
   ch_y = y.b;
   ch_atan = atan_value.b;
   ch_output = ch_atan;
   if(ch_x == 0.0) {
      ch_output = half_pi * sign(ch_y);
   }
   else if(ch_x < 0.0) {
      ch_output = ch_atan + sign(ch_y) * pi;
   }
        result.b = ch_output;
/* channel-a */
   ch_x = x.a;
   ch_y = y.a;
   ch_atan = atan_value.a;
   ch_output = ch_atan;
   if(ch_x == 0.0) {
      ch_output = half_pi * sign(ch_y);
   }
   else if(ch_x < 0.0) {
      ch_output = ch_atan + sign(ch_y) * pi;
   }
        result.a = ch_output;

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
