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

void main() {
   vec4 value = texture2DRect(i1, gl_TexCoord[0].st);
   vec4 rounded = floor(value);
   if(value.r - rounded.r >= 0.5) {
      rounded.r += 1.0;
   }
   if(value.g - rounded.g >= 0.5) {
      rounded.g += 1.0;
   }
   if(value.b - rounded.b >= 0.5) {
      rounded.b += 1.0;
   }
   if(value.a - rounded.a >= 0.5) {
      rounded.a += 1.0;
   }
   gl_FragColor.rgba = rounded;
}
