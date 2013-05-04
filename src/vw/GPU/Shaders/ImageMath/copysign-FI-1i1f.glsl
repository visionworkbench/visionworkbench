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
   vec4 value = vec4(f1);
   vec4 ref = texture2DRect(i1, gl_TexCoord[0].st);

   if(value.r * ref.r < 0.0)
      value.r *= -1.0;
   if(value.g * ref.g < 0.0)
      value.g *= -1.0;
   if(value.b * ref.b < 0.0)
      value.b *= -1.0;
   if(value.a * ref.a < 0.0)
      value.a *= -1.0;

   gl_FragColor.rgba = value;
}


