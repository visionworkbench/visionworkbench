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
uniform float f1;  /* threshold */
uniform float f2;  /* low */
uniform float f3;  /* high */

void main() {
   vec4 value = texture2DRect(i1, gl_TexCoord[0].st);
   if(value.r > f1)
      gl_FragColor.r = f3;
   else
      gl_FragColor.r = f2;

   if(value.g > f1)
      gl_FragColor.g = f3;
   else
      gl_FragColor.g = f2;

   if(value.b > f1)
      gl_FragColor.b = f3;
   else
      gl_FragColor.b = f2;

   if(value.a > f1)
      gl_FragColor.a = f3;
   else
      gl_FragColor.a = f2;
}

