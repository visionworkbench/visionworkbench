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


uniform sampler2DRect image;

void main() {
   vec2 ij = gl_TexCoord[0].st;
   vec2        xy = floor(ij);
   vec2 normxy = ij - xy;
   vec2 st0 = ((2.0 - normxy) * normxy - 1.0) * normxy;
   vec2 st1 = (3.0 * normxy - 5.0) * normxy * normxy + 2.0;
   vec2 st2 = ((4.0 - 3.0 * normxy) * normxy + 1.0) * normxy;
   vec2 st3 = (normxy - 1.0) * normxy * normxy;

   vec4 row0 =
        st0.s * texture2DRect(image, xy + vec2(-1.0, -1.0)) +
        st1.s * texture2DRect(image, xy + vec2(0.0, -1.0)) +
        st2.s * texture2DRect(image, xy + vec2(1.0, -1.0)) +
        st3.s * texture2DRect(image, xy + vec2(2.0, -1.0));

   vec4 row1 =
        st0.s * texture2DRect(image, xy + vec2(-1.0, 0.0)) +
        st1.s * texture2DRect(image, xy + vec2(0.0, 0.0)) +
        st2.s * texture2DRect(image, xy + vec2(1.0, 0.0)) +
        st3.s * texture2DRect(image, xy + vec2(2.0, 0.0));

   vec4 row2 =
        st0.s * texture2DRect(image, xy + vec2(-1.0, 1.0)) +
        st1.s * texture2DRect(image, xy + vec2(0.0, 1.0)) +
        st2.s * texture2DRect(image, xy + vec2(1.0, 1.0)) +
        st3.s * texture2DRect(image, xy + vec2(2.0, 1.0));

   vec4 row3 =
        st0.s * texture2DRect(image, xy + vec2(-1.0, 2.0)) +
        st1.s * texture2DRect(image, xy + vec2(0.0, 2.0)) +
        st2.s * texture2DRect(image, xy + vec2(1.0, 2.0)) +
        st3.s * texture2DRect(image, xy + vec2(2.0, 2.0));

   gl_FragColor.rgba = 0.25 * ((st0.t * row0) + (st1.t * row1) + (st2.t * row2) + (st3.t * row3));
}
