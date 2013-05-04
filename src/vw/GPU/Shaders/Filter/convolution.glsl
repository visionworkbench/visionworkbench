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
uniform sampler2DRect kernel;
uniform float hHalfSize;
uniform float vHalfSize;

void main() {
 int hSize = $1;
 int vSize = $2;
 vec2 startCoord = gl_TexCoord[0].st + vec2(hHalfSize, vHalfSize);
 vec4 sum = vec4(0.0);

 for(int vKernel = 0; vKernel < hSize; vKernel++) {
    for(int hKernel = 0; hKernel < vSize; hKernel++) {
           sum += texture2DRect(kernel, vec2(float(hKernel) + 0.001, float(vKernel) + 0.001)).r
                          * texture2DRect(image, vec2(startCoord.s - (float(hKernel) - 0.001), startCoord.t - (float(vKernel) - 0.001))).rgba;
    }
 }
 gl_FragColor.rgba = sum;
}
