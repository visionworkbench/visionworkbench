// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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
