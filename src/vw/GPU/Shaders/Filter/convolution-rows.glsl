// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


uniform sampler2DRect image;
uniform sampler2DRect kernel;
uniform float halfSize;

void main() {
 int size = $1;
 vec2 startCoord = gl_TexCoord[0].st + vec2(halfSize, 0.0);
 vec4 sum = vec4(0.0);
 for(int vectorPos = 0; vectorPos < size; vectorPos++) {
        float vectorValue = texture2DRect(kernel, vec2(float(vectorPos) + 0.5, 0.0)).r;
        sum += vectorValue * (texture2DRect(image, vec2(startCoord.s - float(vectorPos), startCoord.t))).rgba;
 }
 gl_FragColor.rgba = sum;
}
