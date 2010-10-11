// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


uniform sampler2DRect image;

void main() {
   vec2 coord = gl_TexCoord[0].st;

   vec4 cPixel = texture2DRect(image, coord);
   float sum = cPixel.r;

   gl_FragColor.r = sum;
}
