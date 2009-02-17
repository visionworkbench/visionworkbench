// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


uniform sampler2DRect image;

void main() {
   vec2 coord = gl_TexCoord[0].st;

   vec4 cPixel = texture2DRect(image, coord);
   float cMax = cPixel.r;
   cMax = max(cMax, cPixel.g);
   cMax = max(cMax, cPixel.b);
   cMax = max(cMax, cPixel.a);

   gl_FragColor.r = cMax;
}
