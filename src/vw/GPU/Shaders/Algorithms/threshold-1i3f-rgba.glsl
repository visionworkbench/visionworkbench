// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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

