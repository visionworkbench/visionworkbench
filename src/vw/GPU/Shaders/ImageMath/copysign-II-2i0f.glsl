// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


uniform sampler2DRect i1;
uniform sampler2DRect i2;

void main() {
   vec4 value = texture2DRect(i1, gl_TexCoord[0].st);
   vec4 ref = texture2DRect(i2, gl_TexCoord[0].st);

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


