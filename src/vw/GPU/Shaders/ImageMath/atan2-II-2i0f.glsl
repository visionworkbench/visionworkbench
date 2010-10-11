// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


uniform sampler2DRect i1;
uniform sampler2DRect i2;

void main() {
   float pi = 3.141592654;
   float half_pi = 3.141592654 / 2.0;

   vec4 y = texture2DRect(i1, gl_TexCoord[0].st);
   vec4 x = texture2DRect(i2, gl_TexCoord[0].st);
   vec4 result = vec4(0.0);
   vec4 atan_value = atan(y/x);
/* Per-Channel Variables */
   float ch_x;
   float ch_y;
   float ch_atan;
   float ch_output;
/* channel-r */
   ch_x = x.r;
   ch_y = y.r;
   ch_atan = atan_value.r;
   ch_output = ch_atan;
   if(ch_x == 0.0) {
      ch_output = half_pi * sign(ch_y);
   }
   else if(ch_x < 0.0) {
      ch_output = ch_atan + sign(ch_y) * pi;
   }
        result.r = ch_output;
/* channel-g */
   ch_x = x.g;
   ch_y = y.g;
   ch_atan = atan_value.g;
   ch_output = ch_atan;
   if(ch_x == 0.0) {
      ch_output = half_pi * sign(ch_y);
   }
   else if(ch_x < 0.0) {
      ch_output = ch_atan + sign(ch_y) * pi;
   }
        result.g = ch_output;
/* channel-b */
   ch_x = x.b;
   ch_y = y.b;
   ch_atan = atan_value.b;
   ch_output = ch_atan;
   if(ch_x == 0.0) {
      ch_output = half_pi * sign(ch_y);
   }
   else if(ch_x < 0.0) {
      ch_output = ch_atan + sign(ch_y) * pi;
   }
        result.b = ch_output;
/* channel-a */
   ch_x = x.a;
   ch_y = y.a;
   ch_atan = atan_value.a;
   ch_output = ch_atan;
   if(ch_x == 0.0) {
      ch_output = half_pi * sign(ch_y);
   }
   else if(ch_x < 0.0) {
      ch_output = ch_atan + sign(ch_y) * pi;
   }
        result.a = ch_output;

   gl_FragColor.rgba = result;
}


float sign(float n) {
   if(n < 0.0) {
      return -1.0;
   }
   else {
      return 1.0;
   }
}
