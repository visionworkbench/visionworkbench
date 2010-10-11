// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


uniform sampler2DRect i1;

void main() {
   vec4 value = texture2DRect(i1, gl_TexCoord[0].st);
   vec4 e = vec4(2.7182818284590451);
   vec4 e2x = pow(e, 2.0 * value);
   gl_FragColor.rgba = (e2x - 1.0) / (e2x + 1.0);
}
