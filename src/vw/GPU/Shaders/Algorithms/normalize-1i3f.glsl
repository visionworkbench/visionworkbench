// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


uniform sampler2DRect i1;
uniform float f1;  /* old_min */
uniform float f2;  /* new_min */
uniform float f3;  /* new_to_old_ratio */

void main() {
   vec4 value = texture2DRect(i1, gl_TexCoord[0].st);
   gl_FragColor.rgba = vec4(f2) + vec4(f3) * (value - vec4(f1));
}

