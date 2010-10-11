// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


uniform sampler2DRect i1;

void main() {
   gl_FragColor.rgba = pow(texture2DRect(i1, gl_TexCoord[0].st), vec4(1.0/3.0));
}
