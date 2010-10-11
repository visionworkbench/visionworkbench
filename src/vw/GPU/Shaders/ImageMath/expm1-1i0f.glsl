// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


uniform sampler2DRect i1;

void main() {
   gl_FragColor.rgba = pow(vec4(2.7182818284590451), texture2DRect(i1, gl_TexCoord[0].st)) - 1.0;
}
