// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


uniform sampler2DRect i1;
uniform float xCellOffset;
uniform float yCellOffset;

void main() {
   gl_FragColor.r = texture2DRect(i1, gl_TexCoord[0].st).r;
   gl_FragColor.g = texture2DRect(i1, gl_TexCoord[0].st + vec2(xCellOffset, 0)).r;
   gl_FragColor.b = texture2DRect(i1, gl_TexCoord[0].st + vec2(0, yCellOffset)).r;
   gl_FragColor.a = texture2DRect(i1, gl_TexCoord[0].st + vec2(xCellOffset, yCellOffset)).r;
}
