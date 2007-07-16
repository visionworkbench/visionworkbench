
uniform sampler2DRect i1;

void main() {
   gl_FragColor.rgba = texture2DRect(i1, gl_TexCoord[0].st);
}
