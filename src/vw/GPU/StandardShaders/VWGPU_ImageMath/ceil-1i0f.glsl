
uniform sampler2DRect i1;

void main() {
   gl_FragColor.rgba = ceil(texture2DRect(i1, gl_TexCoord[0].st));
}
 