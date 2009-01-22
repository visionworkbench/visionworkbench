
uniform sampler2DRect i1;
uniform sampler2DRect i2;

void main() {
   gl_FragColor.rgba = pow(texture2DRect(i1, gl_TexCoord[0].st), texture2DRect(i2, gl_TexCoord[0].st));
}