
uniform sampler2DRect i1;

void main() {
   gl_FragColor.rgba = pow(vec4(2.7182818284590451), texture2DRect(i1, gl_TexCoord[0].st));
}
 