
uniform sampler2DRect i1;
uniform float f1;

void main() {
   gl_FragColor.rgba = sqrt(pow(texture2DRect(i1, gl_TexCoord[0].st), vec4(2.0)) + pow(vec4(f1), vec4(2.0)));
}
