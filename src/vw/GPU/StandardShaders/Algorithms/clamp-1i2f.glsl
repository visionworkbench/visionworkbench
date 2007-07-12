
uniform sampler2DRect i1;
uniform float f1;
uniform float f2;


void main() {
   gl_FragColor.rgba = clamp(texture2DRect(i1, gl_TexCoord[0].st), f1, f2);
}