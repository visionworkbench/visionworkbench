
uniform sampler2DRect i1;

void main() {
   vec4 value = texture2DRect(i1, gl_TexCoord[0].st);
   gl_FragColor.rgba = log2(value + sqrt(pow(value, vec4(2.0)) - 1.0)) / 1.442695041;
}
 