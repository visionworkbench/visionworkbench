
uniform sampler2DRect i1;

void main() {
   vec4 value = texture2DRect(i1, gl_TexCoord[0].st);
   vec4 e = vec4(2.7182818284590451);
   gl_FragColor.rgba = (pow(e, value) + pow(e, -1.0 * value)) / 2.0;
}
 