
uniform sampler2DRect i1;

void main() {
   vec4 value = texture2DRect(i1, gl_TexCoord[0].st);
   gl_FragColor.rgba = vec4(0.5) * (log2((1.0 + value) / (1.0 - value)) / vec4(1.442695041));
}
 