
uniform sampler2DRect i1;
uniform float f1;  // threshold
uniform float f2;  // low
uniform float f3;  // high

void main() {
   vec4 value = texture2DRect(i1, gl_TexCoord[0].st);

   if(value.r > f1)
      gl_FragColor.r = f3;
   else
      gl_FragColor.r = f2;

}
  