
uniform sampler2DRect i1;

void main() {
   vec4 value = texture2DRect(i1, gl_TexCoord[0].st);
   vec4 rounded = floor(value);
   if(value.r - rounded.r >= 0.5) {
      rounded.r += 1.0;
   }
   if(value.g - rounded.g >= 0.5) { 
      rounded.g += 1.0;
   }
   if(value.b - rounded.b >= 0.5) {
      rounded.b += 1.0;
   }
   if(value.a - rounded.a >= 0.5) {
      rounded.a += 1.0;
   }
   gl_FragColor.rgba = rounded;
}
 