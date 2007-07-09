
uniform sampler2DRect image;

void main() {
   vec2 coord = gl_TexCoord[0].st;

   vec4 cPixel = texture2DRect(image, coord);
   float cMin = cPixel.r;
   cMin = min(cMin, cPixel.g);
   cMin = min(cMin, cPixel.b);
   cMin = min(cMin, cPixel.a);

   gl_FragColor.r = gl_TexCoord[0].t; //cMin;
}
