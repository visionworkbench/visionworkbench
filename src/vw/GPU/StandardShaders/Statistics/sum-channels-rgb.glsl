
uniform sampler2DRect image;

void main() {
   vec2 coord = gl_TexCoord[0].st;

   vec4 cPixel = texture2DRect(image, coord);
   float sum = cPixel.r;
   sum += cPixel.g;
   sum += cPixel.b;

   gl_FragColor.r = sum;
}
