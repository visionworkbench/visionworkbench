
uniform sampler2DRect image;

void main() {
   vec2 coord = gl_TexCoord[0].st; // - vec2(1.0, 1.0);

   vec4 cPixel = texture2DRect(image, coord);
   float sum = cPixel.r;

   cPixel = texture2DRect(image, coord + vec2(1.0, 0.0));
   sum += cPixel.r;

   cPixel = texture2DRect(image, coord + vec2(0.0, 1.0));
   sum += cPixel.r;

   cPixel = texture2DRect(image, coord + vec2(1.0, 1.0));
   sum += cPixel.r;
		
   gl_FragColor.r = sum;
}
