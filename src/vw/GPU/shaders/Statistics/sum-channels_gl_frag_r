
uniform sampler2DRect image;

void main() {
   vec2 coord = gl_TexCoord[0].st;

   vec4 cPixel = texture2DRect(image, coord);
   float sum = cPixel.r;
		
   gl_FragColor.r = sum;
}
