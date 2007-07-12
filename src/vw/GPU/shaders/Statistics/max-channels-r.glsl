
uniform sampler2DRect image;

void main() {
   vec2 coord = gl_TexCoord[0].st;

   vec4 cPixel = texture2DRect(image, coord);
   float cMax = cPixel.r;

   gl_FragColor.r = cMax;
}
