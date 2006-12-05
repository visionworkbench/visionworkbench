
uniform sampler2DRect image;

void main() {
   vec2 coord = gl_TexCoord[0].st;
   vec2 coord_round = floor(coord);
   gl_FragColor.rgba = texture2DRect(image, coord_round);
}
