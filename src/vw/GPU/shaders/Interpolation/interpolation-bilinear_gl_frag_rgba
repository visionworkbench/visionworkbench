
uniform sampler2DRect image;

void main() {
   vec2 coord = gl_TexCoord[0].st;
  // coord = coord - vec2(0.5, 0.5);
   vec2 coord_fract = fract(coord);
   coord = floor(coord) + vec2(0.5, 0.5);
   gl_FragColor.rgba = 
	texture2DRect(image, coord) * (1.0 - coord_fract.s) * (1.0 - coord_fract.t) +
	texture2DRect(image, coord + vec2(1.0, 0.0)) * coord_fract.s * (1.0 - coord_fract.t) +
	texture2DRect(image, coord + vec2(0.0, 1.0)) * (1.0 - coord_fract.s) * coord_fract.t +
	texture2DRect(image, coord + vec2(1.0, 1.0)) * coord_fract.s * coord_fract.t;
  //gl_FragColor = vec4(gl_TexCoord[0].t);
}
