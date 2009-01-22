
uniform sampler2DRect i1;

void main() {
	vec4 value = texture2DRect(i1, gl_TexCoord[0].st);
	vec4 result = floor(value);
	if(value.r < 0.0) result.r += 1.0;
	if(value.g < 0.0) result.g += 1.0;
	if(value.b < 0.0) result.b += 1.0;
	if(value.a < 0.0) result.a += 1.0;
	gl_FragColor.rgba = result;
}
 