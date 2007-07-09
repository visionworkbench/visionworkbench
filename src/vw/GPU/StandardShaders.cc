#include <vw/GPU/StandardShaders.h>

namespace vw { namespace GPU {

std::map<std::string, char*> standard_shaders_map;

void init_standard_shaders() {
standard_shaders_map["./StandardShaders/Filter/convolution-columns_cg_frag_rgba"] = " \
 \
void main( \
      in float2 texCoord0 : TEXCOORD0, \
      uniform samplerRECT image, \
      uniform samplerRECT kernel, \
      uniform float halfSize, \
      out float4 color : COLOR ) \
  { \
     float2 startCoord = float2(texCoord0[0], texCoord0[1] + halfSize); \
     int size = $1; \
     color = float4(0.0, 0.0, 0.0, 0.0); \
     for(int kernelPos = 0; kernelPos < size; kernelPos++) { \
        float kernelValue = texRECT(kernel, float2(kernelPos, 0)).r; \
	color += kernelValue * texRECT(image, float2(startCoord[0], startCoord[1] - kernelPos)); \
     } \
 } \
 \
";

standard_shaders_map["./StandardShaders/Filter/convolution-columns_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image; \
uniform sampler2DRect kernel; \
uniform float halfSize; \
 \
void main() { \
 int size = $1; \
 vec2 startCoord = gl_TexCoord[0].st + vec2(0.0, halfSize); \
 vec4 sum = vec4(0.0); \
 for(int vectorPos = 0; vectorPos < size; vectorPos++) { \
	float vectorValue = texture2DRect(kernel, vec2(float(vectorPos), 0.0)).r; \
	sum += vectorValue * (texture2DRect(image,  vec2(startCoord.s + 0.0, startCoord.t - float(vectorPos)))).rgba; \
 } \
 gl_FragColor.rgba = sum + 0.1; \
} \
";

standard_shaders_map["./StandardShaders/Filter/convolution-rows_cg_frag_rgba"] = " \
 \
void main( \
      in float2 texCoord0 : TEXCOORD0, \
      uniform samplerRECT image, \
      uniform samplerRECT kernel, \
      uniform float halfSize, \
      out float4 color : COLOR ) \
  { \
     float2 startCoord = float2(texCoord0[0] + halfSize, texCoord0[1]); \
     int size = $1; \
     color = float4(0.0, 0.0, 0.0, 0.0); \
     for(int kernelPos = 0; kernelPos < size; kernelPos++) { \
        float kernelValue = texRECT(kernel, float2(kernelPos, 0)).r; \
	color += kernelValue * texRECT(image, float2(startCoord[0] - kernelPos, startCoord[1])); \
     } \
 } \
 \
";

standard_shaders_map["./StandardShaders/Filter/convolution-rows_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image; \
uniform sampler2DRect kernel; \
uniform float halfSize; \
 \
void main() { \
 int size = $1; \
 vec2 startCoord = gl_TexCoord[0].st + vec2(halfSize, 0.0); \
 vec4 sum = vec4(0.0); \
 for(int vectorPos = 0; vectorPos < size; vectorPos++) { \
	float vectorValue = texture2DRect(kernel, vec2(float(vectorPos), 0.0)).r; \
	sum += vectorValue * (texture2DRect(image, vec2(startCoord.s - float(vectorPos), startCoord.t))).rgba; \
 } \
 gl_FragColor.rgba = sum; \
} \
";

standard_shaders_map["./StandardShaders/Filter/convolution_cg_frag_rgba"] = " \
void main( \
      in float2 texCoord0 : TEXCOORD0, \
      uniform samplerRECT image, \
      uniform samplerRECT kernel, \
      uniform float hHalfSize, \
      uniform float vHalfSize, \
      out float4 color : COLOR ) \
  { \
     float2 startCoord = float2(texCoord0[0] + hHalfSize, texCoord0[1] + vHalfSize); \
     int hSize = $1; \
     int vSize = $2; \
     color = float4(0.0, 0.0, 0.0, 0.0); \
     for(int vKernel = 0; vKernel < vSize; vKernel++) { \
        for(int hKernel = 0; hKernel < hSize; hKernel++) { \
	   float kernelValue = texRECT(kernel, float2(hKernel, vKernel)).r; \
	   color += kernelValue * texRECT(image, float2(startCoord[0] - hKernel, startCoord[1] - vKernel)); \
       } \
    } \
 } \
 \
 \
";

standard_shaders_map["./StandardShaders/Filter/convolution_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image; \
uniform sampler2DRect kernel; \
uniform float hHalfSize; \
uniform float vHalfSize; \
 \
void main() { \
 int hSize = $1; \
 int vSize = $2; \
 vec2 startCoord = gl_TexCoord[0].st + vec2(hHalfSize, vHalfSize); \
 vec4 sum = vec4(0.0); \
 for(int vKernel = 0; vKernel < vSize; vKernel++) { \
    for(int hKernel = 0; hKernel < hSize; hKernel++) { \
	float kernelValue = texture2DRect(kernel, vec2(float(hKernel), float(vKernel))).r; \
	sum += kernelValue * texture2DRect(image, vec2(startCoord.s - float(hKernel), startCoord.t - float(vKernel))).rgba; \
    } \
 } \
 gl_FragColor.rgba = sum; \
} \
";

standard_shaders_map["./StandardShaders/Identity_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT image, \
	out float4 color : COLOR )  \
{ \
   color = texRECT(image, texCoord0); \
} \
";

standard_shaders_map["./StandardShaders/Identity_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image; \
 \
void main() { \
   gl_FragColor.rgba = texture2DRect(image, gl_TexCoord[0].st); \
} \
";

standard_shaders_map["./StandardShaders/Algorithmsclamp-1i2f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform float f1, \
	uniform float f2, \
	out float4 color : COLOR )  \
{ \
   color = clamp(texRECT(i1, texCoord0), f1, f2); \
} \
";

standard_shaders_map["./StandardShaders/Algorithmsclamp-1i2f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
uniform float f1; \
uniform float f2; \
 \
 \
void main() { \
   gl_FragColor.rgba = clamp(texture2DRect(i1, gl_TexCoord[0].st), f1, f2); \
} \
";

standard_shaders_map["./StandardShaders/Algorithmscopy-1i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	out float4 color : COLOR )  \
{ \
   color = texRECT(i1, texCoord0); \
} \
";

standard_shaders_map["./StandardShaders/Algorithmscopy-1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   gl_FragColor.rgba = texture2DRect(i1, gl_TexCoord[0].st); \
} \
";

standard_shaders_map["./StandardShaders/Algorithmsfill-1i4f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform float f1, \
	uniform float f2, \
	uniform float f3, \
	uniform float f4, \
	out float4 color : COLOR )  \
{ \
   color = float4(f1, f2, f3, f4); \
} \
 \
";

standard_shaders_map["./StandardShaders/Algorithmsfill-1i4f_gl_frag_rgba"] = " \
 \
uniform float f1; \
uniform float f2; \
uniform float f3; \
uniform float f4; \
 \
void main() { \
   gl_FragColor.rgba = vec4(f1, f2, f3, f4); \
} \
";

standard_shaders_map["./StandardShaders/Algorithmsnormalize-1i3f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform float f1,  // old_min \
	uniform float f2,  // new_min \
	uniform float f3,  // new_to_old_ratio \
	out float4 color : COLOR )  \
{ \
   color = (texRECT(i1, texCoord0) - f1) * f3 + f2; \
} \
 \
";

standard_shaders_map["./StandardShaders/Algorithmsnormalize-1i3f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
uniform float f1;  // old_min \
uniform float f2;  // new_min \
uniform float f3;  // new_to_old_ratio \
 \
void main() { \
   vec4 value = texture2DRect(i1, gl_TexCoord[0].st); \
   gl_FragColor.rgba = (value - vec4(f1) * vec4(f3)) + vec4(f2); \
} \
  \
";

standard_shaders_map["./StandardShaders/Algorithmsthreshold-1i3f_cg_frag_r"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform float f1,  // threshold \
	uniform float f2,  // low \
	uniform float f3,  // high \
	out float4 color : COLOR )  \
{ \
   float4 value = texRECT(i1, texCoord0); \
 \
   if(value.r > f1) \
      color.r = f3; \
   else \
      color.r = f2; \
} \
";

standard_shaders_map["./StandardShaders/Algorithmsthreshold-1i3f_cg_frag_rgb"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform float f1,  // threshold \
	uniform float f2,  // low \
	uniform float f3,  // high \
	out float4 color : COLOR )  \
{ \
   float4 value = texRECT(i1, texCoord0); \
 \
   if(value.r > f1) \
      color.r = f3; \
   else \
      color.r = f2; \
 \
   if(value.g > f1) \
      color.g = f3; \
   else \
      color.g = f2; \
 \
   if(value.b > f1) \
      color.b = f3; \
   else \
      color.b = f2; \
 \
} \
";

standard_shaders_map["./StandardShaders/Algorithmsthreshold-1i3f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform float f1,  // threshold \
	uniform float f2,  // low \
	uniform float f3,  // high \
	out float4 color : COLOR )  \
{ \
   float4 value = texRECT(i1, texCoord0); \
 \
   if(value.r > f1) \
      color.r = f3; \
   else \
      color.r = f2; \
 \
   if(value.g > f1) \
      color.g = f3; \
   else \
      color.g = f2; \
 \
   if(value.b > f1) \
      color.b = f3; \
   else \
      color.b = f2; \
 \
   if(value.a > f1) \
      color.a = f3; \
   else \
      color.a = f2; \
 \
} \
";

standard_shaders_map["./StandardShaders/Algorithmsthreshold-1i3f_gl_frag_r"] = " \
 \
uniform sampler2DRect i1; \
uniform float f1;  // threshold \
uniform float f2;  // low \
uniform float f3;  // high \
 \
void main() { \
   vec4 value = texture2DRect(i1, gl_TexCoord[0].st); \
 \
   if(value.r > f1) \
      gl_FragColor.r = f3; \
   else \
      gl_FragColor.r = f2; \
 \
} \
   \
";

standard_shaders_map["./StandardShaders/Algorithmsthreshold-1i3f_gl_frag_rgb"] = " \
 \
uniform sampler2DRect i1; \
uniform float f1;  // threshold \
uniform float f2;  // low \
uniform float f3;  // high \
 \
void main() { \
   vec4 value = texture2DRect(i1, gl_TexCoord[0].st); \
 \
   if(value.r > f1) \
      gl_FragColor.r = f3; \
   else \
      gl_FragColor.r = f2; \
 \
   if(value.g > f1) \
      gl_FragColor.g = f3; \
   else \
      gl_FragColor.g = f2; \
 \
   if(value.b > f1) \
      gl_FragColor.b = f3; \
   else \
      gl_FragColor.b = f2; \
 \
} \
  \
";

standard_shaders_map["./StandardShaders/Algorithmsthreshold-1i3f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
uniform float f1;  // threshold \
uniform float f2;  // low \
uniform float f3;  // high \
 \
void main() { \
   vec4 value = texture2DRect(i1, gl_TexCoord[0].st); \
   if(value.r > f1) \
      gl_FragColor.r = f3; \
   else \
      gl_FragColor.r = f2; \
 \
   if(value.g > f1) \
      gl_FragColor.g = f3; \
   else \
      gl_FragColor.g = f2; \
 \
   if(value.b > f1) \
      gl_FragColor.b = f3; \
   else \
      gl_FragColor.b = f2; \
 \
   if(value.a > f1) \
      gl_FragColor.a = f3; \
   else \
      gl_FragColor.a = f2; \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/abs-1i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	out float4 color : COLOR )  \
{ \
   color = abs(texRECT(i1, texCoord0));  \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/abs-1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   gl_FragColor.rgba = abs(texture2DRect(i1, gl_TexCoord[0].st)); \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/abs_1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image1; \
 \
void main() { \
   gl_FragColor.rgba = abs(texture2DRect(image1, gl_TexCoord[0].st)); \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/acos-1i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	out float4 color : COLOR )  \
{ \
   color = acos(texRECT(i1, texCoord0));  \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/acos-1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   gl_FragColor.rgba = acos(texture2DRect(i1, gl_TexCoord[0].st)); \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/acos_1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image1; \
 \
void main() { \
   gl_FragColor.rgba = acos(texture2DRect(image1, gl_TexCoord[0].st)); \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/acosh-1i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	out float4 color : COLOR )  \
{ \
   float4 value = texRECT(i1, texCoord0); \
   color = log(value + sqrt(pow(value, 2.0) - 1.0)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/acosh-1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   vec4 value = texture2DRect(i1, gl_TexCoord[0].st); \
   gl_FragColor.rgba = log2(value + sqrt(pow(value, vec4(2.0)) - 1.0)) / 1.442695041; \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/asin-1i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	out float4 color : COLOR )  \
{ \
   color = asin(texRECT(i1, texCoord0));  \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/asin-1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   gl_FragColor.rgba = asin(texture2DRect(i1, gl_TexCoord[0].st)); \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/asin_1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image1; \
 \
void main() { \
   gl_FragColor.rgba = asin(texture2DRect(image1, gl_TexCoord[0].st)); \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/asinh-1i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	out float4 color : COLOR )  \
{ \
   float4 value = texRECT(i1, texCoord0); \
   color = log(value + sqrt(pow(value, 2.0) + 1.0)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/asinh-1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   vec4 value = texture2DRect(i1, gl_TexCoord[0].st); \
   gl_FragColor.rgba = log2(value + sqrt(pow(value, vec4(2.0)) + 1.0)) / 1.442695041; \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/atan-1i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	out float4 color : COLOR )  \
{ \
   color = atan(texRECT(i1, texCoord0));  \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/atan-1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   gl_FragColor.rgba = atan(texture2DRect(i1, gl_TexCoord[0].st)); \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/atan2-FI-1i1f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform float f1, \
	out float4 color : COLOR )  \
{ \
   color = atan2(float4(f1), texRECT(i1, texCoord0)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/atan2-FI-1i1f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
uniform float f1; \
 \
void main() { \
   float pi = 3.141592654; \
   float half_pi = 3.141592654 / 2.0; \
 \
   vec4 y = vec4(f1); \
   vec4 x = texture2DRect(i1, gl_TexCoord[0].st); \
   vec4 result; \
	 vec4 atan_value = atan(y/x); \
// channel-r \
   if(x.r > 0.0) { \
      result.r = atan_value.r; \
   } \
   else if(x.r < 0.0) { \
      result.r = atan_value.r + sign(y.r) * pi; \
   }	 \
   else if(y.r == 0.0) { \
      result.r = atan_value.r;  \
   }	 \
   else { \
      result.r = half_pi * sign(y.r);	 \
   } \
// channel-g \
   if(x.g > 0.0) { \
      result.g = atan_value.g; \
   } \
   else if(x.g < 0.0) { \
      result.g = atan_value.g + sign(y.g) * pi; \
   }	 \
   else if(y.g == 0.0) { \
      result.g = atan_value.g;  \
   }	 \
   else { \
      result.g = half_pi * sign(y.g);	 \
   } \
// channel-b \
   if(x.b > 0.0) { \
      result.b = atan_value.b; \
   } \
   else if(x.b < 0.0) { \
      result.b = atan_value.b + sign(y.b) * pi; \
   }	 \
   else if(y.b == 0.0) { \
      result.b = atan_value.b;  \
   }	 \
   else { \
      result.b = half_pi * sign(y.b);	 \
   } \
// channel-a \
   if(x.a > 0.0) { \
      result.a = atan_value.a; \
   } \
   else if(x.a < 0.0) { \
      result.a = atan_value.a + sign(y.a) * pi; \
   }	 \
   else if(y.a == 0.0) { \
      result.a = atan_value.a;  \
   }	 \
   else { \
      result.a = half_pi * sign(y.a);	 \
   } \
 \
   gl_FragColor.rgba = result; \
} \
  \
 \
float sign(float n) { \
   if(n < 0.0) { \
      return -1.0; \
   } \
   else { \
      return 1.0; \
   } \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/atan2-IF-1i1f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform float f1, \
	out float4 color : COLOR )  \
{ \
   color = atan2(texRECT(i1, texCoord0), f1); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/atan2-IF-1i1f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
uniform float f1; \
 \
void main() { \
   float pi = 3.141592654; \
   float half_pi = 3.141592654 / 2.0; \
 \
   vec4 y = texture2DRect(i1, gl_TexCoord[0].st); \
   vec4 x = vec4(f1);  \
   vec4 result; \
	 vec4 atan_value = atan(y/x); \
// channel-r \
   if(x.r > 0.0) { \
      result.r = atan_value.r; \
   } \
   else if(x.r < 0.0) { \
      result.r = atan_value.r + sign(y.r) * pi; \
   }	 \
   else if(y.r == 0.0) { \
      result.r = atan_value.r;  \
   }	 \
   else { \
      result.r = half_pi * sign(y.r);	 \
   } \
// channel-g \
   if(x.g > 0.0) { \
      result.g = atan_value.g; \
   } \
   else if(x.g < 0.0) { \
      result.g = atan_value.g + sign(y.g) * pi; \
   }	 \
   else if(y.g == 0.0) { \
      result.g = atan_value.g;  \
   }	 \
   else { \
      result.g = half_pi * sign(y.g);	 \
   } \
// channel-b \
   if(x.b > 0.0) { \
      result.b = atan_value.b; \
   } \
   else if(x.b < 0.0) { \
      result.b = atan_value.b + sign(y.b) * pi; \
   }	 \
   else if(y.b == 0.0) { \
      result.b = atan_value.b;  \
   }	 \
   else { \
      result.b = half_pi * sign(y.b);	 \
   } \
// channel-a \
   if(x.a > 0.0) { \
      result.a = atan_value.a; \
   } \
   else if(x.a < 0.0) { \
      result.a = atan_value.a + sign(y.a) * pi; \
   }	 \
   else if(y.a == 0.0) { \
      result.a = atan_value.a;  \
   }	 \
   else { \
      result.a = half_pi * sign(y.a);	 \
   } \
 \
   gl_FragColor.rgba = result; \
} \
  \
 \
float sign(float n) { \
   if(n < 0.0) { \
      return -1.0; \
   } \
   else { \
      return 1.0; \
   } \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/atan2-II-2i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform samplerRECT i2, \
	out float4 color : COLOR )  \
{ \
   color = atan2(texRECT(i1, texCoord0), texRECT(i2, texCoord0)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/atan2-II-2i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
uniform sampler2DRect i2; \
 \
void main() { \
   float pi = 3.141592654; \
   float half_pi = 3.141592654 / 2.0; \
 \
   vec4 y = texture2DRect(i1, gl_TexCoord[0].st); \
   vec4 x = texture2DRect(i2, gl_TexCoord[0].st); \
   vec4 result; \
	 vec4 atan_value = atan(y/x); \
// channel-r \
   if(x.r > 0.0) { \
      result.r = atan_value.r; \
   } \
   else if(x.r < 0.0) { \
      result.r = atan_value.r + sign(y.r) * pi; \
   }	 \
   else if(y.r == 0.0) { \
      result.r = atan_value.r;  \
   }	 \
   else { \
      result.r = half_pi * sign(y.r);	 \
   } \
// channel-g \
   if(x.g > 0.0) { \
      result.g = atan_value.g; \
   } \
   else if(x.g < 0.0) { \
      result.g = atan_value.g + sign(y.g) * pi; \
   }	 \
   else if(y.g == 0.0) { \
      result.g = atan_value.g;  \
   }	 \
   else { \
      result.g = half_pi * sign(y.g);	 \
   } \
// channel-b \
   if(x.b > 0.0) { \
      result.b = atan_value.b; \
   } \
   else if(x.b < 0.0) { \
      result.b = atan_value.b + sign(y.b) * pi; \
   }	 \
   else if(y.b == 0.0) { \
      result.b = atan_value.b;  \
   }	 \
   else { \
      result.b = half_pi * sign(y.b);	 \
   } \
// channel-a \
   if(x.a > 0.0) { \
      result.a = atan_value.a; \
   } \
   else if(x.a < 0.0) { \
      result.a = atan_value.a + sign(y.a) * pi; \
   }	 \
   else if(y.a == 0.0) { \
      result.a = atan_value.a;  \
   }	 \
   else { \
      result.a = half_pi * sign(y.a);	 \
   } \
 \
   gl_FragColor.rgba = result; \
} \
  \
 \
float sign(float n) { \
   if(n < 0.0) { \
      return -1.0; \
   } \
   else { \
      return 1.0; \
   } \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/atanh-1i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	out float4 color : COLOR )  \
{ \
   float4 value = texRECT(i1, texCoord0); \
   color = 0.5 * log((value + 1.0) / (value - 1.0)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/atanh-1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   vec4 value = texture2DRect(i1, gl_TexCoord[0].st); \
   gl_FragColor.rgba = 0.5 * (log2((value + 1.0) / (value - 1.0)) / 1.442695041); \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/cbrt-1i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	out float4 color : COLOR )  \
{ \
   color = pow(texRECT(i1, texCoord0), 1.0/3.0); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/cbrt-1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   gl_FragColor.rgba = pow(texture2DRect(i1, gl_TexCoord[0].st), vec4(1.0/3.0)); \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/cbrt_1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image1; \
 \
void main() { \
   gl_FragColor.rgba = pow(texture2DRect(image1, gl_TexCoord[0].st), 1.0/3.0); \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/ceil-1i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	out float4 color : COLOR )  \
{ \
   color = ceil(texRECT(i1, texCoord0)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/ceil-1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   gl_FragColor.rgba = ceil(texture2DRect(i1, gl_TexCoord[0].st)); \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/ceil_1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image1; \
 \
void main() { \
   gl_FragColor.rgba = ceil(texture2DRect(image1, gl_TexCoord[0].st)); \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/copysign-FI-1i1f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
uniform float f1; \
 \
void main() { \
   vec4 value = vec4(f1); \
   vec4 ref = texture2DRect(i1, gl_TexCoord[0].st); \
 \
   if(value.r * ref.r < 0.0) \
      value.r *= -1.0; \
   if(value.g * ref.g < 0.0) \
      value.g *= -1.0; \
   if(value.b * ref.b < 0.0) \
      value.b *= -1.0; \
   if(value.a * ref.a < 0.0) \
      value.a *= -1.0; \
 \
   gl_FragColor.rgba = value; \
} \
  \
 \
";

standard_shaders_map["./StandardShaders/ImageMath/copysign-IF-1i1f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform float f1, \
	out float4 color : COLOR )  \
{ \
   float4 value = texRECT(i1, texCoord0); \
   float4 ref = float4(f1); \
 \
   if(value.r * ref.r < 0.0) { \
      value.r *= -1.0; \
   }  \
   if(value.g * ref.g < 0.0) { \
      value.g *= -1.0; \
   }  \
   if(value.b * ref.b < 0.0) { \
      value.b *= -1.0; \
   }  \
   if(value.a * ref.a < 0.0) { \
      value.a *= -1.0; \
   }  \
 \
   color = value; \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/copysign-IF-1i1f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
uniform float f1; \
 \
void main() { \
   vec4 value = texture2DRect(i1, gl_TexCoord[0].st); \
   vec4 ref = vec4(f1); \
 \
   if(value.r * ref.r < 0.0) \
      value.r *= -1.0; \
   if(value.g * ref.g < 0.0) \
      value.g *= -1.0; \
   if(value.b * ref.b < 0.0) \
      value.b *= -1.0; \
   if(value.a * ref.a < 0.0) \
      value.a *= -1.0; \
 \
   gl_FragColor.rgba = value; \
} \
  \
 \
";

standard_shaders_map["./StandardShaders/ImageMath/copysign-II-2i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform samplerRECT i2, \
	out float4 color : COLOR )  \
{ \
   float4 value = texRECT(i1, texCoord0); \
   float4 ref = texRECT(i2, texCoord0); \
 \
   if(value.r * ref.r < 0.0) { \
      value.r *= -1.0; \
   }  \
   if(value.g * ref.g < 0.0) { \
      value.g *= -1.0; \
   }  \
   if(value.b * ref.b < 0.0) { \
      value.b *= -1.0; \
   }  \
   if(value.a * ref.a < 0.0) { \
      value.a *= -1.0; \
   }  \
 \
   color = value; \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/copysign-II-2i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
uniform sampler2DRect i2; \
 \
void main() { \
   vec4 value = texture2DRect(i1, gl_TexCoord[0].st); \
   vec4 ref = texture2DRect(i2, gl_TexCoord[0].st); \
 \
   if(value.r * ref.r < 0.0) \
      value.r *= -1.0; \
   if(value.g * ref.g < 0.0) \
      value.g *= -1.0; \
   if(value.b * ref.b < 0.0) \
      value.b *= -1.0; \
   if(value.a * ref.a < 0.0) \
      value.a *= -1.0; \
 \
   gl_FragColor.rgba = value; \
} \
  \
 \
";

standard_shaders_map["./StandardShaders/ImageMath/cos-1i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	out float4 color : COLOR )  \
{ \
   color = cos(texRECT(i1, texCoord0)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/cos-1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   gl_FragColor.rgba = cos(texture2DRect(i1, gl_TexCoord[0].st)); \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/cos_1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image1; \
 \
void main() { \
   gl_FragColor.rgba = cos(texture2DRect(image1, gl_TexCoord[0].st)); \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/cosh-1i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	out float4 color : COLOR )  \
{ \
   color = cosh(texRECT(i1, texCoord0)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/cosh-1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   vec4 value = texture2DRect(i1, gl_TexCoord[0].st); \
   vec4 e = vec4(2.7182818284590451); \
   gl_FragColor.rgba = (pow(e, value) + pow(e, -1.0 * value)) / 2.0; \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/exp-1i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	out float4 color : COLOR )  \
{ \
   color = exp(texRECT(i1, texCoord0)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/exp-1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   gl_FragColor.rgba = pow(vec4(2.7182818284590451), texture2DRect(i1, gl_TexCoord[0].st)); \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/exp2-1i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	out float4 color : COLOR )  \
{ \
   color = exp2(texRECT(i1, texCoord0)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/exp2-1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   gl_FragColor.rgba = exp2(texture2DRect(i1, gl_TexCoord[0].st)); \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/exp2_1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image1; \
 \
void main() { \
   gl_FragColor.rgba = exp2(texture2DRect(image1, gl_TexCoord[0].st)); \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/exp_1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image1; \
 \
void main() { \
   gl_FragColor.rgba = pow(2.7182818284590451, texture2DRect(image1, gl_TexCoord[0].st)); \
   gl_FragColor.rgba = texture2DRect(image1, gl_TexCoord[0].st); \
 \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/expm1-1i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	out float4 color : COLOR )  \
{ \
   color = exp(texRECT(i1, texCoord0)) - float4(1.0); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/expm1-1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   gl_FragColor.rgba = pow(vec4(2.7182818284590451), texture2DRect(i1, gl_TexCoord[0].st)) - 1.0; \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/expm1_1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image1; \
 \
void main() { \
   gl_FragColor.rgba = pow(2.7182818284590451, texture2DRect(image1, gl_TexCoord[0].st)) - 1.0; \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/fdim-2i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
uniform sampler2DRect i2; \
 \
void main() { \
   gl_FragColor.rgba = max(0.0, texture2DRect(i1, gl_TexCoord[0].st) - texture2DRect(i1, gl_TexCoord[0].st)); \
} \
  \
 \
";

standard_shaders_map["./StandardShaders/ImageMath/fdim-FI-1i1f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform float f1, \
	out float4 color : COLOR )  \
{ \
   color = max(float4(0.0), float4(f1) - texRECT(i1, texCoord0)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/fdim-FI-1i1f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
uniform float f1; \
 \
void main() { \
   gl_FragColor.rgba = max(vec4(0.0),  vec4(f1) - texture2DRect(i1, gl_TexCoord[0].st)); \
} \
  \
 \
";

standard_shaders_map["./StandardShaders/ImageMath/fdim-IF-1i1f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform float f1, \
	out float4 color : COLOR )  \
{ \
   color = max(float4(0.0), texRECT(i1, texCoord0) - float4(f1)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/fdim-IF-1i1f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
uniform float f1; \
 \
void main() { \
   gl_FragColor.rgba = max(vec4(0.0), texture2DRect(i1, gl_TexCoord[0].st) - vec4(f1)); \
} \
  \
 \
";

standard_shaders_map["./StandardShaders/ImageMath/fdim-II-2i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform samplerRECT i2, \
	out float4 color : COLOR )  \
{ \
   color = max(float4(0.0), texRECT(i1, texCoord0) - texRECT(i2, texCoord0)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/fdim-II-2i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
uniform sampler2DRect i2; \
 \
void main() { \
   gl_FragColor.rgba = max(vec4(0.0), texture2DRect(i1, gl_TexCoord[0].st) - texture2DRect(i2, gl_TexCoord[0].st)); \
} \
  \
 \
";

standard_shaders_map["./StandardShaders/ImageMath/floor-1i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	out float4 color : COLOR )  \
{ \
   color = floor(texRECT(i1, texCoord0)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/floor-1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   gl_FragColor.rgba = floor(texture2DRect(i1, gl_TexCoord[0].st)); \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/floor_1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image1; \
 \
void main() { \
   gl_FragColor.rgba = floor(texture2DRect(image1, gl_TexCoord[0].st)); \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/FUNC_1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image1; \
 \
void main() { \
   gl_FragColor.rgba = FUNC(texture2DRect(image1, gl_TexCoord[0].st)); \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/hypot-FI-1i1f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform float f1, \
	out float4 color : COLOR )  \
{ \
   color = sqrt(pow(texRECT(i1, texCoord0), float4(2.0)) + float4(pow(f1, 2.0))); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/hypot-FI-1i1f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
uniform float f1; \
 \
void main() { \
   gl_FragColor.rgba = sqrt(pow(texture2DRect(i1, gl_TexCoord[0].st), 2.0) + pow(vec4(f1), 2.0)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/hypot-IF-1i1f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform float f1, \
	out float4 color : COLOR )  \
{ \
   color = sqrt(float4(pow(f1, 2.0)) + pow(texRECT(i1, texCoord0), float4(2.0))); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/hypot-IF-1i1f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
uniform float f1; \
 \
void main() { \
   gl_FragColor.rgba = sqrt(pow(texture2DRect(i1, gl_TexCoord[0].st), vec4(2.0)) + pow(vec4(f1), vec4(2.0))); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/hypot-II-2i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform samplerRECT i2, \
	out float4 color : COLOR )  \
{ \
   color = sqrt(pow(texRECT(i1, texCoord0), float4(2.0)) + pow(texRECT(i2, texCoord0), float4(2.0))); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/hypot-II-2i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
uniform sampler2DRect i2; \
 \
void main() { \
   gl_FragColor.rgba = sqrt(pow(texture2DRect(i1, gl_TexCoord[0].st), 2.0) + pow(texture2DRect(i2, gl_TexCoord[0].st), 2.0)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/log-1i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	out float4 color : COLOR )  \
{ \
   color = log(texRECT(i1, texCoord0)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/log-1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   gl_FragColor.rgba = log2(texture2DRect(i1, gl_TexCoord[0].st)) / 1.442695041; \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/log10-1i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	out float4 color : COLOR )  \
{ \
   color = log10(texRECT(i1, texCoord0)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/log10-1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   gl_FragColor.rgba = log2(texture2DRect(i1, gl_TexCoord[0].st)) / 3.321928095; \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/log10_1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image1; \
 \
void main() { \
   gl_FragColor.rgba = log2(texture2DRect(image1, gl_TexCoord[0].st)) / 3.321928095; \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/log1p-1i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	out float4 color : COLOR )  \
{ \
   color = log(texRECT(i1, texCoord0)) + 1.0; \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/log1p-1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   gl_FragColor.rgba = (log2(texture2DRect(i1, gl_TexCoord[0].st)) / 1.442695041) + 1.0; \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/log1p_1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image1; \
 \
void main() { \
   gl_FragColor.rgba = (log2(texture2DRect(image1, gl_TexCoord[0].st)) / 1.442695041) + 1.0; \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/log_1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image1; \
 \
void main() { \
   gl_FragColor.rgba = log2(texture2DRect(image1, gl_TexCoord[0].st)) / 1.442695041; \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/pow-FI-1i1f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform float f1, \
	out float4 color : COLOR )  \
{ \
   color = pow(float4(f1), texRECT(i1, texCoord0)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/pow-FI-1i1f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
uniform float f1; \
 \
void main() { \
   gl_FragColor.rgba = pow(vec4(f1), texture2DRect(i1, gl_TexCoord[0].st)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/pow-IF-1i1f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform float f1, \
	out float4 color : COLOR )  \
{ \
   color = pow(texRECT(i1, texCoord0), float4(f1)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/pow-IF-1i1f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
uniform float f1; \
 \
void main() { \
   gl_FragColor.rgba = pow(texture2DRect(i1, gl_TexCoord[0].st), vec4(f1)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/pow-II-2i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform samplerRECT i2, \
	out float4 color : COLOR )  \
{ \
   color = pow(texRECT(i1, texCoord0), texRECT(i2, texCoord0)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/pow-II-2i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image1; \
uniform sampler2DRect image2; \
 \
 \
void main() { \
   gl_FragColor.rgba = pow(texture2DRect(image1, gl_TexCoord[0].st), texture2DRect(image2, gl_TexCoord[0].st)); \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/pow-IS-1i1f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image1; \
uniform float float1; \
 \
void main() { \
   gl_FragColor.rgba = pow(texture2DRect(image1, gl_TexCoord[0].st), vec4(float1)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/pow-SI-1i1f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image1; \
uniform float float1; \
 \
void main() { \
   gl_FragColor.rgba = pow(vec4(float1), texture2DRect(image1, gl_TexCoord[0].st)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/round-1i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	out float4 color : COLOR )  \
{ \
   color = round(texRECT(i1, texCoord0)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/round-1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   vec4 value = texture2DRect(i1, gl_TexCoord[0].st); \
   vec4 rounded = floor(value); \
   if(value.r - rounded.r > 0.5) { \
      rounded.r += 1.0; \
   } \
   if(value.g - rounded.g > 0.5) {  \
      rounded.g += 1.0; \
   } \
   if(value.b - rounded.b > 0.5) { \
      rounded.b += 1.0; \
   } \
   if(value.a - rounded.a > 0.5) { \
      rounded.b += 1.0; \
   } \
   gl_FragColor.rgba = rounded; \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/sin-1i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	out float4 color : COLOR )  \
{ \
   color = sin(texRECT(i1, texCoord0)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/sin-1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   gl_FragColor.rgba = sin(texture2DRect(i1, gl_TexCoord[0].st)); \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/sin_1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image1; \
 \
void main() { \
   gl_FragColor.rgba = sin(texture2DRect(image1, gl_TexCoord[0].st)); \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/sinh-1i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	out float4 color : COLOR )  \
{ \
   color = sinh(texRECT(i1, texCoord0)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/sinh-1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   vec4 value = texture2DRect(i1, gl_TexCoord[0].st); \
   vec4 e = vec4(2.7182818284590451); \
   gl_FragColor.rgba = (pow(e, value) + pow(e, -1.0 * value)) / 2.0; \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/sqrt-1i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	out float4 color : COLOR )  \
{ \
   color = sqrt(texRECT(i1, texCoord0)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/sqrt-1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   gl_FragColor.rgba = sqrt(texture2DRect(i1, gl_TexCoord[0].st)); \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/sqrt_1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image1; \
 \
void main() { \
   gl_FragColor.rgba = sqrt(texture2DRect(image1, gl_TexCoord[0].st)); \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/tan-1i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	out float4 color : COLOR )  \
{ \
   color = tan(texRECT(i1, texCoord0)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/tan-1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   gl_FragColor.rgba = tan(texture2DRect(i1, gl_TexCoord[0].st)); \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/tan_1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image1; \
 \
void main() { \
   gl_FragColor.rgba = tan(texture2DRect(image1, gl_TexCoord[0].st)); \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/tanh-1i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	out float4 color : COLOR )  \
{ \
   color = tanh(texRECT(i1, texCoord0)); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/tanh-1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   vec4 value = texture2DRect(i1, gl_TexCoord[0].st); \
   vec4 e = vec4(2.7182818284590451); \
   vec4 e2x = pow(e, 2.0 * value); \
   gl_FragColor.rgba = (e2x + 1.0) / (e2x - 1.0); \
} \
  \
";

standard_shaders_map["./StandardShaders/ImageMath/add-IF-1i1f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform float f1, \
	out float4 color : COLOR )  \
{ \
   color = texRECT(i1, texCoord0) + f1; \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/add-IF-1i1f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
uniform float f1; \
 \
void main() { \
   gl_FragColor.rgba = texture2DRect(i1, gl_TexCoord[0].st) + f1; \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/add-II-2i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform samplerRECT i2, \
	out float4 color : COLOR )  \
{ \
   color = texRECT(i1, texCoord0) + texRECT(i2, texCoord0); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/add-II-2i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
uniform sampler2DRect i2; \
 \
void main() { \
   gl_FragColor.rgba = texture2DRect(i1, gl_TexCoord[0].st) +  texture2DRect(i2, gl_TexCoord[0].st); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/divide-FI-1i1f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform float f1, \
	out float4 color : COLOR )  \
{ \
   color = f1 / texRECT(i1, texCoord0); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/divide-FI-1i1f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
uniform float f1; \
 \
void main() { \
   gl_FragColor.rgba = f1 / texture2DRect(i1, gl_TexCoord[0].st); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/divide-IF-1i1f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform float f1, \
	out float4 color : COLOR )  \
{ \
   color = texRECT(i1, texCoord0) / f1; \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/divide-IF-1i1f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
uniform float f1; \
 \
void main() { \
   gl_FragColor.rgba = texture2DRect(i1, gl_TexCoord[0].st) / f1; \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/divide-II-2i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform samplerRECT i2, \
	out float4 color : COLOR )  \
{ \
   color = texRECT(i1, texCoord0) / texRECT(i2, texCoord0); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/divide-II-2i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
uniform sampler2DRect i2; \
 \
void main() { \
   gl_FragColor.rgba = texture2DRect(i1, gl_TexCoord[0].st) / texture2DRect(i2, gl_TexCoord[0].st); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/multiply-IF-1i1f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform float f1, \
	out float4 color : COLOR )  \
{ \
   color = texRECT(i1, texCoord0) * f1; \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/multiply-IF-1i1f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
uniform float f1; \
 \
void main() { \
   gl_FragColor.rgba = texture2DRect(i1, gl_TexCoord[0].st) * f1; \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/multiply-II-2i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform samplerRECT i2, \
	out float4 color : COLOR )  \
{ \
   color = texRECT(i1, texCoord0) * texRECT(i2, texCoord0); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/multiply-II-2i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
uniform sampler2DRect i2; \
 \
void main() { \
   gl_FragColor.rgba = texture2DRect(i1, gl_TexCoord[0].st) * texture2DRect(i2, gl_TexCoord[0].st); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/subtract-FI-1i1f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform float f1, \
	out float4 color : COLOR )  \
{ \
   color = f1 - texRECT(i1, texCoord0); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/subtract-FI-1i1f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
uniform float f1; \
 \
void main() { \
   gl_FragColor.rgba = f1 - texture2DRect(i1, gl_TexCoord[0].st); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/subtract-IF-1i1f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform float f1, \
	out float4 color : COLOR )  \
{ \
   color = texRECT(i1, texCoord0) - f1; \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/subtract-IF-1i1f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
uniform float f1; \
 \
void main() { \
   gl_FragColor.rgba = texture2DRect(i1, gl_TexCoord[0].st) - f1; \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/subtract-II-2i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform samplerRECT i2, \
	out float4 color : COLOR )  \
{ \
   color = texRECT(i1, texCoord0) - texRECT(i2, texCoord0); \
} \
";

standard_shaders_map["./StandardShaders/ImageMath/subtract-II-2i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
uniform sampler2DRect i2; \
 \
void main() { \
   gl_FragColor.rgba = texture2DRect(i1, gl_TexCoord[0].st) - texture2DRect(i2, gl_TexCoord[0].st); \
} \
";

standard_shaders_map["./StandardShaders/Statisticsmax-channels_cg_frag_r"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT image, \
	out float4 color : COLOR )  \
{ \
   float4 cPixel = texRECT(image, texCoord0); \
   color = cPixel.r; \
} \
";

standard_shaders_map["./StandardShaders/Statisticsmax-channels_cg_frag_rgb"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT image, \
	out float4 color : COLOR )  \
{ \
   float4 cPixel = texRECT(image, texCoord0); \
   float cMax = cPixel.r; \
   cMax = max(cMax, cPixel.g); \
   cMax = max(cMax, cPixel.b); \
 \
   color = cMax; \
} \
";

standard_shaders_map["./StandardShaders/Statisticsmax-channels_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT image, \
	out float4 color : COLOR )  \
{ \
   float4 cPixel = texRECT(image, texCoord0); \
   float cMax = cPixel.r; \
   cMax = max(cMax, cPixel.g); \
   cMax = max(cMax, cPixel.b); \
   cMax = max(cMax, cPixel.a); \
 \
   color = cMax; \
} \
";

standard_shaders_map["./StandardShaders/Statisticsmax-channels_gl_frag_r"] = " \
 \
uniform sampler2DRect image; \
 \
void main() { \
   vec2 coord = gl_TexCoord[0].st; \
 \
   vec4 cPixel = texture2DRect(image, coord); \
   float cMax = cPixel.r; \
 \
   gl_FragColor.r = cMax; \
} \
";

standard_shaders_map["./StandardShaders/Statisticsmax-channels_gl_frag_rgb"] = " \
 \
uniform sampler2DRect image; \
 \
void main() { \
   vec2 coord = gl_TexCoord[0].st; \
 \
   vec4 cPixel = texture2DRect(image, coord); \
   float cMax = cPixel.r; \
   cMax = max(cMax, cPixel.g); \
   cMax = max(cMax, cPixel.b); \
 \
   gl_FragColor.r = cMax; \
} \
";

standard_shaders_map["./StandardShaders/Statisticsmax-channels_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image; \
 \
void main() { \
   vec2 coord = gl_TexCoord[0].st; \
 \
   vec4 cPixel = texture2DRect(image, coord); \
   float cMax = cPixel.r; \
   cMax = max(cMax, cPixel.g); \
   cMax = max(cMax, cPixel.b); \
   cMax = max(cMax, cPixel.a); \
 \
   gl_FragColor.r = cMax; \
} \
";

standard_shaders_map["./StandardShaders/Statisticsmax-quad_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT image, \
	out float4 color : COLOR )  \
{ \
   float4 cPixel = texRECT(image, texCoord0); \
   float cMax = cPixel.r; \
 \
   cPixel = texRECT(image, texCoord0 + float2(1.0, 0.0)); \
   cMax = max(cMax, cPixel.r); \
 \
   cPixel = texRECT(image, texCoord0 + float2(0.0, 1.0)); \
   cMax = max(cMax, cPixel.r); \
 \
   cPixel = texRECT(image, texCoord0 + float2(1.0, 1.0)); \
   cMax = max(cMax, cPixel.r); \
 \
   color = cMax; \
 \
} \
 \
";

standard_shaders_map["./StandardShaders/Statisticsmax-quad_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image; \
 \
void main() { \
   vec2 coord = gl_TexCoord[0].st - vec2(1.0, 1.0); \
 \
   vec4 cPixel = texture2DRect(image, coord); \
   float cMax = cPixel.r; \
 \
   cPixel = texture2DRect(image, coord + vec2(1.0, 0.0)); \
   cMax = max(cMax, cPixel.r); \
 \
   cPixel = texture2DRect(image, coord + vec2(0.0, 1.0)); \
   cMax = max(cMax, cPixel.r); \
 \
   cPixel = texture2DRect(image, coord + vec2(1.0, 1.0)); \
   cMax = max(cMax, cPixel.r); \
		 \
   gl_FragColor.r = cMax; \
} \
  \
";

standard_shaders_map["./StandardShaders/Statisticsmin-channels_cg_frag_r"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT image, \
	out float4 color : COLOR )  \
{ \
   float4 cPixel = texRECT(image, texCoord0); \
   float cMin = cPixel.r; \
   color = cMin; \
} \
";

standard_shaders_map["./StandardShaders/Statisticsmin-channels_cg_frag_rgb"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT image, \
	out float4 color : COLOR )  \
{ \
   float4 cPixel = texRECT(image, texCoord0); \
   float cMin = cPixel.r; \
   cMin = min(cMin, cPixel.g); \
   cMin = min(cMin, cPixel.b); \
 \
   color = cMin; \
} \
";

standard_shaders_map["./StandardShaders/Statisticsmin-channels_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT image, \
	out float4 color : COLOR )  \
{ \
   float4 cPixel = texRECT(image, texCoord0); \
   float cMin = cPixel.r; \
   cMin = min(cMin, cPixel.g); \
   cMin = min(cMin, cPixel.b); \
   cMin = min(cMin, cPixel.a); \
 \
   color = cMin; \
} \
";

standard_shaders_map["./StandardShaders/Statisticsmin-channels_gl_frag_r"] = " \
 \
uniform sampler2DRect image; \
 \
void main() { \
   vec2 coord = gl_TexCoord[0].st; \
   gl_FragColor.r = texture2DRect(image, coord).r; \
} \
";

standard_shaders_map["./StandardShaders/Statisticsmin-channels_gl_frag_rgb"] = " \
 \
uniform sampler2DRect image; \
 \
void main() { \
   vec2 coord = gl_TexCoord[0].st; \
 \
   vec4 cPixel = texture2DRect(image, coord); \
   float cMin = cPixel.r; \
   cMin = min(cMin, cPixel.g); \
   cMin = min(cMin, cPixel.b); \
   cMin = min(cMin, cPixel.a); \
 \
   gl_FragColor.r = cMin; \
} \
";

standard_shaders_map["./StandardShaders/Statisticsmin-channels_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image; \
 \
void main() { \
   vec2 coord = gl_TexCoord[0].st; \
 \
   vec4 cPixel = texture2DRect(image, coord); \
   float cMin = cPixel.r; \
   cMin = min(cMin, cPixel.g); \
   cMin = min(cMin, cPixel.b); \
   cMin = min(cMin, cPixel.a); \
 \
   gl_FragColor.r = cMin; \
} \
";

standard_shaders_map["./StandardShaders/Statisticsmin-quad_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT image, \
	out float4 color : COLOR )  \
{ \
   float4 cPixel = texRECT(image, texCoord0); \
   float cMin = cPixel.r; \
 \
   cPixel = texRECT(image, texCoord0 + float2(1.0, 0.0)); \
   cMin = min(cMin, cPixel.r); \
 \
   cPixel = texRECT(image, texCoord0 + float2(0.0, 1.0)); \
   cMin = min(cMin, cPixel.r); \
 \
   cPixel = texRECT(image, texCoord0 + float2(1.0, 1.0)); \
   cMin = min(cMin, cPixel.r); \
 \
   color = cMin; \
 \
} \
 \
";

standard_shaders_map["./StandardShaders/Statisticsmin-quad_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image; \
 \
void main() { \
   vec2 coord = gl_TexCoord[0].st - vec2(1.0, 1.0); \
 \
   vec4 cPixel = texture2DRect(image, coord); \
   float cMin = cPixel.r; \
 \
   cPixel = texture2DRect(image, coord + vec2(1.0, 0.0)); \
   cMin = min(cMin, cPixel.r); \
 \
   cPixel = texture2DRect(image, coord + vec2(0.0, 1.0)); \
   cMin = min(cMin, cPixel.r); \
 \
   cPixel = texture2DRect(image, coord + vec2(1.0, 1.0)); \
   cMin = min(cMin, cPixel.r); \
 \
   gl_FragColor.r = cMin; \
} \
  \
";

standard_shaders_map["./StandardShaders/Statisticssum-channels_cg_frag_r"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT image, \
	out float4 color : COLOR )  \
{ \
   float4 cPixel = texRECT(image, texCoord0); \
   color = cPixel.r; \
} \
";

standard_shaders_map["./StandardShaders/Statisticssum-channels_cg_frag_rgb"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT image, \
	out float4 color : COLOR )  \
{ \
   float4 cPixel = texRECT(image, texCoord0); \
   color = cPixel.r + cPixel.g + cPixel.b; \
} \
";

standard_shaders_map["./StandardShaders/Statisticssum-channels_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT image, \
	out float4 color : COLOR )  \
{ \
   float4 cPixel = texRECT(image, texCoord0); \
   color = cPixel.r + cPixel.g + cPixel.b + cPixel.a; \
} \
";

standard_shaders_map["./StandardShaders/Statisticssum-channels_gl_frag_r"] = " \
 \
uniform sampler2DRect image; \
 \
void main() { \
   vec2 coord = gl_TexCoord[0].st; \
 \
   vec4 cPixel = texture2DRect(image, coord); \
   float sum = cPixel.r; \
		 \
   gl_FragColor.r = sum; \
} \
";

standard_shaders_map["./StandardShaders/Statisticssum-channels_gl_frag_rgb"] = " \
 \
uniform sampler2DRect image; \
 \
void main() { \
   vec2 coord = gl_TexCoord[0].st; \
 \
   vec4 cPixel = texture2DRect(image, coord); \
   float sum = cPixel.r; \
   sum += cPixel.g; \
   sum += cPixel.b; \
 \
   gl_FragColor.r = sum; \
} \
";

standard_shaders_map["./StandardShaders/Statisticssum-channels_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image; \
 \
void main() { \
   vec2 coord = gl_TexCoord[0].st; \
 \
   vec4 cPixel = texture2DRect(image, coord); \
   float sum = cPixel.r; \
   sum += cPixel.g; \
   sum += cPixel.b; \
   sum += cPixel.a; \
		 \
   gl_FragColor.r = sum; \
} \
";

standard_shaders_map["./StandardShaders/Statisticssum-quad_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT image, \
	out float4 color : COLOR )  \
{ \
   float cSum = 0; \
   cSum += texRECT(image, texCoord0).r; \
   cSum += texRECT(image, texCoord0 + float2(1.0, 0.0)).r; \
   cSum += texRECT(image, texCoord0 + float2(0.0, 1.0)).r; \
   cSum += texRECT(image, texCoord0 + float2(1.0, 1.0)).r; \
   color = cSum; \
 \
} \
 \
";

standard_shaders_map["./StandardShaders/Statisticssum-quad_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image; \
 \
void main() { \
   vec2 coord = gl_TexCoord[0].st - vec2(1.0, 1.0); \
 \
   vec4 cPixel = texture2DRect(image, coord); \
   float sum = cPixel.r; \
 \
   cPixel = texture2DRect(image, coord + vec2(1.0, 0.0)); \
   sum += cPixel.r; \
 \
   cPixel = texture2DRect(image, coord + vec2(0.0, 1.0)); \
   sum += cPixel.r; \
 \
   cPixel = texture2DRect(image, coord + vec2(1.0, 1.0)); \
   sum += cPixel.r; \
		 \
   gl_FragColor.r = sum; \
} \
";

standard_shaders_map["./StandardShaders/Interpolation/interpolation-bicubic_cg_frag_rgba"] = " \
  \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT image, \
	out float4 color : COLOR )  \
{ \
   float2 ij = texCoord0; \
   float2 xy = floor(ij); \
   float2 normxy = ij - xy; \
   float2 st0 = ((2.0 - normxy) * normxy - 1.0) * normxy; \
   float2 st1 = (3.0 * normxy - 5.0) * normxy * normxy + 2.0; \
   float2 st2 = ((4.0 - 3.0 * normxy) * normxy + 1.0) * normxy; \
   float2 st3 = (normxy - 1.0) * normxy * normxy; \
 \
   float4 row0 =  \
	st0.s * texRECT(image, xy + float2(-1.0, -1.0)) +  \
	st1.s * texRECT(image, xy + float2(0.0, -1.0)) +  \
	st2.s * texRECT(image, xy + float2(1.0, -1.0)) +  \
	st3.s * texRECT(image, xy + float2(2.0, -1.0)); \
 \
   float4 row1 =  \
	st0.s * texRECT(image, xy + float2(-1.0, 0.0)) +  \
	st1.s * texRECT(image, xy + float2(0.0, 0.0)) +  \
	st2.s * texRECT(image, xy + float2(1.0, 0.0)) +  \
	st3.s * texRECT(image, xy + float2(2.0, 0.0)); \
 \
   float4 row2 =  \
	st0.s * texRECT(image, xy + float2(-1.0, 1.0)) +  \
	st1.s * texRECT(image, xy + float2(0.0, 1.0)) +  \
	st2.s * texRECT(image, xy + float2(1.0, 1.0)) +  \
	st3.s * texRECT(image, xy + float2(2.0, 1.0)); \
 \
   float4 row3 =  \
	st0.s * texRECT(image, xy + float2(-1.0, 2.0)) +  \
	st1.s * texRECT(image, xy + float2(0.0, 2.0)) +  \
	st2.s * texRECT(image, xy + float2(1.0, 2.0)) +  \
	st3.s * texRECT(image, xy + float2(2.0, 2.0)); \
 \
   color = 0.25 * ((st0.t * row0) + (st1.t * row1) + (st2.t * row2) + (st3.t * row3));  \
 \
} \
";

standard_shaders_map["./StandardShaders/Interpolation/interpolation-bicubic_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image; \
 \
void main() { \
   vec2 ij = gl_TexCoord[0].st; \
   vec2	xy = floor(ij); \
   vec2 normxy = ij - xy; \
   vec2 st0 = ((2.0 - normxy) * normxy - 1.0) * normxy; \
   vec2 st1 = (3.0 * normxy - 5.0) * normxy * normxy + 2.0; \
   vec2 st2 = ((4.0 - 3.0 * normxy) * normxy + 1.0) * normxy; \
   vec2 st3 = (normxy - 1.0) * normxy * normxy; \
 \
   vec4 row0 =  \
	st0.s * texture2DRect(image, xy + vec2(-1.0, -1.0)) +  \
	st1.s * texture2DRect(image, xy + vec2(0.0, -1.0)) +  \
	st2.s * texture2DRect(image, xy + vec2(1.0, -1.0)) +  \
	st3.s * texture2DRect(image, xy + vec2(2.0, -1.0)); \
 \
   vec4 row1 =  \
	st0.s * texture2DRect(image, xy + vec2(-1.0, 0.0)) +  \
	st1.s * texture2DRect(image, xy + vec2(0.0, 0.0)) +  \
	st2.s * texture2DRect(image, xy + vec2(1.0, 0.0)) +  \
	st3.s * texture2DRect(image, xy + vec2(2.0, 0.0)); \
 \
   vec4 row2 =  \
	st0.s * texture2DRect(image, xy + vec2(-1.0, 1.0)) +  \
	st1.s * texture2DRect(image, xy + vec2(0.0, 1.0)) +  \
	st2.s * texture2DRect(image, xy + vec2(1.0, 1.0)) +  \
	st3.s * texture2DRect(image, xy + vec2(2.0, 1.0)); \
 \
   vec4 row3 =  \
	st0.s * texture2DRect(image, xy + vec2(-1.0, 2.0)) +  \
	st1.s * texture2DRect(image, xy + vec2(0.0, 2.0)) +  \
	st2.s * texture2DRect(image, xy + vec2(1.0, 2.0)) +  \
	st3.s * texture2DRect(image, xy + vec2(2.0, 2.0)); \
 \
   gl_FragColor.rgba = 0.25 * ((st0.t * row0) + (st1.t * row1) + (st2.t * row2) + (st3.t * row3));  \
} \
";

standard_shaders_map["./StandardShaders/Interpolation/interpolation-bilinear_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT image, \
	out float4 color : COLOR )  \
{ \
 \
   float2 coord_fract = frac(texCoord0); \
 \
   color =  \
	texRECT(image, texCoord0) * (1.0 - coord_fract.x) * (1.0 - coord_fract.y) + \
	texRECT(image, texCoord0 + float2(1.0, 0.0)) * coord_fract.x * (1.0 - coord_fract.y) + \
	texRECT(image, texCoord0 + float2(0.0, 1.0)) * (1.0 - coord_fract.x) * coord_fract.y + \
	texRECT(image, texCoord0 + float2(1.0, 1.0)) * coord_fract.x * coord_fract.y; \
} \
";

standard_shaders_map["./StandardShaders/Interpolation/interpolation-bilinear_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image; \
 \
void main() { \
   vec2 coord = gl_TexCoord[0].st; \
   vec2 coord_fract = fract(coord); \
   gl_FragColor.rgba =  \
	texture2DRect(image, coord) * (1.0 - coord_fract.s) * (1.0 - coord_fract.t) + \
	texture2DRect(image, coord + vec2(1.0, 0.0)) * coord_fract.s * (1.0 - coord_fract.t) + \
	texture2DRect(image, coord + vec2(0.0, 1.0)) * (1.0 - coord_fract.s) * coord_fract.t + \
	texture2DRect(image, coord + vec2(1.0, 1.0)) * coord_fract.s * coord_fract.t; \
} \
";

standard_shaders_map["./StandardShaders/Interpolation/interpolation-nearest-pixel_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT image, \
	out float4 color : COLOR )  \
{ \
   color = abs(texRECT(image, round(texCoord0)));  \
} \
";

standard_shaders_map["./StandardShaders/Interpolation/interpolation-nearest-pixel_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image; \
 \
void main() { \
   vec2 coord = gl_TexCoord[0].st; \
   vec2 coord_round = floor(coord); \
   gl_FragColor.rgba = texture2DRect(image, coord_round); \
} \
";

standard_shaders_map["./StandardShaders/Manipulation/PackGrayIntoRGBA_gl_frag_rgba"] = " \
 \
 \
uniform sampler2DRect i1; \
uniform float xCellOffset; \
uniform float yCellOffset; \
 \
void main() { \
   gl_FragColor.r = texture2DRect(i1, gl_TexCoord[0].st).r; \
   gl_FragColor.g = texture2DRect(i1, gl_TexCoord[0].st + vec2(xCellOffset, 0)).r; \
   gl_FragColor.b = texture2DRect(i1, gl_TexCoord[0].st + vec2(0, yCellOffset)).r; \
   gl_FragColor.a = texture2DRect(i1, gl_TexCoord[0].st + vec2(xCellOffset, yCellOffset)).r; \
} \
";

standard_shaders_map["./StandardShaders/Manipulation/SelectChannel-A_gl_frag_rgba"] = " \
 \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   gl_FragColor.r = texture2DRect(i1, gl_TexCoord[0].st).a; \
} \
";

standard_shaders_map["./StandardShaders/Manipulation/SelectChannel-B_gl_frag_rgba"] = " \
 \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   gl_FragColor.r = texture2DRect(i1, gl_TexCoord[0].st).b; \
} \
";

standard_shaders_map["./StandardShaders/Manipulation/SelectChannel-G_gl_frag_rgba"] = " \
 \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   gl_FragColor.r = texture2DRect(i1, gl_TexCoord[0].st).g; \
} \
";

standard_shaders_map["./StandardShaders/Manipulation/SelectChannel-R_gl_frag_rgba"] = " \
 \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   gl_FragColor.r = texture2DRect(i1, gl_TexCoord[0].st).r; \
} \
";

standard_shaders_map["./StandardShaders/Manipulation/transpose-1i0f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT image1, \
	out float4 color : COLOR )  \
{ \
   color = sin(texRECT(image1, texCoord0.yx)); \
} \
";

standard_shaders_map["./StandardShaders/Manipulation/transpose-1i0f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
 \
void main() { \
   gl_FragColor.rgba = texture2DRect(i1, gl_TexCoord[0].ts); \
} \
";

standard_shaders_map["./StandardShaders/Manipulation/UnpackGrayFromRGBA-R_gl_frag_rgba"] = " \
 \
uniform sampler2DRect image; \
 \
void main() { \
   gl_FragColor.r = texture2DRect(image, gl_TexCoord[0].st).r; \
} \
";

standard_shaders_map["./StandardShaders/User/expression1-1i3f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform float f1, \
	uniform float f2,  \
	uniform float f3, \
	out float4 color : COLOR )  \
{ \
   color = acos(sqrt(exp(clamp(sin(pow(cos(texRECT(i1, texCoord0)), f1)), f2, f3)))); \
} \
";

standard_shaders_map["./StandardShaders/User/expression1-1i3f_gl_frag_rgba"] = " \
 \
uniform sampler2DRect i1; \
uniform float f1; \
uniform float f2; \
uniform float f3; \
 \
void main() { \
   vec4 color = texture2DRect(i1, gl_TexCoord[0].st); \
   gl_FragColor.rgba = acos(sqrt(pow(vec4(2.7182818284590451), clamp(sin(pow(cos(color), vec4(f1))), vec4(f2), vec4(f3))))); \
} \
  \
";

standard_shaders_map["./StandardShaders/User/transform-sine-1i2f_cg_frag_rgba"] = " \
 \
void main(in float2 texCoord0 : TEXCOORD0, \
	uniform samplerRECT i1, \
	uniform float f1, // Period \
	uniform float f2, // Amplitude \
	out float4 color : COLOR )  \
{ \
   texCoord0 = float2(texCoord0.x, texCoord0.y + f2 * sin(6.28 * texCoord0.x / f1)); \
 \
   float2 ij = texCoord0; \
   float2 xy = floor(ij); \
   float2 normxy = ij - xy; \
   float2 st0 = ((2.0 - normxy) * normxy - 1.0) * normxy; \
   float2 st1 = (3.0 * normxy - 5.0) * normxy * normxy + 2.0; \
   float2 st2 = ((4.0 - 3.0 * normxy) * normxy + 1.0) * normxy; \
   float2 st3 = (normxy - 1.0) * normxy * normxy; \
 \
   float4 row0 =  \
	st0.s * texRECT(i1, xy + float2(-1.0, -1.0)) +  \
	st1.s * texRECT(i1, xy + float2(0.0, -1.0)) +  \
	st2.s * texRECT(i1, xy + float2(1.0, -1.0)) +  \
	st3.s * texRECT(i1, xy + float2(2.0, -1.0)); \
 \
   float4 row1 =  \
	st0.s * texRECT(i1, xy + float2(-1.0, 0.0)) +  \
	st1.s * texRECT(i1, xy + float2(0.0, 0.0)) +  \
	st2.s * texRECT(i1, xy + float2(1.0, 0.0)) +  \
	st3.s * texRECT(i1, xy + float2(2.0, 0.0)); \
 \
   float4 row2 =  \
	st0.s * texRECT(i1, xy + float2(-1.0, 1.0)) +  \
	st1.s * texRECT(i1, xy + float2(0.0, 1.0)) +  \
	st2.s * texRECT(i1, xy + float2(1.0, 1.0)) +  \
	st3.s * texRECT(i1, xy + float2(2.0, 1.0)); \
 \
   float4 row3 =  \
	st0.s * texRECT(i1, xy + float2(-1.0, 2.0)) +  \
	st1.s * texRECT(i1, xy + float2(0.0, 2.0)) +  \
	st2.s * texRECT(i1, xy + float2(1.0, 2.0)) +  \
	st3.s * texRECT(i1, xy + float2(2.0, 2.0)); \
 \
   color = 0.25 * ((st0.t * row0) + (st1.t * row1) + (st2.t * row2) + (st3.t * row3)); \
} \
";

}
} } // namespace GPU, namespace vw