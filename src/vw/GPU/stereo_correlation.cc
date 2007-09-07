
#include <vw/GPU.h>
using namespace vw;
using namespace GPU;

#define MISSING_PIXEL 1000000

//########################################################################
//#  correlation_iteration - glsl_fragment_strings and a free function                
//########################################################################

char* glsl_string_offset_and_difference = " \
uniform sampler2DRect i1; \
uniform sampler2DRect i2; \
uniform float x_offset; \
uniform float y_offset; \
void main() { \
 vec2 coord = gl_TexCoord[0].st; \
  float v1 = texture2DRect(i1, vec2(coord.s, coord.t)).r; \
  float v2 = texture2DRect(i2, vec2(coord.s + x_offset, coord.t + y_offset)).r; \
  gl_FragColor.r = v1 - v2; \
}";



char* glsl_string_column_sums = " \
uniform sampler2DRect i1; \
void main() { \
 int kHalf = $1; \
 vec2 coord = gl_TexCoord[0].st; \
 float sum = 0; \
 for(int yOffset = -kHalf; yOffset < kHalf+1; yOffset++) { \
     sum += (texture2DRect(i1, vec2(coord.s, coord.t + yOffset))).r; \
 } \
 gl_FragColor.r = sum; \
}";

char* glsl_string_row_sums = " \
uniform sampler2DRect i1; \
void main() { \
 int kHalf = $1; \
 vec2 coord = gl_TexCoord[0].st; \
 float sum = 0; \
 for(int xOffset = -kHalf; xOffset < kHalf+1; xOffset++) { \
     sum += (texture2DRect(i1, vec2(coord.s + xOffset, coord.t))).r; \
 } \
 gl_FragColor.r = sum; \
}";

char* glsl_string_correlation_iteration = " \
uniform sampler2DRect inSums; \
uniform sampler2DRect inBestValues; \
uniform float dx; \
uniform float dy; \
void main() { \
  vec2 coord = gl_TexCoord[0].st; \
  float sum = texture2DRect(inSums, coord).r; \
  vec3 old_bests = texture2DRect(inBestValues, coord).rgb; \
  bool isBest = sum < old_bests.r; \
  if(isBest) { \
    gl_FragColor.rgb = vec3(sum, dx, dy); \
 } \
 else { \
    gl_FragColor.rgb = old_bests; \
 } \
}";




GPUImageBase correlation_iteration(int dx, 
				   int dy, 
				   int kernalHalfSize, 
				   const GPUImageBase& in_left,
				   const GPUImageBase& in_right,
				   GPUImageBase& out_left_bests,  
				   GPUImageBase& out_right_bests)
{
  int width = in_left.width();
  int height = in_left.height();
  // Temp Images  
  GPUImage<PixelGray<float> > temp1(width, height);
  GPUImage<PixelGray<float> > temp2(width, height);
  GPUImageBase temp_bests;
  temp_bests.copy_attributes(out_left_bests);
  
  // *********  STAGE 1 - Offset and Difference **********
  // Setup
  ShaderInvocation_SetupGLState(width, height);
  // Program
  static GPUProgram* program = NULL;
  if(!program)
    program = create_gpu_program_glsl_string(glsl_string_offset_and_difference);
  program->install();
  // Output
  ShaderInvocation_SetOutputImage(temp1);
  // Input
  program->set_uniform_texture("i1", 0, in_left);
  program->set_uniform_texture("i2", 0, in_right);
  program->set_uniform_float("x_offset", -dx);
  program->set_uniform_float("y_offset", -dy);   
  // Drawing
  ShaderInvocation_DrawRectOneTexture(temp1);
  // *********  STAGE 2 - Column Sums ***********
  {
    // Setup
    ShaderInvocation_SetupGLState(width, height);
    // Program - Because we are using the gpu_program_string function, we must do our own caching!
    static std::map<int, GPUProgram*> program_map;
    GPUProgram* program;
    std::map<int, GPUProgram*>::iterator iter = program_map.find(kernalHalfSize);
    if(iter == program_map.end()) {
      vector<int> frag_attributes;
      frag_attributes.push_back(kernalHalfSize);
      program = create_gpu_program_glsl_string(glsl_string_column_sums, frag_attributes);
      program_map[kernalHalfSize] = program;
    }
    else {
      program = (*iter).second;
    }
    program->install();
    // Output
    ShaderInvocation_SetOutputImage(temp2);
    // Input
    program->set_uniform_texture("i1", 0, temp1);
    // Drawing
    ShaderInvocation_DrawRectOneTexture(temp2);  
  }
  // *********  STAGE 3 - Row Sums ***********
  {
    // Setup
    ShaderInvocation_SetupGLState(width, height);
    // Program - Because we are using the gpu_program_string function, we must do our own caching!
    static std::map<int, GPUProgram*> program_map;
    GPUProgram* program;
    std::map<int, GPUProgram*>::iterator iter = program_map.find(kernalHalfSize);
    if(iter == program_map.end()) {
      vector<int> frag_attributes;
      frag_attributes.push_back(kernalHalfSize);
      program = create_gpu_program_glsl_string(glsl_string_row_sums, frag_attributes);
      program_map[kernalHalfSize] = program;
    }
    else {
      program = (*iter).second;
    }
    program->install();
    // Output
    ShaderInvocation_SetOutputImage(temp1);
    // Input
    program->set_uniform_texture("i1", 0, temp2);
    // Drawing
    ShaderInvocation_DrawRectOneTexture(temp1);  
  }
  // *********  STAGE 4A/B - Update Best Values for Left and Right ***********
  {
    static GPUProgram* program = NULL;
    if(!program)
      program = create_gpu_program_glsl_string(glsl_string_correlation_iteration);
    // LEFT SIDE: Setup
    ShaderInvocation_SetupGLState(width, height);
    // Program
    program->install();
    // Output
    ShaderInvocation_SetOutputImage(temp_bests);
    // Input
    program->set_uniform_texture("inSums", 0, temp1);
    program->set_uniform_texture("inBestValues", 1, out_left_bests);
    program->set_uniform_float("dx", dx);
    program->set_uniform_float("dy", dy);
    // Drawing
    ShaderInvocation_DrawRectOneTexture(temp_bests);
    // Copy to output
    out_left_bests = copy(temp_bests);

    // RIGHT SIDE: Setup
    ShaderInvocation_SetupGLState(width, height);
    // Program
    program->install();
    // Output
    ShaderInvocation_SetOutputImage(temp_bests);
    // Input
    program->set_uniform_texture("inSums", 0, temp1);
    program->set_uniform_texture("inBestValues", 1, out_right_bests);
    program->set_uniform_float("dx", -dx);
    program->set_uniform_float("dy", -dy);
    // Drawing
    ShaderInvocation_DrawRectOneTexture(temp_bests);
    // Copy to output
    out_right_bests = temp_bests;
  }			
}

//########################################################################
//#  correlation_cross_check                                             #
//########################################################################

char* glsl_string_correlation_cross_check = " \
uniform sampler2DRect inLeftBestValues; \
uniform sampler2DRect inRightBestValues; \
uniform float crossCheckThreshold; \
uniform float missingPixel; \
void main() { \
 vec3 leftValues = texture2DRect(inLeftBestValues, gl_TexCoord[0].st).rgb; \
 vec2 rightCorresponding = texture2DRect(inRightBestValues, vec2(gl_TexCoord[0].s + leftValues.g, gl_TexCoord[0].t + leftValues.b )).gb; \
 float xAbsDifference = abs(leftValues.g + rightCorresponding.r); \
 float yAbsDifference = abs(leftValues.b + rightCorresponding.g); \
 bool withinThreshhold = xAbsDifference < crossCheckThreshold && yAbsDifference < crossCheckThreshold; \
 if(withinThreshhold) { \
    gl_FragColor.rgb = leftValues.rgb; \
 } \
 else { \
    gl_FragColor.rgb = vec3(0.0, missingPixel, missingPixel); \
 } \
}";

GPUImageBase correlation_cross_check(float threshold, GPUImageBase best_values_left, GPUImageBase best_values_right) {
  // Setup
  ShaderInvocation_SetupGLState(best_values_left.width(), best_values_left.height());
  // Program
  static GPUProgram* program = NULL;
  if(!program)
    program = create_gpu_program_glsl_string(glsl_string_correlation_cross_check);
  program->install();
  // Output
  GPUImageBase output;
  output.copy_attributes(best_values_left);
  ShaderInvocation_SetOutputImage(output);
  // Input
  program->set_uniform_texture("inLeftBestValues", 0, best_values_left);
  program->set_uniform_texture("inRightBestValues", 1, best_values_right);
  program->set_uniform_float("crossCheckThreshold", threshold);
  program->set_uniform_float("missingPixel", MISSING_PIXEL);
  // Drawing
  ShaderInvocation_DrawRectOneTexture(output);
  // Return
  return output;
}


//########################################################################
//#  stereo_correlation                   
//########################################################################

GPUImage<PixelRGB<float> > stereo_correlation(const GPUImageBase &leftImage,
				const GPUImageBase &rightImage,
				GPUImageBase* outCorrelationScores,
				int kernalHalfSize,
				int minDX,
				int maxDX,
				int minDY,
				int maxDY,
				float crossCorrThreshold)  
{
  int nShifts = (1 + maxDX - minDX) * (1 + maxDY - minDY);
  int nPixels = leftImage.width() * leftImage.height();  
  // Convolution Matrices
  float matrix1[] = 
    {   1, 1, 1,
	1, 1, 1,
	1, 1, 1,
    }; 
  GPUImageBase gaussian_matrix(3, 3, TEX_R, TEX_FLOAT32, TEX_R, TEX_FLOAT32, matrix1);

  float matrix2[] = 
    {   0, -1, 0,
	-1, 4, -1,
	0, -1, 0,
    };
  GPUImageBase laplacian_matrix(3, 3, TEX_R, TEX_FLOAT32, TEX_R, TEX_FLOAT32, matrix2);
  // Left Image - Input / SLOG
  GPUImageBase temp_left = convolution_filter(leftImage, gaussian_matrix);
  temp_left = convolution_filter(temp_left, laplacian_matrix);
  temp_left = threshold(temp_left, 0, -1, 1);
  // Right Image - Input / SLOG  
  GPUImageBase temp_right = convolution_filter(rightImage, gaussian_matrix);
  temp_right = convolution_filter(temp_right, laplacian_matrix);
  temp_right = threshold(temp_right, 0, -1, 1);
  // Best Scores
  GPUImage<PixelRGB<float> > best_values_left(leftImage.width(), leftImage.height());
  fill(best_values_left, 1000000, 0, 0, 0); 

  GPUImage<PixelRGB<float> > best_values_right(leftImage.width(), leftImage.height());
  fill(best_values_right, 1000000, 0, 0, 0);
  // LOOP
  int c=0;
  for(int DX = minDX; DX < maxDX+1; DX++) {
    for(int DY = minDY; DY < maxDY+1; DY++) {       
      correlation_iteration(DX, 
			    DY, 
			    kernalHalfSize, 
			    temp_left, 
			    temp_right,
			    best_values_left,
			    best_values_right);

    }
  }
  
  // CROSS CHECK
  GPUImageBase best_values_cross = 
    correlation_cross_check(crossCorrThreshold, best_values_left, best_values_right);

  return best_values_cross;

  // Read Back Output
  // TO DO - Need to decide if it should return it as GPUImage or ImageView or DisparityMap
  /*
    int oWidth = bestSOADCross_Texture.Width();
    int oHeight = bestSOADCross_Texture.Height();

    vil_image_view<float> DHImage_FView(oWidth, oHeight);
    bestSOADCross_Texture.Read(TEX_G, TEX_FLOAT32, DHImage_FView.top_left_ptr());

    vil_image_view<float> DVImage_FView(oWidth, oHeight);
    bestSOADCross_Texture.Read(TEX_B, TEX_FLOAT32, DVImage_FView.top_left_ptr());

    DisparityMap disparityMap(oWidth, oHeight);
    for(int y=0; y < DHImage_FView.nj(); y++) {
    for(int x=0; x < DHImage_FView.ni(); x++) {
    disparityMap(x, y, 0) = DHImage_FView(x, y);
    disparityMap(x, y, 1) = DVImage_FView(x, y);
    }
    }
  
    printf("returning...\n");
    return DisparityMap();
  */
}

int main(int argc, char *argv[]) {
  // Parse args: { input image, output image }
  if(argc < 4) {
    printf("Error: You must provide at least 3 arguments. Args List:\n");
    printf("  * left_image \n  * right_image \n  * output_image \n  * kernal_size = 21\n  * min_dx = 0 \n\
  * max_dx = 30 \n  * min_dy = -5 \n  * max_dy = 5 \n  * error_threshold = 1.0\n");
    return -1;
  }
  string left_path = argv[1];
  string right_path = argv[2];
  string output_path = argv[3];

  int kernalSize = (argc >= 5) ? atoi(argv[4]) : 21;
  int minDX = (argc >= 6) ? atoi(argv[5]) : 0;
  int maxDX = (argc >= 7) ? atoi(argv[6]) : 30;
  int minDY = (argc >= 8) ? atoi(argv[7]) : -5;
  int maxDY = (argc >= 9) ? atoi(argv[8]) : 5;
  float crossCorrThreshold = (argc >= 10) ? atof(argv[9]) : 1.0;

  printf("  * left_image = %s\n  * right_image = %s\n  * output_image = %s\n  * kernal_size = %i\n  * min_dx = %i \n\
  * max_dx = %i \n  * min_dy = %i \n  * max_dy = %i \n  * error_threshold = %f\n",
	 left_path.c_str(), right_path.c_str(), output_path.c_str(), kernalSize, minDX, maxDX, minDY, maxDY, crossCorrThreshold);

  // Init
  gpu_init(true, true);
  set_shader_language_choice(SHADER_LANGUAGE_CHOICE_GLSL);
    // Read image
  ImageView<PixelGray<float> > left_image, right_image;
  read_image(left_image, left_path);
  read_image(right_image, right_path);

  if(!left_image.cols() || !right_image.cols()) {
    printf("Error: Image input error.\n");
    return -1;
  }
  //
  GPUImage<PixelRGB<float> > output = stereo_correlation((GPUImage<PixelGray<float> >) left_image,
							 (GPUImage<PixelGray<float> >) right_image,
							 NULL,
							 (int) floorf(kernalSize / 2.0),
							 minDX,
							 maxDX,
							 minDY,
							 maxDY,
							 crossCorrThreshold);

  ImageView<PixelGray<float> > dx_image(output.cols(), output.rows());
  output.read(TEX_G, TEX_FLOAT32, &(dx_image(0,0)));
  normalize(dx_image, 0, maxDX-minDX);
  write_image(output_path, dx_image);
}
