// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/GPU.h>
using namespace vw;
using namespace GPU;

using std::vector;
using std::string;

//########################################################################
//#  correlation_iteration - glsl_frag strings and a free function
//########################################################################

const char* glsl_frag_offset_and_difference = " \
uniform sampler2DRect i1; \
uniform sampler2DRect i2; \
uniform float x_offset; \
uniform float y_offset; \
void main() { \
 vec2 coord = gl_TexCoord[0].st; \
  float v1 = texture2DRect(i1, vec2(coord.s, coord.t)).r; \
  float v2 = texture2DRect(i2, vec2(coord.s + x_offset, coord.t + y_offset)).r; \
  gl_FragColor.r = abs(v1 - v2); \
}";



const char* glsl_frag_column_sums = " \
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

const char* glsl_frag_row_sums = " \
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

const char* glsl_frag_correlation_iteration = " \
uniform sampler2DRect inSums; \
uniform sampler2DRect inBestValues; \
uniform float dx; \
uniform float dy; \
uniform float xOffset; \
uniform float yOffset; \
void main() { \
  vec2 coord = gl_TexCoord[0].st; \
  float sum = texture2DRect(inSums, vec2(coord.s + xOffset, coord.t + yOffset)).r; \
  vec3 old_bests = texture2DRect(inBestValues, coord).rgb; \
  bool isBest = sum < old_bests.r; \
  if(isBest) { \
    gl_FragColor.rgb = vec3(sum, dx, dy); \
 } \
 else { \
    gl_FragColor.rgb = old_bests; \
 } \
}";

void correlation_iteration(int dx,
                           int dy,
                           int kernalHalfSize,
                           const GPUImage<PixelGray<float> >& in_left,
                           const GPUImage<PixelGray<float> >& in_right,
                           GPUImageBase& out_left_bests,
                           GPUImageBase& out_right_bests)
{
  int width = in_left.width();
  int height = in_left.height();
  // Temporary Images
  GPUImage<PixelGray<float> > temp1(width, height);
  GPUImage<PixelGray<float> > temp2(width, height);
  GPUImageBase temp_bests;
  temp_bests.copy_attributes(out_left_bests);
  // *********  STAGE 1 - Offset and Absolute Difference **********
  {
    ShaderInvocation_SetupGLState(width, height);
    static GPUProgram* program = NULL;
    if(!program)
      program = create_gpu_program_glsl_string(glsl_frag_offset_and_difference);
    program->install();
    ShaderInvocation_SetOutputImage(temp1);
    program->set_input_image("i1", in_left);
    program->set_input_image("i2", in_right);
    program->set_input_float("x_offset", -dx);
    program->set_input_float("y_offset", -dy);
    ShaderInvocation_DrawRectOneTexture(temp1);
  }
  // write_image("OUTPUT_TEST_Difference.png", (GPUImage<PixelGray<float> >) temp1);

  // *********  STAGE 2 - Column Sum ***********
  {
    ShaderInvocation_SetupGLState(width, height);
    // Get Program and cache the program based on kernal size (string creation function doesn't cache)
    static std::map<int, GPUProgram*> program_map;
    GPUProgram* program;
    std::map<int, GPUProgram*>::iterator iter = program_map.find(kernalHalfSize);
    if(iter == program_map.end()) {
      vector<int> frag_attributes;
      frag_attributes.push_back(kernalHalfSize);
      program = create_gpu_program_glsl_string(glsl_frag_column_sums, frag_attributes);
      program_map[kernalHalfSize] = program;
    }
    else {
      program = (*iter).second;
    }
    // Execute
    program->install();
    ShaderInvocation_SetOutputImage(temp2);
    program->set_input_image("i1", temp1);
    ShaderInvocation_DrawRectOneTexture(temp2);
  }
  // *********  Stage 3 - Row Sum ***********
  {
    ShaderInvocation_SetupGLState(width, height);
    // Get Program and cache the program based on kernal size (string creation function doesn't cache)
    static std::map<int, GPUProgram*> program_map;
    GPUProgram* program;
    std::map<int, GPUProgram*>::iterator iter = program_map.find(kernalHalfSize);
    if(iter == program_map.end()) {
      vector<int> frag_attributes;
      frag_attributes.push_back(kernalHalfSize);
      program = create_gpu_program_glsl_string(glsl_frag_row_sums, frag_attributes);
      program_map[kernalHalfSize] = program;
    }
    else {
      program = (*iter).second;
    }
    // Execute
    program->install();
    ShaderInvocation_SetOutputImage(temp1);
    program->set_input_image("i1", temp2);
    ShaderInvocation_DrawRectOneTexture(temp1);
  }
  // *********  STAGE 4 L/R - Update Best Values for Left and Right ***********
  {
    static GPUProgram* program = NULL;
    if(!program)
      program = create_gpu_program_glsl_string(glsl_frag_correlation_iteration);
    GPUProgram* program_copy = create_gpu_program("VWGPU_Algorithms/copy-1i0f");

    // LEFT SIDE:
    // Copy invalid region:
    ShaderInvocation_SetupGLState(width, height);
    program_copy->install();
    ShaderInvocation_SetOutputImage(temp_bests);
    program_copy->set_input_image("i1", out_left_bests);

    int left_bound = dx;
    glBegin(GL_QUADS);
    glTexCoord2f(0, 0);         glVertex2f(0, 0);
    glTexCoord2f(left_bound, 0);      glVertex2f(left_bound, 0);
    glTexCoord2f(left_bound, height);  glVertex2f(left_bound, height);
    glTexCoord2f(0, height);     glVertex2f(0, height);
    glEnd();
    // Update best values in valid region
    ShaderInvocation_SetupGLState(width, height);
    program->install();
    ShaderInvocation_SetOutputImage(temp_bests);
    program->set_input_image("inSums", temp1);
    program->set_input_image("inBestValues", out_left_bests);
    program->set_input_float("dx", dx);
    program->set_input_float("dy", dy);
    program->set_input_float("xOffset", 0);
    program->set_input_float("yOffset", 0);

    glBegin(GL_QUADS);
    glTexCoord2f(left_bound, 0);
    glVertex2f(left_bound, 0);
    glTexCoord2f(width, 0);
    glVertex2f(width, 0);
    glTexCoord2f(width, height);
    glVertex2f(width, height);
    glTexCoord2f(left_bound, height);
    glVertex2f(left_bound, height);
    glEnd();

    out_left_bests = copy(temp_bests);
    // RIGHT SIDE:
    // Copy invalid region
    ShaderInvocation_SetupGLState(width, height);
    program_copy->install();
    ShaderInvocation_SetOutputImage(temp_bests);
    program_copy->set_input_image("i1", out_right_bests);

    int right_bound = width - dx - kernalHalfSize;
    glBegin(GL_QUADS);
    glTexCoord2f(right_bound, 0);         glVertex2f(right_bound, 0);
    glTexCoord2f(width, 0);      glVertex2f(width, 0);
    glTexCoord2f(width, height);  glVertex2f(width, height);
    glTexCoord2f(right_bound, height);     glVertex2f(right_bound, height);
    glEnd();
    // Update best values in valid region
    ShaderInvocation_SetupGLState(width, height);
    program->install();
    ShaderInvocation_SetOutputImage(temp_bests);
    program->set_input_image("inSums", temp1);
    program->set_input_image("inBestValues", out_right_bests);
    program->set_input_float("dx", dx);
    program->set_input_float("dy", dy);
    program->set_input_float("xOffset", dx);
    program->set_input_float("yOffset", dy);

    glBegin(GL_QUADS);
    glTexCoord2f(0, 0);
    glVertex2f(0, 0);
    glTexCoord2f(right_bound, 0);
    glVertex2f(right_bound, 0);
    glTexCoord2f(right_bound, height);
    glVertex2f(right_bound, height);
    glTexCoord2f(0, height);
    glVertex2f(0, height);
    glEnd();

    out_right_bests = temp_bests;
  }
}

//########################################################################
//#  correlation_cross_check                                             #
//########################################################################

#define MISSING_PIXEL 10000

const char* glsl_frag_correlation_cross_check = " \
uniform sampler2DRect inLeftBestValues; \
uniform sampler2DRect inRightBestValues; \
uniform float crossCheckThreshold; \
uniform float missingPixel; \
void main() { \
 vec3 leftValues = texture2DRect(inLeftBestValues, gl_TexCoord[0].st).rgb; \
 vec3 rightValues = texture2DRect(inRightBestValues, vec2(gl_TexCoord[0].s - leftValues.g, gl_TexCoord[0].t - leftValues.b )).rgb; \
 float xAbsDifference = abs(leftValues.g - rightValues.g); \
 float yAbsDifference = abs(leftValues.b - rightValues.b); \
 bool withinThreshhold = xAbsDifference < crossCheckThreshold && yAbsDifference < crossCheckThreshold; \
 if(withinThreshhold) { \
    gl_FragColor.rgb = leftValues.rgb; \
 } \
 else { \
    gl_FragColor.rgb = vec3(missingPixel, missingPixel, missingPixel); \
 } \
}";

GPUImageBase correlation_cross_check(float threshold, GPUImageBase best_values_left, GPUImageBase best_values_right) {
  ShaderInvocation_SetupGLState(best_values_left.width(), best_values_left.height());
  static GPUProgram* program = NULL;
  if(!program)
    program = create_gpu_program_glsl_string(glsl_frag_correlation_cross_check);
  program->install();

  GPUImageBase output;
  output.copy_attributes(best_values_left);
  ShaderInvocation_SetOutputImage(output);
  program->set_input_image("inLeftBestValues", best_values_left);
  program->set_input_image("inRightBestValues", best_values_right);
  program->set_input_float("crossCheckThreshold", threshold);
  program->set_input_float("missingPixel", MISSING_PIXEL);
  ShaderInvocation_DrawRectOneTexture(output);

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
  // Convolution Matrices
  float matrix1[] =
    {   1, 2, 1,
        2, 4, 2,
        1, 2, 1,
    };
  GPUImageBase gaussian_matrix(3, 3, GPU_RED, GPU_FLOAT32, GPU_RED, GPU_FLOAT32, matrix1);

  float matrix2[] =
    {   0, -1, 0,
        -1, 4, -1,
        0, -1, 0,
    };
  GPUImageBase laplacian_matrix(3, 3, GPU_RED, GPU_FLOAT32, GPU_RED, GPU_FLOAT32, matrix2);
  // Left Image - SLOG
  GPUImageBase temp_left = convolution_filter(leftImage, gaussian_matrix);
  temp_left = convolution_filter(temp_left, laplacian_matrix);
  temp_left = threshold(temp_left, 0, 0, 1);
  // Right Image - SLOG
  GPUImageBase temp_right = convolution_filter(rightImage, gaussian_matrix);
  temp_right = convolution_filter(temp_right, laplacian_matrix);
  temp_right = threshold(temp_right, 0, 0, 1);
  // Init Best Scores
  GPUImage<PixelRGB<float> > best_values_left(leftImage.width(), leftImage.height());
  fill(best_values_left, 1000000, 1000000, 1000000, 0);

  GPUImage<PixelRGB<float> > best_values_right(leftImage.width(), leftImage.height());
  fill(best_values_right, 1000000, 1000000, 1000000, 0);
  // Main Loop
  for(int dx = minDX; dx <= maxDX; dx++) {
    for(int dy = minDY; dy <= maxDY; dy++) {
      correlation_iteration(dx,
                            dy,
                            kernalHalfSize,
                            temp_left,
                            temp_right,
                            best_values_left,
                            best_values_right);

    }
  }
  // Perform Cross-Check and Return
  GPUImageBase best_values_cross =
    correlation_cross_check(crossCorrThreshold, best_values_left, best_values_right);
  return best_values_cross;
}


//########################################################################
//#  main
//########################################################################

int main(int argc, char *argv[]) {
  // Parse arguments
  if(argc < 3) {
    printf("Error: The first two arguments are required. \nArgument list:   [ NAME (TYPE, DEFAULT) ]\n");
    printf("  * left_image (string, REQUIRED) \n  * right_image (string, REQUIRED) \n  * kernal_size (int, 21)\n  * min_dx (int, 0) \n\
  * max_dx (int, 100) \n  * min_dy (int, 0) \n  * max_dy (int, 0) \n  * error_threshold (float, 1.0)\n\
  * output_image_dx (string, \'./Output_DX.png\') - Normalized from [mix_dx to max_dx]\n\
  * output_image_dy (string, '') - Normalized from [mix_dy to max_dy]\n\
  * output_image_error (string, '') - Normalized from [0 to highest possible error score]\n");
    return -1;
  }
  string left_path = argv[1];
  string right_path = argv[2];
  int kernalSize = (argc >= 4) ? atoi(argv[3]) : 21;
  int minDX = (argc >= 5) ? atoi(argv[4]) : 0;
  int maxDX = (argc >= 6) ? atoi(argv[5]) : 30;
  int minDY = (argc >= 7) ? atoi(argv[6]) : -5;
  int maxDY = (argc >= 8) ? atoi(argv[7]) : 5;
  float crossCorrThreshold = (argc >= 9) ? atof(argv[8]) : 1.0;
  string output_path_dx = (argc >= 10) ? argv[9] : "Output_DX.png";
  string output_path_dy = (argc >= 11) ? argv[10] : "";
  string output_path_score = (argc >= 12) ? argv[11] : "";

  printf("  * left_image = \'%s\'\n  * right_image = \'%s\'\n  * kernal_size = %i\n  * min_dx = %i \n\
  * max_dx = %i \n  * min_dy = %i \n  * max_dy = %i \n  * error_threshold = %f\n\
  * output_image_dx = \'%s\' \n  * output_image_dy = \'%s\' \n  * output_image_error = \'%s\' \n",
         left_path.c_str(), right_path.c_str(), kernalSize, minDX, maxDX, minDY, maxDY, crossCorrThreshold,
         output_path_dx.c_str(), output_path_dy.c_str(), output_path_score.c_str());

  // Init VWGPU
  gpu_init("Log_VWGPU.txt");
  set_shader_language_choice(SHADER_LANGUAGE_CHOICE_GLSL);
  set_gpu_memory_recycling(true);
  // Read images
  ImageView<PixelGray<float> > left_image, right_image;
  read_image(left_image, left_path);
  read_image(right_image, right_path);

  if(!left_image.cols() || !right_image.cols()) {
    printf("ERROR: Image input error.\n");
    return -1;
  }
  // Stereo Correlation
  GPUImage<PixelRGB<float> > best_values = stereo_correlation((GPUImage<PixelGray<float> >) left_image,
                                                              (GPUImage<PixelGray<float> >) right_image,
                                                              NULL,
                                                              (int) floorf(kernalSize / 2.0),
                                                              minDX,
                                                              maxDX,
                                                              minDY,
                                                              maxDY,
                                                              crossCorrThreshold);
  // Slice requested channels and return as image files (RED = best_score, GREEN = best_dx, BLUE = best_dy)
  ImageView<PixelGray<float> > output_image(best_values.cols(), best_values.rows());

  if(!output_path_dx.empty()) {
    best_values.read(GPU_GREEN, GPU_FLOAT32, &(output_image(0,0)));
    output_image = output_image - PixelGray<float>(minDX);
    output_image = output_image / PixelGray<float>(maxDX - minDX);
    try {
      write_image(output_path_dx, output_image);
    }
    catch(...) {
      printf("ERROR: Problem writing output_image_dx\n");
      return -1;
    }
  }
  if(!output_path_dy.empty()) {
    best_values.read(GPU_BLUE, GPU_FLOAT32, &(output_image(0,0)));
    output_image = output_image - PixelGray<float>(minDY);
    output_image = output_image / PixelGray<float>(maxDY - minDY);
    try {
      write_image(output_path_dy, output_image);
    }
    catch(...) {
      printf("ERROR: Problem writing output_image_dy\n");
      return -1;
    }
  }
  if(!output_path_score.empty()) {
    int max_score = (int) powf(1 + 2 * floorf(kernalSize / 2.0), 2);
    best_values.read(GPU_RED, GPU_FLOAT32, &(output_image(0,0)));
    output_image = output_image / PixelGray<float>(max_score);
    try {
      write_image(output_path_score, output_image);
    }
    catch(...) {
      printf("ERROR: Problem writing output_image_error\n");
      return -1;
    }
  }
  // Clean up and return
  gpu_cleanup();
  return 0;
}
