// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef EXPRESSION_H
#define EXPRESSION_H

#include <list>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <memory>

#include <vw/GPU/Setup.h>
#include <vw/Math/Matrix.h>

namespace vw { namespace GPU {

  static const char* expressionShaderPath = "expression_shaders/";

  /*
    Does Rasterize() method delete the object and return a new one?
  */

  class GPUImageBase;
  typedef GPUImageBase ShaderNode_Image;

// ########################################################################################################
//                                           ShaderBindings struct
// ########################################################################################################


 class ShaderBindings {
 public:
   std::list<float> floats;
   std::list<std::vector<float> > vectors;

   std::list<int> imageBindingList;
   int floatBindingStartIndex;
   int vectorBindingStartIndex;
   //list<matrix<float, 4, 4> > matrices;
 };

// ########################################################################################################
//                                           ShaderBuilder
// ########################################################################################################

class ShaderBuilder {
  bool isCG;
  int cFloatIndex;
  int cVectorIndex;
  int cImageIndex;
  int cFunctionIndex;

  std::list<float> floatBindingList;
  std::list<std::vector<float> > vectorBindingList;
  std::list<GPUImageBase> imageBindingList;
  std::list<std::pair<bool, Matrix<float> > > homographyBindingList;
  std::map<int, int> textureBindingMap;

  std::string bodyText;
  std::string functionsText;

 public:
  ShaderBuilder() {
    cImageIndex = -1;
    cFloatIndex = -1;
    cVectorIndex = -1;
    cFunctionIndex = -1;
    isCG = false;
  }

  void AppendBody(const std::string& str) {
    bodyText += str;
  }

  void RemoveLastBodyChar() {
    bodyText.resize(bodyText.size() - 1);
  }

  void AppendFunctions(const std::string& str) {
    functionsText += str;
    functionsText += "\n\n";
  }

  void AddFunctionFromFile(const std::string& shaderBaseName) {
    // make path
    std::string path = expressionShaderPath;
    path += shaderBaseName;
    if(!isCG)
      path += "_gl_";
    else
      path += "_cg_";
    path += "rgba";
    // read file
    std::string content;
    if(!ReadFileAsString(path, content)) {
      printf("[ShaderBuilder::AddFunctionFromFile] Couldn't find file: %s\n", path.c_str());
      exit(0);
    }
    // change name
    // add function
    AppendFunctions(content);
  }

  bool Complete() {
    std::stringstream shader;
    if(!isCG) {
      for(int i=0; i <= cImageIndex; i++)
        shader << "sampler2DRect f" << i << ";\n";
      for(int i=0; i <= cFloatIndex; i++)
        shader << "float f" << i << ";\n";
      for(int i=0; i <= cVectorIndex; i++)
        shader << "vec4 f" << i << ";\n";
      shader << "main()\n" << bodyText << ";\n" << "}\n\n" << functionsText;
    }
    printf("\n****************************************************\n%s\n\n", shader.str().c_str());
    return true;
  }

  // BINDINGS


  int BindFloat(float inFloat) { floatBindingList.push_back(inFloat); return ++cFloatIndex; }

  int BindVector(std::vector<float>& inVector) { vectorBindingList.push_back(inVector); return ++cVectorIndex; }

  int BindFunction() { return ++cFunctionIndex; }

  int BindImage(GPUImageBase& image) {
    /*
    int name = image.name();
    std::map<int, int>::iterator findIter = textureBindingMap.find(name);
    if(findIter == textureBindingMap.end()) {
      cImageIndex++;
      textureBindingMap[name] = cImageIndex;
      return cImageIndex;
    }
    else {
      return (*findIter).second;
      } */
  }

  void BindGroup(ShaderBindings& shaderBindings) {
    int i = 0;
    for(std::list<float>::iterator iter = shaderBindings.floats.begin(); iter != shaderBindings.floats.end(); iter++, i++) {
      int index = BindFloat(*iter);
      if(i == 0)
        shaderBindings.floatBindingStartIndex = index;
    }

    i = 0;
    for(std::list<std::vector<float> >::iterator iter = shaderBindings.vectors.begin(); iter != shaderBindings.vectors.end(); iter++, i++) {
      int index = BindVector(*iter);
      if(i == 0)
        shaderBindings.vectorBindingStartIndex = index;
    }
  }

  bool BindingComplete() {
    return true;
  }

  ShaderNode_Image* Evaluate() {
    return NULL;
  }

};

 class ShaderNode_Accessor;

//// Utility
 ShaderNode_Accessor* ShaderNodeCreateDefaultInterpolator();


// ########################################################################################################
//                                           ShaderNode_Base
// ########################################################################################################

 class ShaderNode_Base {
 public:
  std::auto_ptr<ShaderNode_Base> previous;
  ShaderBindings shaderBindings;
  int m_width;
  int m_height;
 public:
  int width() { return m_width; }
  int height() { return m_height; }
  void SetSize(int w, int h) { m_width = w; m_height = h; }
  void Disconnect() { previous.release(); }
  ShaderNode_Image* Evaluate();
  virtual ~ShaderNode_Base() { }
  virtual void AddToShader(ShaderBuilder& shaderBuilder) { }
  virtual int CalculateGraphImageCount() { return 0; }
  virtual int CalculateGraphFloatCount() { return 0; }
  virtual int RecursiveBind(ShaderBuilder& shaderBuilder) { return 0; }
  virtual void RecursiveCompose(ShaderBuilder& shaderBuilder) { }
};


// ########################################################################################################
//                                           BaseTransform - Linear / Nonlinear Subclasses
// ########################################################################################################

class ShaderNode_BaseTransform : public ShaderNode_Base {
 public:
  std::auto_ptr<ShaderNode_Accessor> m_interpolator;
 public:
  virtual ~ShaderNode_BaseTransform() { }

  ShaderNode_Base* AppendInterpolateNode();

};


class ShaderNode_LinearTransform : public ShaderNode_BaseTransform {
 public:
    Matrix<float> m_homography;
 public:

  ShaderNode_LinearTransform(ShaderNode_Base* node, Matrix<float>& inHomography, ShaderNode_Accessor* interpolator = NULL);

  ~ShaderNode_LinearTransform() { }

  int RecursiveBind(ShaderBuilder& shaderBuilder);

  void RecursiveCompose(ShaderBuilder& shaderBuilder);
};

class ShaderNode_NonlinearTransform : public ShaderNode_BaseTransform {
  std::string m_shader;
  ShaderBindings m_bindings;
 public:
  ShaderNode_NonlinearTransform(ShaderNode_Base* node, const std::string& str, const ShaderBindings& bindings, ShaderNode_Accessor* interpolator = NULL);

  ~ShaderNode_NonlinearTransform() { }

  int RecursiveBind(ShaderBuilder& shaderBuilder);

  void RecursiveCompose(ShaderBuilder& shaderBuilder);

};




// ########################################################################################################
//                                           Accessor
// ########################################################################################################

class ShaderNode_Accessor : public ShaderNode_Base {
 public:
  std::list<ShaderNode_Base*> m_input_nodes;
  std::string m_shader;
  ShaderBindings m_bindings;
  std::auto_ptr<ShaderNode_Accessor> m_interpolator;
  std::list<int> m_image_indices;
  int m_quality;
 public:
  ShaderNode_Accessor(std::list<ShaderNode_Base*>& input_nodes, const std::string& shader, const ShaderBindings& bindings, bool forceRasterize = false, bool isInterpolateNode = false);

  ~ShaderNode_Accessor() { for(std::list<ShaderNode_Base*>::iterator iter = m_input_nodes.begin(); iter != m_input_nodes.end(); iter++) delete *iter; }

  int RecursiveBind(ShaderBuilder& shaderBuilder);

  void RecursiveCompose(ShaderBuilder& shaderBuilder);

};

// ########################################################################################################
//                                           PixelModifier
// ########################################################################################################

class ShaderNode_PixelModifier : public ShaderNode_Base {
 public:
  std::list<ShaderNode_Base*> m_input_nodes;
  std::string m_shader;
  ShaderBindings m_bindings;
  std::auto_ptr<ShaderNode_Accessor> m_interpolator;
 public:
  ShaderNode_PixelModifier(std::list<ShaderNode_Base*>& input_nodes, const std::string& shader, const ShaderBindings& bindings);

  ~ShaderNode_PixelModifier() { for(std::list<ShaderNode_Base*>::iterator iter = m_input_nodes.begin(); iter != m_input_nodes.end(); iter++) delete *iter; }

  int RecursiveBind(ShaderBuilder& shaderBuilder);

  void RecursiveCompose(ShaderBuilder& shaderBuilder);

};


// ########################################################################################################
//                                           test function
// ########################################################################################################

/*

 ShaderNode_LinearTransform* scale(ShaderNode_Base* image, float factor) { // Linear Transform
   Matrix<float> homography(3, 3);
   homography.set_identity();
   homography(0, 0) = factor;
   homography(1, 1) = factor;
   return new ShaderNode_LinearTransform(image, homography);
 }

 ShaderNode_LinearTransform* translate(ShaderNode_Base* image, float x, float y) { // Linear Transform
   Matrix<float> homography(3, 3);
   homography.set_identity();
   homography(0, 2) = x;
   homography(1, 2) = y;
   return new ShaderNode_LinearTransform(image, homography);
 }

 ShaderNode_NonlinearTransform* radial(ShaderNode_Base* image, float period, float amplitude) { // Nonlinear Transform
   ShaderBindings bindings;
   bindings.floats.push_back(period);
   bindings.floats.push_back(amplitude);
   return new ShaderNode_NonlinearTransform(image, "NonlinearTransform/sine_wave", bindings);
 }
*/
 /*
 ShaderNode_Base* operator+(ShaderNode_Base* image, float scalar) {  // PixelModifier
   std::list<ShaderNode_Base*> nodeList;
   nodeList.push_back(image);
   ShaderBindings bindings;
   bindings.floats.push_back(scalar);
   return new ShaderNode_PixelModifier(nodeList, "PixelModifier/sum-IF", bindings);
 }
 ShaderNode_PixelModifier* operator*(ShaderNode_Base* image1, ShaderNode_Base* image2) {  // PixelModifier
   std::list<ShaderNode_Base*> nodeList;
   nodeList.push_back(image1);
   nodeList.push_back(image2);
   ShaderBindings bindings;
   return new ShaderNode_PixelModifier(nodeList, "PixelModifier/multiply-II", ShaderBindings());
 }

 */




 /*

  // Nonlinear Transform
  ShaderNode_LinearTransform(ShaderNode_NonlinearTransform* node, Matrix<float>& inHomography) {
    if(interpolator)
      m_interpolator = interpolator;
    else
      m_interpolator = ShaderNodeCreateDefaultInterpolator();
    if(node->m_interpolator->m_quality > m_interpolator->m_quality)
      m_interpolator.reset(node->m_interpolator);
    previous.reset(node);
    m_homography = inHomography;
  }
  // Accessor
  ShaderNode_LinearTransform(ShaderNode_Accessor* node, Matrix<float>& inHomography, ShaderNode_Accessor* interpolator = NULL) {
    if(interpolator)
      m_interpolator = interpolator;
    else
      m_interpolator = ShaderNodeCreateDefaultInterpolator();
    previous.reset(node->Rasterize());
    homography = inHomography;
  }
  // PixelModifier
  ShaderNode_LinearTransform(ShaderNode_PixelModifier* node, Matrix<float>& inHomography, ShaderNode_Accessor* interpolator = NULL) {
    if(interpolator)
      m_interpolator = interpolator;
    else
      m_interpolator = ShaderNodeCreateDefaultInterpolator();
    previous.reset(node->Rasterize());
    homography = inHomography;
  }

 */



  /*
class ShaderNode_NonlinearTransform : public ShaderNode_Base {
 public:
  std::string m_shader;
  ShaderBindings m_bindings;
  std::auto_ptr<ShaderNode_Accessor> m_interpolator;
 public:
  // Image
  ShaderNode_NonlinearTransform(ShaderNode_Image* node, std::string shader, ShaderBindings& bindings, ShaderNode_Accessor* interpolator = NULL) {
    if(interpolator)
      m_interpolator = interpolator;
    else
      m_interpolator = ShaderNodeCreateDefaultInterpolator();
    previous.reset(node);
    m_shader = shader;
    m_bindings = bindings;
  }
  // LinearTransform
  ShaderNode_NonlinearTransform(ShaderNode_LinearTransform* node, std::string shader, ShaderBindings& bindings, ShaderNode_Accessor* interpolator = NULL) {
    if(interpolator)
      m_interpolator = interpolator;
    else
      m_interpolator = ShaderNodeCreateDefaultInterpolator();
    if(node->m_interpolator->m_quality > m_interpolator->m_quality)
      m_interpolator.reset(node->m_interpolator);
    previous.reset(node);
    m_shader = shader;
    m_bindings = bindings;
  }
  // NonlinearTransform
  ShaderNode_NonlinearTransform(ShaderNode_NonlinearTransform* node, std::string shader, ShaderBindings& bindings, ShaderNode_Accessor* interpolator = NULL) {
    if(interpolator)
      m_interpolator = interpolator;
    else
      m_interpolator = ShaderNodeCreateDefaultInterpolator();
    if(node->m_interpolator->m_quality > m_interpolator->m_quality)
      m_interpolator.reset(node->m_interpolator);
    previous.reset(node);
    m_shader = shader;
    m_bindings = bindings;
  }
  // Base: Accessor & PixelModifier
  ShaderNode_NonlinearTransform(ShaderNode_Base* node, std::string shader, ShaderBindings& bindings, ShaderNode_Accessor* interpolator = NULL) {
    if(interpolator)
      m_interpolator = interpolator;
    else
      m_interpolator = ShaderNodeCreateDefaultInterpolator();
    if(node->m_interpolator->m_quality > m_interpolator->m_quality)
      m_interpolator.reset(node->m_interpolator);
    m_shader = shader;
    m_bindings = bindings;
    previous.reset(node->Rasterize());
  }

};
  */

} } // namespace GPU, namespace vw

#endif
