// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#if(1)

#include <vw/GPU/GPUImage.h>

#include <vw/GPU/Expressions.h>

using std::vector;
using std::string;
using std::stringstream;
using std::list;

namespace vw {
  namespace GPU {

 ShaderNode_Accessor* ShaderNodeCreateDefaultInterpolator() {
   //new ShaderNode_Accessor();
 }

    //typedef GPUImageBase ShaderNode_Image;

// ########################################################################################################
//                                           Base
// ########################################################################################################

    ShaderNode_Image* ShaderNode_Base::Evaluate() {
      ShaderBuilder shaderBuilder;
      RecursiveBind(shaderBuilder);

      if(!shaderBuilder.BindingComplete()) {
        printf("[ShaderNode_BaseL::Evaluate] Shader Binding failed.\n");
        exit(0);
      }
      RecursiveCompose(shaderBuilder);
      if(!shaderBuilder.Complete()) {
        printf("[ShaderNode_BaseL::Evaluate] Shader Composition failed.\n");
        exit(0);
      }
      return shaderBuilder.Evaluate();
    }


// ########################################################################################################
//                                           BaseTransform
// ########################################################################################################

    ShaderNode_Base* ShaderNode_BaseTransform::AppendInterpolateNode() {
      ShaderNode_Accessor* interpolator = m_interpolator.get();
      m_interpolator.release();

      interpolator->m_input_nodes.push_back(this);
      return interpolator;
    }


// ########################################################################################################
//                                           LinearTransform
// ########################################################################################################

  ShaderNode_LinearTransform::ShaderNode_LinearTransform(ShaderNode_Base* node, Matrix<float>& inHomography, ShaderNode_Accessor* interpolator) {
    // COMMON
    if(interpolator)
      m_interpolator.reset(interpolator);
    else
      m_interpolator.reset(ShaderNodeCreateDefaultInterpolator());
    SetSize(node->width(), node->height());
    // LinearTransform: Combine homographies, use previous interpolation mode if its quality is better
    ShaderNode_LinearTransform* linearTransform = dynamic_cast<ShaderNode_LinearTransform*>(node);
    if(linearTransform && linearTransform->m_interpolator->m_quality > m_interpolator->m_quality) {
       m_interpolator.reset(linearTransform->m_interpolator.get());
       linearTransform->m_interpolator.release();
       m_homography = inHomography * linearTransform->m_homography;
       previous.reset(linearTransform->previous.get());
       linearTransform->Disconnect();
       delete node;
       return;
    }
    // NonlinearTransform: use previous interpolation mode if its quality is better
    ShaderNode_NonlinearTransform* nonlinearTransform = dynamic_cast<ShaderNode_NonlinearTransform*>(node);
    if(nonlinearTransform && nonlinearTransform->m_interpolator->m_quality > m_interpolator->m_quality) {
       m_interpolator.reset(nonlinearTransform->m_interpolator.get());
       nonlinearTransform->m_interpolator.release();
       m_homography = inHomography;
       return;
    }
    // Accessor or PixelModifier: prerasterize
    else if(dynamic_cast<ShaderNode_Accessor*>(node) || dynamic_cast<ShaderNode_PixelModifier*>(node)) {
      previous.reset(reinterpret_cast<ShaderNode_Base*>(node->Evaluate()));
      m_homography = inHomography;
    }
  }


  int ShaderNode_LinearTransform::RecursiveBind(ShaderBuilder& shaderBuilder) {
    int imageBindingIndex = previous->RecursiveBind(shaderBuilder);
    // If previous is ImageNode, then add homography and don't bind anything. Otherwise we have to bind an input matrix.
    ShaderNode_Image* imageNode = dynamic_cast<ShaderNode_Image*>(previous.get());
    if(imageNode) {
      //imageNode->SetHomography(m_homography);
    }
    else {
      vector<float> vec(4);
      vec[0] = m_homography(0, 0);
      vec[1] = m_homography(0, 1);
      vec[2] = m_homography(0, 2);
      shaderBindings.vectorBindingStartIndex = shaderBuilder.BindVector(vec);
      vec[0] = m_homography(1, 0);
      vec[1] = m_homography(1, 1);
      vec[2] = m_homography(1, 2);
      shaderBuilder.BindVector(vec);
      vec[0] = m_homography(2, 0);
      vec[1] = m_homography(2, 1);
      vec[2] = m_homography(2, 2);
      shaderBuilder.BindVector(vec);
    }
    return imageBindingIndex;
  }

  void ShaderNode_LinearTransform::RecursiveCompose(ShaderBuilder& shaderBuilder) {
    ShaderNode_Image* imageNode = dynamic_cast<ShaderNode_Image*>(previous.get());
    if(imageNode) {
      previous->RecursiveCompose(shaderBuilder);
    }
    else {
      int functionIndex = shaderBuilder.BindFunction();
      // Add Function
      shaderBuilder.AddFunctionFromFile("NonlinearTransform/homography");
      // function call
      stringstream str;
      str << "transform" << functionIndex;
      shaderBuilder.AppendBody(str.str());
      //
      previous->RecursiveCompose(shaderBuilder);
      shaderBuilder.AppendBody(",");
      // vectors
      int vectorStartIndex = shaderBindings.vectorBindingStartIndex;
      for(int i = 0; i < 3; i++) {
        stringstream str;
        str << "v" << vectorStartIndex + i << ",";
        shaderBuilder.AppendBody(str.str());
      }
      // end
      shaderBuilder.RemoveLastBodyChar();
      shaderBuilder.AppendBody(")");
    }
  }



// ########################################################################################################
//                                           NonlinearTransform
// ########################################################################################################

  ShaderNode_NonlinearTransform::ShaderNode_NonlinearTransform(ShaderNode_Base* node, const string& shader, const ShaderBindings& bindings, ShaderNode_Accessor* interpolator) {
    // COMMON
    if(interpolator)
      m_interpolator.reset(interpolator);
    else
      m_interpolator.reset(ShaderNodeCreateDefaultInterpolator());
    m_shader = shader;
    m_bindings = bindings;
    SetSize(node->width(), node->height());
    // LinearTransform:  use previous interpolation mode if its quality is better
    ShaderNode_LinearTransform* linearTransform = dynamic_cast<ShaderNode_LinearTransform*>(node);
    if(linearTransform && linearTransform->m_interpolator->m_quality > m_interpolator->m_quality) {
       m_interpolator.reset(linearTransform->m_interpolator.get());
       linearTransform->m_interpolator.release();
       previous.reset(node);
       return;
    }
    // NonlinearTransform:  use previous interpolation mode if its quality is better
    ShaderNode_NonlinearTransform* nonlinearTransform = dynamic_cast<ShaderNode_NonlinearTransform*>(node);
    if(nonlinearTransform && nonlinearTransform->m_interpolator->m_quality > m_interpolator->m_quality) {
       m_interpolator.reset(linearTransform->m_interpolator.get());
       linearTransform->m_interpolator.release();
       previous.reset(node);
       return;
    }
    // Accessor or PixelModifier: prerasterize
    if(dynamic_cast<ShaderNode_Accessor*>(node) || dynamic_cast<ShaderNode_PixelModifier*>(node)) {
      previous.reset(reinterpret_cast<ShaderNode_Base*>(node->Evaluate()));
      return;
    }
    // Other
    previous.reset(node);
  }

 int ShaderNode_NonlinearTransform::RecursiveBind(ShaderBuilder& shaderBuilder) {
    int imageBindingIndex = previous->RecursiveBind(shaderBuilder);
    return imageBindingIndex;
 }

 void ShaderNode_NonlinearTransform::RecursiveCompose(ShaderBuilder& shaderBuilder) {
    int functionIndex = shaderBuilder.BindFunction();
    // Add Function
    shaderBuilder.AddFunctionFromFile(m_shader);
    // function call
    stringstream string;
    string << "nonlinear_transform" << functionIndex;
    shaderBuilder.AppendBody(string.str());
    //
    previous->RecursiveCompose(shaderBuilder);
    shaderBuilder.AppendBody(",");
    // floats
    int floatStartIndex = shaderBindings.floatBindingStartIndex;
    for(int i = 0; i < shaderBindings.floats.size(); i++) {
      stringstream string;
      string << "f" << floatStartIndex + i << ",";
      shaderBuilder.AppendBody(string.str());
    }
    // vectors
    int vectorStartIndex = shaderBindings.vectorBindingStartIndex;
    for(int i = 0; i < shaderBindings.vectors.size(); i++) {
      stringstream string;
      string << "v" << vectorStartIndex + i << ",";
      shaderBuilder.AppendBody(string.str());
    }
    // end
    shaderBuilder.RemoveLastBodyChar();
    shaderBuilder.AppendBody(")");

 }


// ########################################################################################################
//                                           Accessor
// ########################################################################################################

  ShaderNode_Accessor::ShaderNode_Accessor(list<ShaderNode_Base*>& input_nodes, const string& shader, const ShaderBindings& bindings, bool forceRasterize, bool isInterpolateNode)   {
    int i = 0;
    for(list<ShaderNode_Base*>::iterator iter = input_nodes.begin(); iter != input_nodes.end(); iter++, i++) {
      // If input is a transform, and this isn't an interpolate node, then we need to call AppendInterpolateNode.
      if(!isInterpolateNode) {
        ShaderNode_BaseTransform* transform = dynamic_cast<ShaderNode_BaseTransform*>(*iter);
        if(transform) {
          *iter = transform->AppendInterpolateNode();
        }
      }

      // Accessor is always prerasterized, transorms are only if forceRasterize is true
      if(dynamic_cast<ShaderNode_Accessor*>(*iter) || dynamic_cast<ShaderNode_Image*>(*iter)
         || (forceRasterize && dynamic_cast<ShaderNode_BaseTransform*>(*iter))) {
        m_input_nodes.push_back(reinterpret_cast<ShaderNode_Base*>((*iter)->Evaluate()));
        //delete (*iter);
      }
      else
        m_input_nodes.push_back(reinterpret_cast<ShaderNode_Base*>((*iter)->Evaluate()));

    }
    m_shader = shader;
    m_bindings = bindings;
  }

  int ShaderNode_Accessor::RecursiveBind(ShaderBuilder& shaderBuilder) {
    for(list<ShaderNode_Base*>::iterator iter = m_input_nodes.begin(); iter != m_input_nodes.end(); iter++) {
      int imageBindingIndex = (*iter)->RecursiveBind(shaderBuilder);
      m_image_indices.push_back(imageBindingIndex);
    }
    shaderBuilder.BindGroup(m_bindings);
    return 0;
  }

  void ShaderNode_Accessor::RecursiveCompose(ShaderBuilder& shaderBuilder) {
    bool needsLastCharRemoved = false;
    int functionIndex = shaderBuilder.BindFunction();
    // Add Function
    shaderBuilder.AddFunctionFromFile(m_shader);
    // function call
    stringstream string;
    string << "pixelmodifier" << functionIndex;
    shaderBuilder.AppendBody(string.str());
    // images
    int i = 0;
    for(list<ShaderNode_Base*>::iterator iter = m_input_nodes.begin(); iter != m_input_nodes.end(); iter++, i++) {
      (*iter)->RecursiveCompose(shaderBuilder);
      shaderBuilder.AppendBody(",");
      needsLastCharRemoved = true;
    }
    // floats
    int floatStartIndex = shaderBindings.floatBindingStartIndex;
    for(int i = 0; i < shaderBindings.floats.size(); i++) {
      stringstream string;
      string << "f" << floatStartIndex + i << ",";
      needsLastCharRemoved = true;
    }
    // vectors
    int vectorStartIndex = shaderBindings.vectorBindingStartIndex;
    for(int i = 0; i < shaderBindings.vectors.size(); i++) {
      stringstream string;
      string << "v" << vectorStartIndex + i << ",";
      needsLastCharRemoved = true;
    }
    // end
    if(needsLastCharRemoved)
      shaderBuilder.RemoveLastBodyChar();
    shaderBuilder.AppendBody(")");
  }


// ########################################################################################################
//                                           PixelModifier
// ########################################################################################################

  ShaderNode_PixelModifier::ShaderNode_PixelModifier(list<ShaderNode_Base*>& input_nodes, const string& shader, const ShaderBindings& bindings) {
    for(list<ShaderNode_Base*>::iterator iter = input_nodes.begin(); iter != input_nodes.end(); iter++) {
      // If input is a transform then we need to call AppendInterpolateNode.
      ShaderNode_BaseTransform* transform = dynamic_cast<ShaderNode_BaseTransform*>(*iter);
      if(transform) {
        *iter = transform->AppendInterpolateNode();
      }
      m_input_nodes.push_back(reinterpret_cast<ShaderNode_Base*>((*iter)->Evaluate()));
    }
    m_input_nodes = input_nodes;
    m_shader = shader;
    m_bindings = bindings;
  }


  int ShaderNode_PixelModifier::RecursiveBind(ShaderBuilder& shaderBuilder) {
    for(list<ShaderNode_Base*>::iterator iter = m_input_nodes.begin(); iter != m_input_nodes.end(); iter++) {
      (*iter)->RecursiveBind(shaderBuilder);
    }
    shaderBuilder.BindGroup(m_bindings);
    return 0;
  }

  void ShaderNode_PixelModifier::RecursiveCompose(ShaderBuilder& shaderBuilder) {
    bool needsLastCharRemoved = false;
    int functionIndex = shaderBuilder.BindFunction();
    // Add Function
    shaderBuilder.AddFunctionFromFile(m_shader);
    // function call
    stringstream string;
    string << "pixelmodifier" << functionIndex;
    shaderBuilder.AppendBody(string.str());
    // images
    for(list<ShaderNode_Base*>::iterator iter = m_input_nodes.begin(); iter != m_input_nodes.end(); iter++) {
      (*iter)->RecursiveCompose(shaderBuilder);
      shaderBuilder.AppendBody(",");
      needsLastCharRemoved = true;
    }
    // floats
    int floatStartIndex = shaderBindings.floatBindingStartIndex;
    for(int i = 0; i < shaderBindings.floats.size(); i++) {
      stringstream string;
      string << "f" << floatStartIndex + i << ",";
      needsLastCharRemoved = true;
    }
    // vectors
    int vectorStartIndex = shaderBindings.vectorBindingStartIndex;
    for(int i = 0; i < shaderBindings.vectors.size(); i++) {
      stringstream string;
      string << "v" << vectorStartIndex + i << ",";
      needsLastCharRemoved = true;
    }
    // end
    if(needsLastCharRemoved)
      shaderBuilder.RemoveLastBodyChar();
    shaderBuilder.AppendBody(")");

  }




  } } // namespace GPU // namespace vw

#endif
