
#ifndef SeparableConvolution_H
#define SeparableConvolution_H

#include "OpenGLManager.h"
#include "TexRef.h"
#include "GPUProgram.h"

TexRef SeparableConvolution(TexRef& image, 
							TexRef& hVector, 
							TexRef& vVector, 
							float factor, 
							float bias);

#endif

