
#ifndef SeperableConvolution_H
#define SeperableConvolution_H

#include "OpenGLManager.h"
#include "TexRef.h"
#include "GPUProgram.h"

TexRef SeperableConvolution(TexRef& image, 
							TexRef& hVector, 
							TexRef& vVector, 
							float factor, 
							float bias);

#endif

