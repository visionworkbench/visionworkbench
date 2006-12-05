/*
 *  SeperableConvolution.cpp
 *  VWGPU
 *
 *  Created by Ian Saxton on 5/17/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "SeperableConvolution.h"



//########################################################################
//#    RowConvolution <int $SIZE>                             
//#    Channel Variants: [r, rgba]
//#	   
//########################################################################





//########################################################################
//#    SeperableConvolution                         
//########################################################################


#include "OpenGLManager.h"
#include "TexRef.h"
#include "SeperableConvolution.h"
#include "Timer.h"


TexRef SeperableConvolution(TexRef& image, 
							TexRef& hVector, 
							TexRef& vVector, 
							float factor, 
							float bias) 
{
// Static
	static vector<int> fAttributes(2);
	static GPUProgramSet_GLSL programSet_Rows;
	static GPUProgramSet_GLSL programSet_Columns;
	static bool needsInit = true;
	if(needsInit) {
		programSet_Rows.SetBasePaths("", "ConvolutionRows");
		programSet_Columns.SetBasePaths("", "ConvolutionColumns");
		needsInit = false;
	}
	vector<int> emptyVector;
// GLState
	glEnable(GL_TEXTURE_RECTANGLE_ARB); 
	glPolygonMode(GL_FRONT,GL_FILL);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, g_framebuffer);
// Variables
	int width = image.Width();
	int height = image.Height();
	Tex_Format format = image.Format();
	int nComps;
	if(format != TEX_RGBA && format != TEX_RGB) 
		nComps = 4;
	else 
		nComps = 4;
// Viewport
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0,width,0.0, height);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glViewport(0,0,width, height);
// Temp Texture	
	TexRef temp(width, height, image.Format(), image.Type());
// *********** STAGE 1 - Row Convolution *******************
	//temp.BindToFramebufferAsColor0();
// Program Install
	fAttributes[0] = nComps;
	fAttributes[1] = hVector.Width();
	GPUProgram* program_row = programSet_Rows.GetProgram(emptyVector, fAttributes, true);
	program_row->Install();
// OUTPUT
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, g_framebuffer);
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, temp.Target(), temp.Name(), 0);	
// INPUT
	program_row->SetUniformTexture("image", 0, image);
	program_row->SetUniformTexture("vector", 1, hVector);
	
	program_row->SetUniform1f("factor", factor);
	program_row->SetUniform1f("bias", bias);
	program_row->SetUniform1f("halfSize", (hVector.Width() / 2) - 1);	
// DRAW
	float imageX = image.X(width);
	float imageY = image.Y(height);
	
	glBegin(GL_QUADS);							  
	glTexCoord2f(0.0, 0.0); 
	glVertex2f(0.0, 0.0);
	
	glTexCoord2f(imageX, 0.0); 
	glVertex2f(width, 0.0);
	
	glTexCoord2f(imageX, imageY);
	glVertex2f(width, height);
	
	glTexCoord2f(0.0, imageY); 
	glVertex2f(0.0, height);
	glEnd();

// *********** STAGE 2 - Column Convolution *******************
	TexRef temp2(width, height, image.Format(), image.Type());

	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, temp2.Target(), temp2.Name(), 0);	
// Program Set Up 
	fAttributes[0] = nComps;
	fAttributes[1] = vVector.Height();
	GPUProgram* program_column = programSet_Columns.GetProgram(emptyVector, fAttributes, false);
	program_column->Install();

	program_column->SetUniformTexture("image", 0, temp);
	program_column->SetUniformTexture("vector", 1, vVector);
	program_column->SetUniform1f("factor", factor);
	program_column->SetUniform1f("bias", bias);
	program_column->SetUniform1f("halfSize", (vVector.Height() / 2) - 1);

// Draw
	glBegin(GL_QUADS);							  
	glTexCoord2f(0.0, 0.0); 
	glVertex2f(0.0, 0.0);
	
	glTexCoord2f(imageX, 0.0); 
	glVertex2f(width, 0.0);
	
	glTexCoord2f(imageX, imageY);
	glVertex2f(width, height);
	
	glTexCoord2f(0.0, imageY); 
	glVertex2f(0.0, height);
	glEnd();
// Restore State
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
	
	return temp2;
}

