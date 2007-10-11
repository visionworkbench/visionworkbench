/**
 ****************************************************************************
 * <P> XML.c - XML parser test example - wchar_t* version </P>
 *
 * @version     V2.22
 * @author      Frank Vanden Berghen
 *
 * BSD license:
 * Copyright (c) 2002, Frank Vanden Berghen
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Frank Vanden Berghen nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE REGENTS AND CONTRIBUTORS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************
 */
#ifdef WIN32
#define _CRT_SECURE_NO_DEPRECATE
#endif

#include <stdio.h>
#include "xmlParser.h"

int main(int argc, char **argv)
{
    /*********************************************************************************
     *                                                                               *
     *  Example 1: Basic operations to parse and collect data from a XML file        *
     *                                                                               *
     *********************************************************************************/

    // this open and parse the XML file:
    XMLNode xMainNode=XMLNode::openFileHelper(_T("PMMLModel.xml"),_T("PMML"));

    // this prints "RANK For <you>":
    XMLNode xNode=xMainNode.getChildNode(_T("Header"));
    printf("Application Name is: '%S'\n", xNode.getChildNode(_T("Application")).getAttribute(_T("name")));

    // this prints "Hello World!"
    printf("Text inside Header tag is :'%S'\n", xNode.getText());

    // this gets the number of "NumericPredictor" tags:
    xNode=xMainNode.getChildNode(_T("RegressionModel")).getChildNode(_T("RegressionTable"));
    int n=xNode.nChildNode(_T("NumericPredictor"));

    // this prints the "coefficient" value for all the "NumericPredictor" tags:
    int i,myIterator=0;
    for (i=0; i<n; i++)
    printf("coeff %i=%S\n",i+1,xNode.getChildNode(_T("NumericPredictor"),&myIterator).getAttribute(_T("coefficient")));

    // this create a file named "testUnicode.xml" based on the content of the first "Extension" tag of the XML file:
    xMainNode.getChildNode(_T("Extension")).writeToFile(_T("testUnicode.xml"));

    printf("The content of the clear tag is:%S",xMainNode.getChildNode(_T("html_page")).getClear().lpszValue);

    /****************************************************************************************
     *                                                                                      *
     *  Example 2: memory management: when to use the 'stringDup' and the 'free' functions  *
     *                                                                                      *
     ****************************************************************************************/

    // compare these 4 lines ...
    wchar_t *t=stringDup(xMainNode.getAttribute(_T("version")));       // get version number
    xMainNode=XMLNode::emptyXMLNode;                          // free from memory the top of the xml Tree
    printf("PMML Version :%S\n",t);                           // print version number
    free(t);                                                  // free version number

    // ... with the following 3 lines (currently commented, because of error):
    //  t=xMainNode.getAttribute(_T("version"));      // get version number (note that there is no 'stringDup')
    //  xMainNode=XMLNode::emptyXMLNode;              // free from memory the top of the xml Tree AND the version number inside 't' var
    //  printf("PMML Version :%S\n",t);               // since the version number in 't' has been free'd this will not work

    /**********************************************************
     *                                                        *
     *  Example 3: constructing & updating a tree of XMLNode  *
     *                                                        *
     **********************************************************/

    // We create in memory from scratch the following XML structure:
    //  <?xml version="1.0"?>
    //      <body color="#FFFFFF"> Hello "universe". </body>
    // ... and we transform it into a standard C string that is printed on screen.
    xMainNode=XMLNode::createXMLTopNode(_T("xml"),TRUE);
    xMainNode.addAttribute(_T("version"),_T("1.0"));
    xNode=xMainNode.addChild(_T("body"));
    xNode.addText(_T("Hello \"univ\"!"));
    xNode.deleteText();
    xNode.addText(_T("Hello \"universe\"!"));
    xNode.addAttribute(_T("color"),_T("#wrongcolor"));
    xNode.updateAttribute(_T("#FFFFFF"),NULL,_T("color"));

    t=xMainNode.createXMLString(false);
    printf("XMLString created from scratch:\n%S",t);
    free(t);

    // we delete some parts:
    xNode.deleteAttribute(_T("color"));
    t=xMainNode.createXMLString(false);
    printf("\nWith the \"color\" attribute deleted:\n%S\n\n",t);
    free(t);

    /*********************************************************************************************************
     *                                                                                                       *
     *  Example 4: by default, the XML parser is "forgiving" with respect to errors inside XML strings&files *
     *                                                                                                       *
     *********************************************************************************************************/

    // By default, the XML parser is "forgiving":
    // (You can de-activate this behavior: see the header of xmlParser.cpp file)
    wchar_t *t2=(wchar_t*)_T("<a><b>some text</b><b>other text    </a>");
    XMLResults xe;
    xMainNode=XMLNode::parseString(t2,NULL,&xe);
    t=xMainNode.createXMLString(false);
    printf("The following XML: %S\n  ...is parsed as: %S\nwith the following info: '%S'\n",t2,t?t:_T("(null)"),XMLNode::getError(xe.error));
    free(t);

    /*******************************************************
     *                                                     *
     *  Example 5: deleting a part of the tree of XMLNode  *
     *                                                     *
     *******************************************************/

    // this deletes the "<b>other text</b>" subtree part:
    xMainNode.getChildNode(_T("b"),1).deleteNodeContent();

    // To perform the same "delete" as above, we can also do:
    // xNode=xMainNode.getChildNode(_T("a")).getChildNode(_T("b"),1); xNode.deleteNodeContent(); xNode=XMLNode::emptyXMLNode;
    // If you forget the last part of the delete ("xNode=XMLNode::emptyXMLNode"), then the XMLNode will NOT be deleted:
    // As long as there exists a reference to an XMLNode, the smartPointer mechanism prevent the node to be deleted.

    t=xMainNode.createXMLString(false);
    printf("\n...with the wrong node deleted: %S\n",t);
    free(t);

    /************************************************************************************************************
    *                                                                                                          *
    *  Example 5: inserting (and moving) a new XMLNode in the middle of an already existing XMLNode structure  *
    *                                                                                                          *
    ************************************************************************************************************/

    // This creates a XMLNode 'a' that is "<a><b>some text</b><b>other text</b></a>":
    xMainNode=XMLNode::parseString(t2);
    // This creates a XMLNode 'c' that is "<c>hello</c>":
    xNode=XMLNode::parseString(_T("<c>hello</c>"));

    xMainNode.addChild(xNode,0);
    t=xMainNode.createXMLString(false);
    printf("\nWe inserted a new node 'c' as the first tag inside 'a':\n       %S",t);
    free(t);

    xMainNode.addChild(xNode,xMainNode.positionOfChildNode(_T("b"),1));
    t=xMainNode.createXMLString(false);
    printf("\nWe moved the node 'c' at the position of the second 'b' tag:\n       %S\n",t);
    free(t);

    /*******************************************
     *                                         *
     *  Example 6: base 64 encoding/decoding   *
     *                                         *
     *******************************************/

    unsigned char *originalBinaryData=(unsigned char *)"this is binary data.";
    XMLParserBase64Tool b64;
    t=b64.encode(originalBinaryData,21);
    printf(
      "\nTo be able to include any binary data into an xml file, some Base64 conversion"
      "\nfunctions (binary data <--> ascii text) are provided:\n"
      "  original binary data   : %s\n"
      "  encoded as text        : %S\n",originalBinaryData,t);
    printf("  decoded as binary again: %s\n",b64.decode(t));

    /******************************************************
     *                                                    *
     *  Example 7: usage of the "getParentNode()" method  *
     *                                                    *
     ******************************************************/

    printf("\nTwo examples of usage of the \"getParentNode()\" method:\n");
    // let's consider these 2 examples (each example on a separate line):
    xMainNode=XMLNode::parseString(t2);     xNode=xMainNode.getChildNode();         xNode=xNode.getParentNode(); t=(wchar_t*)    xNode.getName(); printf(" Ex1: Name of top node; '%S'\n",t?t:_T("null"));
    xMainNode=XMLNode::parseString(t2); xMainNode=xMainNode.getChildNode(); xMainNode=xMainNode.getParentNode(); t=(wchar_t*)xMainNode.getName(); printf(" Ex2: Name of top node; '%S'\n",t?t:_T("null"));
    // In these two examples, I create a tree of XMLNode based on the string
    // "<a><b>some text</b><b>other text</b></a>". After parsing this string
    // I get a XMLNode that represents the <a> tag. Thereafter I "go down" one
    // level, using getChildNode: I now have a XMLNode that represents the <b> tag.
    // Thereafter I "go up" one level, using getParentNode(): I now have once again
    // a XMLNode that represents the <a> tag. Thereafter, I print the name ('a') of
    // this last XMLNode. The first example is working as intended (it prints 'a'
    // on the screen). However, the second example prints "null" because when we
    // do xMainNode=xMainNode.getChildNode()" we lost all references to the
    // top node and thus it's automatically "garbage collected" (free memory).

    return 0;
}
