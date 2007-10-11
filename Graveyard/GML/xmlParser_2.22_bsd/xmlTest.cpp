/**
 ****************************************************************************
 * <P> XML.c - XML parser test example - char* version </P>
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

void myfree(void *t); // {free(t);}

int main(int argc, char **argv)
{
    /*********************************************************************************
     *                                                                               *
     *  Example 1: Basic operations to parse and collect data from a XML file        *
     *                                                                               *
     *********************************************************************************/

    // this open and parse the XML file:
    XMLNode xMainNode=XMLNode::openFileHelper("PMMLModel.xml","PMML");

    // this prints "RANK For <you>":
    XMLNode xNode=xMainNode.getChildNode("Header");
    printf("Application Name is: '%s' (note that &lt; has been replaced by '<')\n", xNode.getChildNode("Application").getAttribute("name"));

    // this prints "Hello World!"
    printf("Text inside Header tag is :'%s'\n", xNode.getText());

    // this gets the number of "NumericPredictor" tags:
    xNode=xMainNode.getChildNode("RegressionModel").getChildNode("RegressionTable");
    int n=xNode.nChildNode("NumericPredictor");

    // this prints the "coefficient" value for all the "NumericPredictor" tags:
    int i,myIterator=0;
    for (i=0; i<n; i++)
        printf("coeff %i=%s\n",i+1,xNode.getChildNode("NumericPredictor",&myIterator).getAttribute("coefficient"));

    // this create a file named "test.xml" based on the content of the first "Extension" tag of the XML file:
    xMainNode.getChildNode(_T("Extension")).writeToFile("test.xml","ISO-8859-1");

    printf("The content of the clear tag is:%s",xMainNode.getChildNode("html_page").getClear().lpszValue);

    /****************************************************************************************
     *                                                                                      *
     *  Example 2: memory management: when to use the 'stringDup' and the 'free' functions  *
     *                                                                                      *
     ****************************************************************************************/

    // compare these 4 lines ...
    char *t=stringDup(xMainNode.getAttribute("version")); // get version number
    xMainNode=XMLNode::emptyXMLNode;                      // free from memory the top of the xml Tree
    printf("PMML Version :%s\n",t);                       // print version number
    myfree(t);                                              // free version number

    // ... with the following 3 lines (currently commented, because of error):
    //  t=xMainNode.getAttribute("version");      // get version number (note that there is no 'stringDup')
    //  xMainNode=XMLNode::emptyXMLNode;          // free from memory the top of the xml Tree AND the version number inside 't' var
    //  printf("PMML Version :%s\n",t);           // since the version number in 't' has been free'd, this will not work

    /**********************************************************
     *                                                        *
     *  Example 3: constructing & updating a tree of XMLNode  *
     *                                                        *
     **********************************************************/

    // We create in memory from scratch the following XML structure:
    //  <?xml version="1.0"?>
    //      <body color="FFFFFF"> Hello universe. </body>
    // ... and we transform it into a standard C string that is printed on screen.
    xMainNode=XMLNode::createXMLTopNode("xml",TRUE);
    xMainNode.addAttribute("version","1.0");
    xNode=xMainNode.addChild("body");
    xNode.addText("Hello \"univ\"!");
    xNode.deleteText();
    xNode.addText("Hello \"universe\"!");
    xNode.addAttribute("color","#wrongcolor");
    xNode.updateAttribute("#FFFFFF",NULL,"color");

    t=xMainNode.createXMLString(false);
    printf("XMLString created from scratch:\n%s",t);
    myfree(t);

    // we delete some parts:
    xNode.deleteAttribute("color");
    t=xMainNode.createXMLString(false);
    printf("\nWith the \"color\" attribute deleted:\n%s\n\n",t);
    myfree(t);

    /*********************************************************************************************************
     *                                                                                                       *
     *  Example 4: by default, the XML parser is "forgiving" with respect to errors inside XML strings&files *
     *                                                                                                       *
     *********************************************************************************************************/

    // By default, the XML parser is "forgiving":
    // (You can de-activate this behavior: see the header of the xmlParser.cpp file)
    const char *t2="<a><b>some text</b><b>other text    </a>";
    XMLResults xe;
    xMainNode=XMLNode::parseString(t2,NULL,&xe);
    t=xMainNode.createXMLString(false);
    printf("The following XML: %s\n  ...is parsed as: %s\nwith the following info: '%s'\n",t2,t?t:"(null)",XMLNode::getError(xe.error));
    myfree(t);

    /*******************************************************
     *                                                     *
     *  Example 5: deleting a part of the tree of XMLNode  *
     *                                                     *
     *******************************************************/

    // this deletes the "<b>other text</b>" subtree part:
    xMainNode.getChildNode("b",1).deleteNodeContent();

    // To perform the same "delete" as above, we can also do:
    // xNode=xMainNode.getChildNode("a").getChildNode("b",1); xNode.deleteNodeContent(); xNode=XMLNode::emptyXMLNode;
    // If you forget the last part of the delete ("xNode=XMLNode::emptyXMLNode"), then the XMLNode will NOT be deleted:
    // In this case, as long as there exists a reference to the XMLNode, the smartPointer mechanism prevent the node to be deleted.

    // To perform the same "delete" as above, we can also do:
    // xNode=xMainNode.getChildNode("a").getChildNode("b",1); xNode.deleteNodeContent(true);
    // The "true" parameter will force the deletion, even if there still exists some references to the XMLNode.

    t=xMainNode.createXMLString(false);
    printf("\n...with the wrong node deleted: %s\n",t);
    myfree(t);

    /************************************************************************************************************
     *                                                                                                          *
     *  Example 5: inserting (and moving) a new XMLNode in the middle of an already existing XMLNode structure  *
     *                                                                                                          *
     ************************************************************************************************************/

    // This creates a XMLNode 'a' that is "<a><b>some text</b><b>other text</b></a>":
    xMainNode=XMLNode::parseString(t2);
    // This creates a XMLNode 'c' that is "<c>hello</c>":
    xNode=XMLNode::parseString("<c>hello</c>");

    xMainNode.addChild(xNode,0);
    t=xMainNode.createXMLString(false);
    printf("\nWe inserted a new node 'c' as the first tag inside 'a':\n       %s",t);
    myfree(t);

    xMainNode.addChild(xNode,xMainNode.positionOfChildNode("b",1));
    t=xMainNode.createXMLString(false);
    printf("\nWe moved the node 'c' at the position of the second 'b' tag:\n       %s\n",t);
    myfree(t);

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
      "\nfunctions (binary data <--> ascii/utf8 text) are provided:\n"
      "  original binary data   : %s\n"
      "  encoded as text        : %s\n",originalBinaryData,t);
    printf("  decoded as binary again: %s\n",b64.decode(t));

    /***************************************************************
    *                                                              *
    *  Example 7: demonstration of multi-lingual XML file parsing  *
    *                                                              *
    ****************************************************************/

    printf("\nProcessing XML file containing chinese, cyrilic and other extended characters.\n");
    xMainNode=XMLNode::openFileHelper("utf8test.xml");
    xMainNode.writeToFile("outputTestUTF8.xml");
    printf("... resulting multi-lingual file is 'outputTestUTF8.xml'.\n");

    /******************************************************
     *                                                    *
     *  Example 8: usage of the "getParentNode()" method  *
     *                                                    *
     ******************************************************/

    printf("\nTwo examples of usage of the \"getParentNode()\" method:\n");
    // let's consider these 2 examples (each example on a separate line):
    xMainNode=XMLNode::parseString(t2);     xNode=xMainNode.getChildNode();         xNode=xNode.getParentNode(); t=(char*)    xNode.getName(); printf(" Ex1: Name of top node; '%s'\n",t?t:"null");
    xMainNode=XMLNode::parseString(t2); xMainNode=xMainNode.getChildNode(); xMainNode=xMainNode.getParentNode(); t=(char*)xMainNode.getName(); printf(" Ex2: Name of top node; '%s'\n",t?t:"null");
    // In these two examples, I create a tree of XMLNode based on the string
    // "<a><b>some text</b><b>other text</b></a>". After parsing this string
    // I get a XMLNode that represents the <a> tag. Thereafter I "go down" one
    // level, using getChildNode: I now have a XMLNode that represents the <b> tag.
    // Thereafter I "go up" one level, using getParentNode(): I now have once again
    // a XMLNode that represents the <a> tag. Thereafter, I print the name ('a') of
    // this last XMLNode. The first example is working as intended (it prints 'a'
    // on the screen). However, the second example prints "null" because when we
    // did "xMainNode=xMainNode.getChildNode()" we lost all references to the
    // top node and thus it's automatically "garbage collected" (free memory).

    return 0;
}

#ifdef _USE_XMLPARSER_DLL
    // We are using the DLL version of the XMLParser library.
    // NOTE: With visual studio .NET, you can always use the standard "free()" function: You don't
    //       need a special "DLL free" version.
    void myfree(void *t){free_XMLDLL(t);}
#else
    // we are using the normal, classical version of the XMLParser library (directly from C++ sources)
    void myfree(void *t){free(t);}
#endif
