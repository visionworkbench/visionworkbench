Todd R. Templeton <ttemplet@email.arc.nasa.gov>, 03/08/2007
NOTE: This is a slightly modified version. A 'lib' target has been added to the makefile; this creates the static library lib/libxmlparser.a with header file include/xmlParser.h.




XMLParser v2.22
===============

The library is composed by two files: "xmlParser.cpp" and "xmlParser.h".
These are the only 2 files that you need when using the library inside your own projects.
All the functions of the library are documented inside the comments of the file "xmlParser.h".

To use the unicode windows version, you need to define the "_UNICODE" preprocessor
definition inside your project definition file.

Some small test examples are also given: see the files "xmlTest.cpp" and "xmlTestUnicode.cpp".
The examples are described inside the html file "xmlParser.html".

To build the examples:
- linux/unix: type "make"
- solaris: type "make -f makefile.solaris"
- windows: Visual Studio: double-click on xmlParser.dsw
  (under Visual Studio .NET, the .dsp and .dsw files will be automatically converted to .vcproj and .sln files)

In order to build the examples you need some project files:
- linux/unix: makefile
- solaris: makefile.solaris
- windows: Visual Studio: *.dsp, xmlParser.dsw

WINDOWS:
	Inside Visual C++, the "debug versions" of the memory allocation functions are 
	very slow: Do not forget to compile in "release mode" to get maximum speed. 
	When I have to debug a software that is using the XMLParser Library, it's usually
	a nightmare because the library is sooOOOoooo slow in debug mode. To solve this
	problem, during all the debugging session, I use a very fast DLL version of the 
	XMLParser Library (the DLL is compiled in release mode). Using the DLL version of 
	the XMLParser Library allows me to have lightening XML parsing speed even in debug!
	Other than that, the DLL version is useless: In the release version of my tool,
	I always use the normal, ".cpp"-based, XMLParser Library.
LINUX:
	The speed of the debug version of the XMLParser library is tolerable.

Change Log
----------

* V1.00: February 20, 2002: initial version.
* V1.01: February 13, 2005: first bug-free release.
* V1.02: March 6, 2005: 2 minor changes:
   o "parseString" function declaration changed to allow easy parsing from memory buffer
   o Minor changes to allow easy compilation under old GCC under QNX
* V1.03: April 2,2005: 3 minors changes:
   o When parsing from a user-supplied memory buffer, the library was previously modifying the content of the memory buffer. This is not the case anymore
   o Non-unicode Windows version: You can now work with unicode XML files: They are converted to ANSI charset before being processed
   o Added Visual Studio 6.0 project files
* V1.04: May 16, 2005: 3 minors changes, 1 bug fix:
   o FIX: When creating an xml string with no formatting, the formatting did not work always (due to an un-initialized variable)
   o Improved parsing speed (try increasing the constant "memoryIncrease" if you need more speed)
   o Minor changes to allow easy compilation under MSYS/MINGW under Windows
   o Added more character entities
* V1.05: May 31, 2005: 2 minors changes:
   o Changed some "char *" to "const char *"
   o Improved robustness against badly formed xml strings
* V1.06: July 11, 2005: 1 change, 1 bug fix: 
   o FIX: Some character entities were not previously correctly processed.
   o Major speed improvement. The library is now at least 10 times faster. (Try increasing the constant "memoryIncrease" if you need more speed)
   o moved the log file out of the HTML file
* V1.07: July 25, 2005: 1 change
   o Added a pre-compiler directive named "APPROXIMATE_PARSING". See header of xmlParser.cpp for more info.
* V1.08: September 8,2005: 1 bug fix: 
   o FIX: on special cases, non-matching quotes were causing malfunction
* V1.09: November 22, 2005: 1 addition
   o Added some new functions to be able to easily create a XML structure in memory
* V1.10: December 29, 2005: 2 minor change.
   o Changed some formatting when rendering a XML tree to a string 
   o added the STRICT_PARSING option
* V1.11: December 31, 2005: 1 bug fix: 
   o FIX: reduced memory consumption.
* V1.12: January 4, 2006: 1 addition. 
   o added the function "removeNodeContent" to delete a subtree
* V1.13: February 25, 2006: 1 addition. 
   o added a primitive UNICODE support under linux (thanks to Sunny Bains)
* V1.14: April 24, 2006: 1 bug fix: 
   o FIX: memory allocation errors when the XML tree is created from scratch using "addChild" method.
* V1.15: April 28, 2006: 2 additions
   o added some methods to delete attributes,clearTags and textFields from an XMLNode tree.
   o added the "addChild(XMLNode x)" method
* V1.16: May 17, 2006: 1 bug fix: 
   o FIX: memory allocation errors under linux
* V1.17: May 28, 2006: 1 bug fix, 2 additions:
   o FIX: character entities not always processed inside text block
   o position of the eXMLErrorMissingEndTag error is computed
   o added the eXMLErrorUnknownEscapeSequence
* V1.18: June 8, 2006: 1 bug fix, minors changes
   o FIX: the 'eXMLErrorFirstTagNotFound' error was not reported.
   o changed license to BSD - added some examples of usage.
* V1.19: July 4, 2006: 3 addition.
   o added automatic convertion from/to UNICODE/ANSI in linux (this was already done in windows)
   o added getChildNodeWithAttribute()
   o added support for SOLARIS unicode (Thanks to Joseph Vijay!).
   o added support for 32 bit unicode (so that the library works on Redhat Enterprise v4 EMT64).
* V1.20: July 22, 2006: 13 additions.
   o added 9 "update" functions (like updateAttribute(LPCTSTR lpszNewValue, LPCTSTR lpszNewName=NULL,LPCTSTR lpszOldName);)
   o added 4 functions that allows you to include any binary data (images, sounds,...) into an XML file or string using "Base64 encoding".
* V2.01: July 24, 2006: 1 major change, 2 minor change, 3 additions
   o added extended support for strict UTF-8 character encoding (The characters in UTF-8 have different lengths (from 1 byte to 4 bytes)).
     Chinese, cyrilic and other extended characters are now processed successfully.
   o new character entities: This is now equivalent: '&#x41;' or '&#65;' or 'A' 
     (The ascci code of 'A' is 65 in decimal and 0x41 in hexadecimal).
   o added a function that try to guess if the encoding is UTF-8.
   o the code has been modified in order to allow easy inclusion of new entities and new clearTags (minor change).
   o the "updateAttribute" function is now adding a new attribute if the one to update is missing.
     (same behavior for "updateText" and "updateClear").
   o no more "stringDup" required for functions like "addText", "addAttribute",... 
     The old behavior is still accessible through functions like "addText_WOSD", "addAttribute_WOSD",...
	 ("_WSOD" stands for "WithOut StringDup").
     This change greatly simplifies the user's code (major update). 
     Unfortunately, old user's code must be updated to work with the new version.
     Fortunately, all the user's code used to READ the content of an XML file is left unchanged:
     Only the "creation of XML" and the "update of XML" user's code require a little updating work.
* V2.02: July 25, 2006: 1 minor change
   o changed the function "createXMLTopNode()" to "createXMLTopNode(LPCTSTR lpszName, int isDeclaration=FALSE);".
* V2.03: July 28, 2006: 1 minor change
   o changed LPTSTR to XMLSTR to avoid name-clash with the definitions in <winnt.h>
* V2.04: August 6, 2006: 1 addition
   o added one heuristic inside the function "guessUTF8ParsingParameterValue".
* V2.05: August 15, 2006: 1 addition
   o now displays the error message inside the method "openFileHelper" in a MessageBox window (WIN32 only).
* V2.06: August 16, 2006: 2 additions
   o added the method XMLNode::writeToFile to make it easier to write the content of an XMLNode to a file.
   o added support for Byte-order marks (or BOM).
* V2.07: August 22, 2006: 1 additions
   o added _XMLUNICODE preprocessor variable to make it easy to force the library into either utf16-mode or utf8-mode.
* V2.08: August 22, 2006: 1 bug fix
   o inside the tag content, the ">" and "/>" strings are not reported as errors anymore.
* V2.09: August 31, 2006: 1 bug fix
   o the character entities of type &#x04B; were not working properly (thanks to José Carlos Medeiros for notifying me!).
* V2.10: September 21, 2006: 1 bug fix
   o two consecutive calls to the deleteNodeContent() function on the same node has now no effect (as it should be).
     (Thanks to Hartmut Lemmel for notifying me!)
   o improved compatibility to Borland C++
* V2.11: October 24, 2006: 3 additions, 1 bug fix.
   o added the function getParentNode(). Thanks to Jakub Siudzinski for notifying me a good way to do it easily.
   o added one parameter to the deleteNodeContent() function to force the deletion of the underlying XMLNode tree. 
     This will release all the memory occupied by the XMLNode tree even if there still exist references to some part
     of the tree.
   o changed the usage of the base64Encode() function to reduce the number of malloc's (speed increase).
   o FIX: when parsing an XML string, if the TOP tag has no closing tag, the library now correctly 
     reports "eXMLErrorMissingEndTag".
* V2.12: October 25, 2006: 2 additions
   o refactoring of the Base64 functions to make things easier to use
   o added the _XMLPARSER_NO_MESSAGEBOX_ preprocessor variable (see header of xmlParser.cpp for explanation)
* V2.13: October 31, 2006: 1 minor change, 1 bug fix
   o changed the signature of _strnicmp to allow easy compilation under linux .
   o FIX: size of buffer for the convertion from ascii/utf8 to utf-16 was incorrect.
* V2.14: November 13, 2006: 1 minor change, 1 bug fix
   o changed the parseFile,openFileHelper,writeToFile functions so that the filename parameter is widechar when UNICODE=1
   o fixed a bug in openFileHelper when sizeof(wchar_t)=4
* V2.15: December 22, 2006: 2 additions
   o added the parameter 'pos' to the addChild,addText,addClear methods to allow insertion of new components anywhere in an already existing XMLNode structure
   o added 'postionOf*' methods.
* V2.16: December 27, 2006: 1 minor change
   o removed the un-necessary method "firstPosition()" & some code re-structuration.
* V2.17: January 9, 2007: 1 addition, 1 minor change
   o added the preprocessor variable "XML_NO_WIDE_CHAR" to allow easy compilation on exotic compilers
   o added the "const" method qualifier to some methods.
* V2.18: January 15, 2007: 1 bug fix
   o FIX: addChild(XMLNode x,int pos) was sometime inserting at the wrong position when pos!=-1
* V2.19: January 30, 2007: 1 bug fix, 3 additions.
   o FIX: Unknown Character Entities are now always reported correctly. Thanks to Vincent Vanhoucke.
   o The XML specification indicates that no white spaces should be lost when parsing the file. This is now possible setting the
     new global parameter "dropWhiteSpaces" to false.
   o The library now works under Windows CE 4.2, Windows Mobile (PPC) 2003(5) (xscale) and Mac OS X Tiger. Thanks to Zdenek Nemec.
   o The "<!DOCTYPE" tag is now always handled properly. Thanks to Zdenek Nemec.
* V2.20: February 17, 2007: 1 addition
   o added a Visual Studio projet file to build a DLL version of the library.
     Under Windows, when I have to debug a software that is using the XMLParser Library, 
     it's usually a nightmare because the library is sooOOOoooo slow in debug mode. To 
     solve this problem, during all the debugging session, I use a very fast DLL version 
     of the XMLParser Library (the DLL is compiled in release mode). Using the DLL version 
     of the XMLParser Library allows me to have lightening XML parsing speed, even in 
     debug mode! Other than that, the DLL version is useless: In the release version 
     of my tool, I always use the normal, ".cpp"-based, XMLParser Library.
* V2.21: Mars 1, 2007: 1 minor change, 1 bug fix
   o changed to a better algorithm to handle the "<!DOCTYPE" tag.
   o FIX: under SPARC processor, the heuristic that is used to check if the xml file is wchar_t* was generating a BUS error.
* V2.22: Mars 6, 2007: 1 bug fix
   o FIX: the 'tag' parameter was not always working properly in parseString
