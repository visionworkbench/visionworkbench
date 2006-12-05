
#include "TokenReplacer.h"
#include <vector>
#include <list>

namespace vw { namespace GPU {


int TokenReplacerDemo(int argc, const char** argv) {
	string testString("This is a test. Size = $SIZE_TYPE; Iterations = $ITERATIONS. Done.");
	string replacedString;
	
	TokenReplacer svr;
	svr.AddVariable("SIZE_TYPE", "5");
	svr.AddVariable("ITERATIONS", "100");
	
	svr.Replace(testString, replacedString);
	printf("ORIGINAL: %s\n", testString.c_str());
	printf("REPLACED: %s\n", replacedString.c_str());
	
}

void TokenReplacer::AddVariable(string& variableName, string& variableValue) {
	stringMap[variableName] = variableValue;
}

void TokenReplacer::AddVariable(const char* variableName, const char*  variableValue) {
	stringMap[string(variableName)] = string(variableValue);
}


void TokenReplacer::Replace(string& inText, string& outText) {
	int inLength = inText.size();
	outText.resize(0);
	outText.reserve((int) (inLength * 1.5));
	string tempToken;
	int i=0;
	
	while(i < inLength) {
		char cur = inText[i];
		if(cur == '$') {
			tempToken.resize(0);
			tempToken.reserve(32);
			i++;
			while(i < inLength) {
				cur = inText[i];
				if((cur > 47 && cur < 58) || (cur > 64 && cur < 91) || (cur > 96 && cur < 123) || cur == 95) {
					tempToken.push_back(cur);
					i++;
				}
				else {
					break;
				}
			}
			if(tempToken.size()) {
				map<string, string>::iterator iter = stringMap.find(tempToken);
				if(iter != stringMap.end()) {
					outText.append((*iter).second);  // OK to push string?
				}
			}
		}
		else {
			outText.push_back(cur);
			i++;
		}
	}
}


} } // namespaces GPU, vw
