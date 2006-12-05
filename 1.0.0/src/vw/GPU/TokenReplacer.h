
#ifndef TokenReplacer_H
#define TokenReplacer_H

#include <map>
#include <string>
using namespace std;


namespace vw { namespace GPU {


class TokenReplacer {
	map<string, string> stringMap;
public:
	void AddVariable(string& variableName, string& variableValue);
	void AddVariable(const char* variableName, const char*  variableValue);
	void Replace(string& inText, string& outText);
// INLINE
	void Clear() {
		stringMap.clear();
	};
};


} } // namespaces GPU, vw


#endif
