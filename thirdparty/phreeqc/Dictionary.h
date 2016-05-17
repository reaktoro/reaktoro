#if !defined(DICTIONARY_H_INCLUDED)
#define DICTIONARY_H_INCLUDED
#include <iostream>
#include <map>
#include <vector>
#include <sstream>
class Phreeqc;
class Dictionary
{
public:
	Dictionary(void);
	Dictionary(std::string & words_string);
	~Dictionary(void);
	int Find(std::string str);
	int MapSize() {return (int) this->dictionary_map.size();}
	int OssSize() {return (int) this->dictionary_oss.str().size();}
	std::ostringstream &GetDictionaryOss() {return this->dictionary_oss;}
	std::vector<std::string> &GetWords() {return this->words;}

protected:
	std::map<std::string, int> dictionary_map;
	std::vector<std::string> words;
	std::ostringstream dictionary_oss;

};

#endif // !defined(DICTIONARY_H_INCLUDED)

