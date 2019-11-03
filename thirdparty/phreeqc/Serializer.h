#if !defined(SERIALIZER_H_INCLUDED)
#define SERIALIZER_H_INCLUDED
#include <iostream>
#include <map>
#include <vector>
#include <sstream>
#include "PHRQ_base.h"
#include "Dictionary.h"
class Phreeqc;
class Serializer : public PHRQ_base
{
public:
	Serializer(PHRQ_io *io = NULL);
	~Serializer(void);

	enum PACK_TYPE
	{
		PT_SOLUTION     = 0,
		PT_EXCHANGE     = 1,
		PT_GASPHASE     = 2,
		PT_KINETICS     = 3,
		PT_PPASSEMBLAGE = 4,
		PT_SSASSEMBLAGE = 5,
		PT_SURFACES     = 6,
		PT_TEMPERATURE  = 7,
		PT_PRESSURE     = 8
	};
	bool Serialize(Phreeqc &phreeqc_ptr, int start, int end, bool include_t, bool include_p, PHRQ_io *io = NULL);
	bool Deserialize(Phreeqc &phreeqc_ptr, Dictionary &dictionary, std::vector<int> &ints, std::vector<double> &doubles);
	Dictionary &GetDictionary(void) {return this->dictionary;}
	std::vector<int> &GetInts(void) {return this->ints;}
	std::vector<double> &GetDoubles(void) {return this->doubles;}
	//std::string &GetWordsString(void) {return this->words_string;}


protected:
	std::vector<int> ints;
	std::vector<double> doubles;
	//std::string words_string;
	Dictionary dictionary;
};
#endif // !defined(SERIALIZER_H_INCLUDED)
