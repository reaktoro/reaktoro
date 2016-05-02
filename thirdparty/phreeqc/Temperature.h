#if !defined(TEMPERATURE_H_INCLUDED)
#define TEMPERATURE_H_INCLUDED

#include <cassert>				// assert
#include <map>					// std::map
#include <string>				// std::string
#include <list>					// std::list
#include <vector>				// std::vector

#include "NumKeyword.h"

class cxxTemperature:public cxxNumKeyword
{

  public:
	cxxTemperature(PHRQ_io *io=NULL);
	~cxxTemperature();

	//void dump_xml(std::ostream& os, unsigned int indent = 0)const;

	void dump_raw(std::ostream & s_oss, unsigned int indent, int *n_out=NULL) const;

	void read_raw(CParser & parser, bool check = false);
	int read(CParser & parser);
	LDBLE Temperature_for_step(int step_number);
	std::vector<LDBLE> & Get_temps(void) {return temps;}
	const std::vector<LDBLE> & Get_temps(void)const {return temps;}
	int Get_countTemps(void) const;
	void Set_countTemps(int i) {countTemps = i;}
	bool Get_equalIncrements(void) const {return equalIncrements;}
	void Set_equalIncrements(bool tf) {equalIncrements = tf;}
	void Serialize(Dictionary & dictionary, std::vector < int >&ints, std::vector < double >&doubles);
	void Deserialize(Dictionary & dictionary, std::vector < int >&ints, std::vector < double >&doubles, int &ii, int &dd);
	
protected:
	std::vector < LDBLE >temps;
	int countTemps;
	bool equalIncrements;
	const static std::vector < std::string > vopts;
};

#endif // !defined(TEMPERATURE_H_INCLUDED)
