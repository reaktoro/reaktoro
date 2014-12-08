#if !defined(REACTION_H_INCLUDED)
#define REACTION_H_INCLUDED

#include <cassert>				// assert
#include <map>					// std::map
#include <string>				// std::string
#include <list>					// std::list
#include <vector>				// std::vector

#include "NumKeyword.h"
#include "NameDouble.h"

class cxxReaction:public cxxNumKeyword
{

  public:
	cxxReaction(PHRQ_io *io = NULL);
	 ~cxxReaction();

	//void dump_xml(std::ostream& os, unsigned int indent = 0)const;

	void dump_raw(std::ostream & s_oss, unsigned int indent, int *n_out=NULL) const;

	void read_raw(CParser & parser, bool check=true);
	const cxxNameDouble &Get_elementList(void) const {return this->elementList;}
	void Set_elementList(cxxNameDouble nd) {this->elementList = nd;}
	cxxNameDouble &Get_reactantList(void) {return this->reactantList;}
	const cxxNameDouble &Get_reactantList(void)const {return this->reactantList;}
	std::vector < LDBLE > &Get_steps(void) {return this->steps;}
	const std::vector < LDBLE > &Get_steps(void)const {return this->steps;}
	void Set_steps(std::vector<LDBLE> &v) {steps = v;}
	int Get_reaction_steps(void) const;
	int Get_countSteps(void) const {return this->countSteps;}
	void Set_countSteps(int i) {countSteps = i;}
	bool Get_equalIncrements(void) const {return this->equalIncrements;}
	void Set_equalIncrements(bool tf) {equalIncrements = tf;}
	const std::string &Get_units(void) const {return this->units;}

	void Set_units(const char * s)
	{
		if (s != NULL)
			this->units = std::string(s);
		else
			this->units.clear();
	}

protected:
	cxxNameDouble reactantList;
	cxxNameDouble elementList;
	std::vector < LDBLE >steps;
	int countSteps;
	bool equalIncrements;
	std::string units;
	const static std::vector < std::string > vopts;
};

#endif // !defined(REACTION_H_INCLUDED)
