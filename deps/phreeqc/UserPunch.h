#if !defined(USERPUNCH_H_INCLUDED)
#define USERPUNCH_H_INCLUDED
#include <vector>
#include <string>				// std::string
#include "NumKeyword.h"
class Phreeqc;
class UserPunch:public cxxNumKeyword
{
public:
	UserPunch(int n=1, PHRQ_io *io=NULL);
	~UserPunch(void);

	// headings
	//
	std::vector <std::string> &Get_headings() {return this->headings;}
	const std::vector <std::string> &Get_headings()const {return this->headings;}

	void Set_headings(std::vector <std::string> & h) {this->headings = h;}

	// phreeqc pointer
	//
	Phreeqc * Get_PhreeqcPtr() {return this->PhreeqcPtr;}
	const Phreeqc * Get_PhreeqcPtr()const {return this->PhreeqcPtr;}

	void Set_PhreeqcPtr(Phreeqc * p) {this->PhreeqcPtr = p;}

	// rate
	//
	struct rate * Get_rate() {return this->rate;}
	const struct rate * Get_rate()const {return this->rate;}

	void Set_rate(struct rate * r) {this->rate = r;}

protected:
	std::vector <std::string> headings;
	struct rate *             rate;
	Phreeqc *                 PhreeqcPtr;
};
#endif // !defined(USERPUNCH_H_INCLUDED)