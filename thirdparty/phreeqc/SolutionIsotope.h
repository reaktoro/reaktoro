#if !defined(SOLUTIONISOTOPE_H_INCLUDED)
#define SOLUTIONISOTOPE_H_INCLUDED

#include <ostream>				// std::ostream
#include <string>				// std::string
#include <list>					// std::list
#include "Parser.h"
class Dictionary;

class cxxSolutionIsotope: public PHRQ_base
{
  public:
	cxxSolutionIsotope(PHRQ_io *io=NULL);
	cxxSolutionIsotope(struct isotope *isotope_ptr, PHRQ_io *io=NULL);
	virtual ~cxxSolutionIsotope(void);

	void dump_xml(std::ostream & os, unsigned int indent) const;
	void dump_raw(std::ostream & os, unsigned int indent) const;

	//CParser::STATUS_TYPE read_raw(CParser & parser, std::istream::pos_type next_char);
	void read_raw(CParser & parser, bool check);

	LDBLE Get_isotope_number() const
	{
		return this->isotope_number;
	}
	void Set_isotope_number(LDBLE d)
	{
		this->isotope_number = d;
	}
	const std::string &Get_elt_name() const	{return this->elt_name;}
	void Set_elt_name(const char *cstring)
	{
		if (cstring != NULL)
			this->elt_name = std::string(cstring);
		else
			this->elt_name.clear();
	}

	const std::string &Get_isotope_name() const	{return this->isotope_name;}
	void Set_isotope_name(const char *cstring)
	{
		if (cstring != NULL)
			this->isotope_name = std::string(cstring);
		else
			this->isotope_name.clear();
	}

	LDBLE Get_total() const	{return this->total;}
	void Set_total(LDBLE d)	{this->total = d;}

	LDBLE Get_ratio() const	            {return this->ratio;}
	void Set_ratio(LDBLE t)             {this->ratio = t;}

	LDBLE Get_ratio_uncertainty() const	{return this->ratio_uncertainty;}
	void Set_ratio_uncertainty(LDBLE t) {this->ratio_uncertainty = t;}

	LDBLE Get_x_ratio_uncertainty() const	{return this->x_ratio_uncertainty;}
	void Set_x_ratio_uncertainty(LDBLE t)   {this->x_ratio_uncertainty = t;}

	bool Get_ratio_uncertainty_defined() const  {return this->ratio_uncertainty_defined;}
	void Set_ratio_uncertainty_defined(bool tf) {this->ratio_uncertainty_defined = tf;}
	LDBLE Get_coef() const	{return this->coef;}
	void Set_coef(LDBLE t)   {this->coef = t;}

	bool operator<(const cxxSolutionIsotope & conc) const;
	void add(const cxxSolutionIsotope & isotope_ptr, LDBLE intensive,
			 LDBLE extensive);
	void multiply(LDBLE extensive);
	void Serialize(Dictionary & dictionary, std::vector < int >&ints, std::vector < double >&doubles);
	void Deserialize(Dictionary & dictionary, std::vector < int >&ints, std::vector < double >&doubles, int &ii, int &dd);

  protected:
	LDBLE isotope_number;
	std::string elt_name;
	std::string isotope_name;
	LDBLE total;
	LDBLE ratio;
	LDBLE ratio_uncertainty;
	bool ratio_uncertainty_defined;
 	LDBLE x_ratio_uncertainty;
 	LDBLE coef;					/* coefficient of element in phase */
	const static std::vector < std::string > vopts;
};
#endif // SOLUTIONISOTOPE_H_INCLUDED
