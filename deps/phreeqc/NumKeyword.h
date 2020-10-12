#if !defined(NUMKEYWORD_H_INCLUDED)
#define NUMKEYWORD_H_INCLUDED

#include <ostream>				// std::ostream
#include <string>				// std::string
#include "PHRQ_base.h"
class CParser;

class cxxNumKeyword: public PHRQ_base
{
  public:
	cxxNumKeyword(PHRQ_io *io=NULL);
	virtual ~ cxxNumKeyword();
	const std::string &Get_description() const
	{
		return this->description;
	}
	void Set_description(std::string str)
	{
		this->description = str;
	}
	void Set_description(const char *str)
	{
		if (str != NULL)
			this->description = str;
	}

	int Get_n_user() const
	{
		return this->n_user;
	}
	void Set_n_user(int user)
	{
		this->n_user = user;
	}

	int Get_n_user_end() const
	{
		return this->n_user_end;
	}
	void Set_n_user_end(int user_end)
	{
		this->n_user_end = user_end;
	}
	void Set_n_user_both(int user_end)
	{
		this->n_user = this->n_user_end = user_end;
	}

	bool operator<(const cxxNumKeyword & key) const
	{
		return (this->n_user < key.n_user);
	}

	virtual void dump_xml(std::ostream & os, unsigned int indent = 0) const;

	void read_number_description(CParser & parser);
	void read_number_description(const std::string & line_in);
  protected:
	int n_user;
	int n_user_end;
	std::string description;

  protected:
	void read_number_description(std::istream & is);
};
#endif // !defined(NUMKEYWORD_H_INCLUDED)
