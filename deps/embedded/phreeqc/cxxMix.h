#if !defined(CXXMIX_H_INCLUDED)
#define CXXMIX_H_INCLUDED

#include <cassert>				// assert
#include <map>					// std::map
#include <string>				// std::string
#include <list>					// std::list
#include <vector>				// std::vector
#include "NumKeyword.h"
#include "PHRQ_base.h"
#include "phrqtype.h"

class cxxMix:public cxxNumKeyword
{

  public:
	cxxMix(PHRQ_io *io=NULL);
	 ~cxxMix();

	//void dump_xml(std::ostream& os, unsigned int indent = 0)const;

	void dump_raw(std::ostream & s_oss, unsigned int indent, int *n_out=NULL) const;

	void read_raw(CParser & parser);

	void Add(int n, LDBLE f)
	{
		if (this->mixComps.find(n) != this->mixComps.end())
		{
			mixComps[n] += f;
		}
		else
		{
			mixComps[n] = f;
		}
	}
	void Multiply(LDBLE f)
	{
		for (std::map < int, LDBLE >::iterator it = this->mixComps.begin();
			it != this->mixComps.end(); it++)
		{
			it->second *= f;
		}
	}

	const std::map < int, LDBLE > & Get_mixComps() const
	{
		return mixComps;
	}
	void Vectorize(std::vector<int> &n, std::vector<LDBLE> &f);
  protected:
	friend class cxxStorageBin;
	std::map < int, LDBLE >mixComps;
	const static std::vector < std::string > vopts;

};

#endif // !defined(CXXMIX_H_INCLUDED)
