#if !defined(STORAGEBINLIST_H_INCLUDED)
#define STORAGEBINLIST_H_INCLUDED
#include <set>                  // std::set
#include <string>               // std::string
#include <list>                 // std::list
#include <vector>               // std::vector
#include "PHRQ_base.h"
class CParser;

class StorageBinListItem
{
public:
	StorageBinListItem(void);
	StorageBinListItem(CParser & parser);
	~StorageBinListItem(void);
	void Set_defined(bool tf) { this->defined = tf; };
	bool Get_defined(void)const { return(this->defined); };
	void Augment(std::string token);
	void Augment(int i);
	std::set < int > &Get_numbers(void) { return(this->numbers); };
	const std::set < int > &Get_numbers(void)const { return(this->numbers); };
	void Clear(void) { this->numbers.clear(); };
protected:
	std::set < int > numbers;
	bool defined;
};
class StorageBinList: public PHRQ_base
{
public:
	StorageBinList(PHRQ_io *io=NULL);
	StorageBinList(CParser & parser, PHRQ_io *io=NULL);
	virtual ~StorageBinList(void);
	bool Read(CParser & parser);
	void SetAll(bool tf);
	void TransferAll(StorageBinListItem &source);
	std::set<StorageBinListItem *> GetAllItems(void);

	StorageBinListItem & Get_solution(void)      { return(this->solution); };
	StorageBinListItem & Get_pp_assemblage(void) { return(this->pp_assemblage); };
	StorageBinListItem & Get_exchange(void)      { return(this->exchange); };
	StorageBinListItem & Get_surface(void)       { return(this->surface); };
	StorageBinListItem & Get_ss_assemblage(void) { return(this->ss_assemblage); };
	StorageBinListItem & Get_gas_phase(void)     { return(this->gas_phase); };
	StorageBinListItem & Get_kinetics(void)      { return(this->kinetics); };
	StorageBinListItem & Get_mix(void)           { return(this->mix); };
	StorageBinListItem & Get_reaction(void)      { return(this->reaction); };
	StorageBinListItem & Get_temperature(void)   { return(this->temperature); };
	StorageBinListItem & Get_pressure(void)      { return(this->pressure); };
	StorageBinListItem & Get_cell(void)          { return(this->cell); };

	const StorageBinListItem & Get_solution(void)const      { return(this->solution); };
	const StorageBinListItem & Get_pp_assemblage(void)const { return(this->pp_assemblage); };
	const StorageBinListItem & Get_exchange(void)const      { return(this->exchange); };
	const StorageBinListItem & Get_surface(void)const       { return(this->surface); };
	const StorageBinListItem & Get_ss_assemblage(void)const { return(this->ss_assemblage); };
	const StorageBinListItem & Get_gas_phase(void)const     { return(this->gas_phase); };
	const StorageBinListItem & Get_kinetics(void)const      { return(this->kinetics); };
	const StorageBinListItem & Get_mix(void)const           { return(this->mix); };
	const StorageBinListItem & Get_reaction(void)const      { return(this->reaction); };
	const StorageBinListItem & Get_temperature(void)const   { return(this->temperature); };
	const StorageBinListItem & Get_pressure(void)const      { return(this->pressure); };
	const StorageBinListItem & Get_cell(void)const          { return(this->cell); };
protected:
	// update GetAllItems() if StorageBinListItem is added/removed
	StorageBinListItem solution;
	StorageBinListItem pp_assemblage;
	StorageBinListItem exchange;
	StorageBinListItem surface;
	StorageBinListItem ss_assemblage;
	StorageBinListItem gas_phase;
	StorageBinListItem kinetics;
	StorageBinListItem mix;
	StorageBinListItem reaction;
	StorageBinListItem temperature;
	StorageBinListItem pressure;
	const static std::vector < std::string > vopts;
	StorageBinListItem cell; // not included in GetAllItems
};


#endif // !defined(STORAGEBINLIST_H_INCLUDED)
