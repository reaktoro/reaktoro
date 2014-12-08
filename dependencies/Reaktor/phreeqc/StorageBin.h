#if !defined(STORAGEBIN_H_INCLUDED)
#define STORAGEBIN_H_INCLUDED
#include <cassert>				// assert
#include <map>					// std::map
#include <string>				// std::string
#include <list>					// std::list
#include <vector>				// std::vector

#include "System.h"
#include "PHRQ_io.h"
#include "PHRQ_base.h"
class cxxSolution;
class cxxExchange;
class cxxGasPhase;
class cxxKinetics;
class cxxPPassemblage;
class cxxSSassemblage;
class cxxSurface;
class cxxReaction;
class cxxTemperature;
class cxxUse;

class cxxStorageBin: public PHRQ_base
{

  public:
	cxxStorageBin(PHRQ_io *io=NULL);
	cxxStorageBin(cxxUse &use_ref, PHRQ_io *io=NULL);
	virtual ~cxxStorageBin();

	void Copy(int destination, int source);
	void Remove(int n);
	void Clear(void);

	cxxSolution *Get_Solution(int n_user);
	void Set_Solution(int n_user, cxxSolution * entity);
	void Set_Solution(int n_user, cxxSolution & entity);
	void Remove_Solution(int n_user);

	cxxExchange *Get_Exchange(int n_user);
	void Set_Exchange(int n_user, cxxExchange * entity);
	void Set_Exchange(int n_user, cxxExchange & entity);
	void Remove_Exchange(int n_user);

	cxxPPassemblage *Get_PPassemblage(int n_user);
	void Set_PPassemblage(int n_user, cxxPPassemblage * entity);
	void Set_PPassemblage(int n_user, cxxPPassemblage & entity);
	void Remove_PPassemblage(int n_user);

	cxxGasPhase *Get_GasPhase(int n_user);
	void Set_GasPhase(int n_user, cxxGasPhase * entity);
	void Set_GasPhase(int n_user, cxxGasPhase & entity);
	void Remove_GasPhase(int n_user);

	cxxSSassemblage *Get_SSassemblage(int n_user);
	void Set_SSassemblage(int n_user, cxxSSassemblage * entity);
	void Set_SSassemblage(int n_user, cxxSSassemblage & entity);
	void Remove_SSassemblage(int n_user);

	cxxKinetics *Get_Kinetics(int n_user);
	void Set_Kinetics(int n_user, cxxKinetics * entity);
	void Set_Kinetics(int n_user, cxxKinetics & entity);
	void Remove_Kinetics(int n_user);

	cxxSurface *Get_Surface(int n_user);
	void Set_Surface(int n_user, cxxSurface * entity);
	void Set_Surface(int n_user, cxxSurface & entity);
	void Remove_Surface(int n_user);

	cxxMix *Get_Mix(int n_user);
	void Set_Mix(int n_user, cxxMix * entity);
	void Set_Mix(int n_user, cxxMix & entity);
	void Remove_Mix(int n_user);

	cxxReaction *Get_Reaction(int n_user);
	void Set_Reaction(int n_user, cxxReaction * entity);
	void Set_Reaction(int n_user, cxxReaction & entity);
	void Remove_Reaction(int n_user);

	cxxTemperature *Get_Temperature(int n_user);
	void Set_Temperature(int n_user, cxxTemperature * entity);
	void Set_Temperature(int n_user, cxxTemperature & entity);
	void Remove_Temperature(int n_user);

	cxxPressure *Get_Pressure(int n_user);
	void Set_Pressure(int n_user, cxxPressure * entity);
	void Set_Pressure(int n_user, cxxPressure & entity);
	void Remove_Pressure(int n_user);

	cxxSystem &Get_System(void);
	void Set_System(cxxUse *use_ptr);
	void Set_System(int i);

	void dump_raw(std::ostream & s_oss, unsigned int indent) const;

	void dump_raw(std::ostream & s_oss, int i, unsigned int indent, int *n_out=NULL);

	void read_raw(CParser & parser);
	int read_raw_keyword(CParser & parser);

	void Add(cxxStorageBin &src, int n);

	//cxxSolution *mix_cxxSolutions(cxxMix &mixmap);
	cxxExchange *mix_cxxExchange(cxxMix & mixmap);

	std::map < int, cxxSolution > &Get_Solutions();
	std::map < int, cxxExchange > &Get_Exchangers();
	std::map < int, cxxGasPhase > &Get_GasPhases();
	std::map < int, cxxKinetics > &Get_Kinetics();
	std::map < int, cxxPPassemblage > &Get_PPassemblages();
	std::map < int, cxxSSassemblage > &Get_SSassemblages();
	std::map < int, cxxSurface > &Get_Surfaces();
	std::map < int, cxxMix > &Get_Mixes();
	std::map < int, cxxReaction > &Get_Reactions();
	std::map < int, cxxTemperature > &Get_Temperatures();
	std::map < int, cxxPressure > &Get_Pressures();

	cxxSystem & Get_system(void) {return system;};

  protected:
	// Tidied classes
	std::map < int, cxxSolution > Solutions;
	std::map < int, cxxExchange > Exchangers;
	std::map < int, cxxGasPhase > GasPhases;
	std::map < int, cxxKinetics > Kinetics;
	std::map < int, cxxPPassemblage > PPassemblages;
	std::map < int, cxxSSassemblage > SSassemblages;
	std::map < int, cxxSurface > Surfaces;

	// Reaction classes
	std::map < int, cxxMix > Mixes;
	std::map < int, cxxReaction > Reactions;
	std::map < int, cxxTemperature > Temperatures;
	std::map < int, cxxPressure > Pressures;
	cxxSystem system;

};

#endif // !defined(STORAGEBIN_H_INCLUDED)
