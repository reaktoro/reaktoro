// StorageBin.cxx: implementation of the cxxStorageBin class.
//
//////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
#pragma warning(disable : 4786)	// disable truncation warning (Only used by debugger)
#endif
#include <fstream>
#include <iostream>				// std::cout std::cerr
#include <cassert>				// assert
#include <algorithm>			// std::sort

#include "Utils.h"				// define first
#include "Phreeqc.h"
#include "NameDouble.h"
#include "StorageBin.h"
#include "SSassemblage.h"
#include "Solution.h"
#include "Exchange.h"
#include "GasPhase.h"
#include "cxxKinetics.h"
#include "PPassemblage.h"
#include "SS.h"
#include "SSassemblage.h"
#include "Surface.h"
#include "cxxMix.h"
#include "Reaction.h"
#include "Temperature.h"
#include "phqalloc.h"
#include "Use.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
cxxStorageBin::cxxStorageBin(PHRQ_io *io)
:
PHRQ_base(io)
{
	// default constructor for cxxStorageBin 
	this->system.Set_io(io);
	this->system.Initialize();
}

cxxStorageBin::cxxStorageBin(cxxUse &use_ref, PHRQ_io *io)
:
PHRQ_base(io)
{
	this->system.Set_io(io);
	this->system.Initialize();
	// Solution
	if (use_ref.Get_solution_ptr() != NULL)
	{
		this->Set_Solution(use_ref.Get_solution_ptr()->Get_n_user(), use_ref.Get_solution_ptr());
	}

	// Exchange
	if (use_ref.Get_exchange_ptr() != NULL)
	{
		this->Set_Exchange(use_ref.Get_exchange_ptr()->Get_n_user(), use_ref.Get_exchange_ptr());
	}

	// gas_phase
	if (use_ref.Get_gas_phase_ptr() != NULL)
	{
		this->Set_GasPhase(use_ref.Get_gas_phase_ptr()->Get_n_user(), use_ref.Get_gas_phase_ptr());
	}

	// kinetics
	if (use_ref.Get_kinetics_ptr() != NULL)
	{
		this->Set_Kinetics(use_ref.Get_kinetics_ptr()->Get_n_user(), use_ref.Get_kinetics_ptr());
	}

	// pp_assemblage
	if (use_ref.Get_pp_assemblage_ptr() != NULL)
	{
		this->Set_PPassemblage(use_ref.Get_pp_assemblage_ptr()->Get_n_user(), use_ref.Get_pp_assemblage_ptr());
	}

	// ss_assemblage
	if (use_ref.Get_ss_assemblage_ptr() != NULL)
	{
		this->Set_SSassemblage(use_ref.Get_ss_assemblage_ptr()->Get_n_user(), use_ref.Get_ss_assemblage_ptr());
	}

	// surface
	if (use_ref.Get_surface_ptr() != NULL)
	{
		this->Set_Surface(use_ref.Get_surface_ptr()->Get_n_user(), use_ref.Get_surface_ptr());
	}

	// mix
	if (use_ref.Get_mix_ptr() != NULL)
	{
		this->Set_Mix(use_ref.Get_mix_ptr()->Get_n_user(), use_ref.Get_mix_ptr());
	}

	// reaction
	if (use_ref.Get_reaction_ptr() != NULL)
	{
		this->Set_Reaction(use_ref.Get_reaction_ptr()->Get_n_user(), use_ref.Get_reaction_ptr());
	}

	// reaction temperature
	if (use_ref.Get_temperature_ptr() != NULL)
	{
		this->Set_Temperature(use_ref.Get_temperature_ptr()->Get_n_user(), use_ref.Get_temperature_ptr());
	}

	// reaction pressure
	if (use_ref.Get_pressure_ptr() != NULL)
	{
		this->Set_Pressure(use_ref.Get_pressure_ptr()->Get_n_user(), use_ref.Get_pressure_ptr());
	}
}
cxxStorageBin::~cxxStorageBin()
{
}
void
cxxStorageBin::Add(cxxStorageBin &src, int n)
{
	// Solution
	if (src.Get_Solution(n) != NULL)
	{
		this->Set_Solution(n, src.Get_Solution(n));
	}

	// Exchange
	if (src.Get_Exchange(n) != NULL)
	{
		this->Set_Exchange(n, src.Get_Exchange(n));
	}

	// gas_phase
	if (src.Get_GasPhase(n) != NULL)
	{
		this->Set_GasPhase(n, src.Get_GasPhase(n));
	}

	// kinetics
	if (src.Get_Kinetics(n) != NULL)
	{
		this->Set_Kinetics(n, src.Get_Kinetics(n));
	}

	// pp_assemblage
	if (src.Get_PPassemblage(n) != NULL)
	{
		this->Set_PPassemblage(n, src.Get_PPassemblage(n));
	}

	// ss_assemblage
	if (src.Get_SSassemblage(n) != NULL)
	{
		this->Set_SSassemblage(n, src.Get_SSassemblage(n));
	}

	// surface
	if (src.Get_Surface(n) != NULL)
	{
		this->Set_Surface(n, src.Get_Surface(n));
	}

	// mix
	if (src.Get_Mix(n) != NULL)
	{
		this->Set_Mix(n, src.Get_Mix(n));
	}

	// reaction
	if (src.Get_Reaction(n) != NULL)
	{
		this->Set_Reaction(n, src.Get_Reaction(n));
	}

	// reaction temperature
	if (src.Get_Temperature(n) != NULL)
	{
		this->Set_Temperature(n, src.Get_Temperature(n));
	}

	// reaction pressure
	if (src.Get_Pressure(n) != NULL)
	{
		this->Set_Pressure(n, src.Get_Pressure(n));
	}
}
void
cxxStorageBin::Copy(int destination, int source)
{
	if (destination == source)
		return;
	this->Remove(destination);
	// Solution
	{
		std::map < int, cxxSolution >::iterator it = this->Solutions.find(source);
		if (it != this->Solutions.end())
		{
			this->Set_Solution(destination, &(it->second));
		}
	}

	// Exchange
	{
		std::map < int, cxxExchange >::iterator it = this->Exchangers.find(source);
		if (it != this->Exchangers.end())
		{
			this->Set_Exchange(destination, &(it->second));
		}
	}

	// gas_phase
	{
		std::map < int, cxxGasPhase >::iterator it = this->GasPhases.find(source);
		if (it != this->GasPhases.end())
		{
			this->Set_GasPhase(destination, &(it->second));
		}
	}
	// kinetics
	{
		std::map < int, cxxKinetics >::iterator it = this->Kinetics.find(source);
		if (it != this->Kinetics.end())
		{
			this->Set_Kinetics(destination, &(it->second));
		}
	}
	// pp_assemblage
	{
		std::map < int, cxxPPassemblage >::iterator it = this->PPassemblages.find(source);
		if (it != this->PPassemblages.end())
		{
			this->Set_PPassemblage(destination, &(it->second));
		}
	}
	// ss_assemblage
	{
		std::map < int, cxxSSassemblage >::iterator it = this->SSassemblages.find(source);
		if (it != this->SSassemblages.end())
		{
			this->Set_SSassemblage(destination, &(it->second));
		}
	}
	// surface
	{
		std::map < int, cxxSurface >::iterator it = this->Surfaces.find(source);
		if (it != this->Surfaces.end())
		{
			this->Set_Surface(destination, &(it->second));
		}
	}
	// mix
	{
		std::map < int, cxxMix >::iterator it =	this->Mixes.find(source);
		if (it != this->Mixes.end())
		{
			this->Set_Mix(destination, &(it->second));
		}
	}
	// reaction
	{
		std::map < int, cxxReaction >::iterator it = this->Reactions.find(source);
		if (it != this->Reactions.end())
		{
			this->Set_Reaction(destination, &(it->second));
		}
	}
	// reaction temperature
	{
		std::map < int, cxxTemperature >::iterator it = this->Temperatures.find(source);
		if (it != this->Temperatures.end())
		{
			this->Set_Temperature(destination, &(it->second));
		}
	}

	// reaction pressure
	{
		this->Set_Pressure(destination, Utilities::Rxn_find(this->Pressures, source));
	}
}

cxxSolution *
cxxStorageBin::Get_Solution(int n_user)
{
	if (this->Solutions.find(n_user) != this->Solutions.end())
	{
		return (&(this->Solutions.find(n_user)->second));
	}
	return (NULL);
}
void 
cxxStorageBin::Set_Solution(int n_user, cxxSolution * entity)
{
	if (entity == NULL)
		return;
	Solutions[n_user] = *entity;
	Solutions.find(n_user)->second.Set_n_user_both(n_user);
}
void 
cxxStorageBin::Set_Solution(int n_user, cxxSolution & entity)
{
	Solutions[n_user] = entity;
	Solutions.find(n_user)->second.Set_n_user_both(n_user);
}
void 
cxxStorageBin::Remove_Solution(int n_user)
{
	Solutions.erase(n_user);
}

cxxExchange *
cxxStorageBin::Get_Exchange(int n_user)
{
	if (this->Exchangers.find(n_user) != this->Exchangers.end())
	{
		return (&(this->Exchangers.find(n_user)->second));
	}
	return (NULL);
}
void 
cxxStorageBin::Set_Exchange(int n_user, cxxExchange * entity)
{
	if (entity == NULL)
		return;
	Exchangers[n_user] = *entity;
	Exchangers.find(n_user)->second.Set_n_user_both(n_user);
}
void 
cxxStorageBin::Set_Exchange(int n_user, cxxExchange & entity)
{
	Exchangers[n_user] = entity;
	Exchangers.find(n_user)->second.Set_n_user_both(n_user);
}
void 
cxxStorageBin::Remove_Exchange(int n_user)
{
	Exchangers.erase(n_user);
}

cxxPPassemblage *
cxxStorageBin::Get_PPassemblage(int n_user)
{
	if (this->PPassemblages.find(n_user) != this->PPassemblages.end())
	{
		return (&(this->PPassemblages.find(n_user)->second));
	}
	return (NULL);
}
void 
cxxStorageBin::Set_PPassemblage(int n_user, cxxPPassemblage * entity)
{
	if (entity == NULL)
		return;
	PPassemblages[n_user] = *entity;
	PPassemblages.find(n_user)->second.Set_n_user_both(n_user);
}
void 
cxxStorageBin::Set_PPassemblage(int n_user, cxxPPassemblage & entity)
{
	PPassemblages[n_user] = entity;
	PPassemblages.find(n_user)->second.Set_n_user_both(n_user);
}
void 
cxxStorageBin::Remove_PPassemblage(int n_user)
{
	PPassemblages.erase(n_user);
}

cxxGasPhase *
cxxStorageBin::Get_GasPhase(int n_user)
{
	if (this->GasPhases.find(n_user) != this->GasPhases.end())
	{
		return (&(this->GasPhases.find(n_user)->second));
	}
	return (NULL);
}
void 
cxxStorageBin::Set_GasPhase(int n_user, cxxGasPhase * entity)
{
	if (entity == NULL)
		return;
	GasPhases[n_user] = *entity;
	GasPhases.find(n_user)->second.Set_n_user_both(n_user);
}
void 
cxxStorageBin::Set_GasPhase(int n_user, cxxGasPhase & entity)
{
	GasPhases[n_user] = entity;
	GasPhases.find(n_user)->second.Set_n_user_both(n_user);
}
void 
cxxStorageBin::Remove_GasPhase(int n_user)
{
	GasPhases.erase(n_user);
}

cxxSSassemblage *
cxxStorageBin::Get_SSassemblage(int n_user)
{
	if (this->SSassemblages.find(n_user) != this->SSassemblages.end())
	{
		return (&(this->SSassemblages.find(n_user)->second));
	}
	return (NULL);
}
void 
cxxStorageBin::Set_SSassemblage(int n_user, cxxSSassemblage * entity)
{
	if (entity == NULL)
		return;
	SSassemblages[n_user] = *entity;
	SSassemblages.find(n_user)->second.Set_n_user_both(n_user);
}
void 
cxxStorageBin::Set_SSassemblage(int n_user, cxxSSassemblage & entity)
{
	SSassemblages[n_user] = entity;
	SSassemblages.find(n_user)->second.Set_n_user_both(n_user);
}
void 
cxxStorageBin::Remove_SSassemblage(int n_user)
{
	SSassemblages.erase(n_user);
}

cxxKinetics *
cxxStorageBin::Get_Kinetics(int n_user)
{
	if (this->Kinetics.find(n_user) != this->Kinetics.end())
	{
		return (&(this->Kinetics.find(n_user)->second));
	}
	return (NULL);
}
void 
cxxStorageBin::Set_Kinetics(int n_user, cxxKinetics * entity)
{
	if (entity == NULL)
		return;
	Kinetics[n_user] = *entity;
	Kinetics.find(n_user)->second.Set_n_user_both(n_user);
}
void 
cxxStorageBin::Set_Kinetics(int n_user, cxxKinetics & entity)
{
	Kinetics[n_user] = entity;
	Kinetics.find(n_user)->second.Set_n_user_both(n_user);
}
void 
cxxStorageBin::Remove_Kinetics(int n_user)
{
	Kinetics.erase(n_user);
}

cxxSurface *
cxxStorageBin::Get_Surface(int n_user)
{
	if (this->Surfaces.find(n_user) != this->Surfaces.end())
	{
		return (&(this->Surfaces.find(n_user)->second));
	}
	return (NULL);
}
void 
cxxStorageBin::Set_Surface(int n_user, cxxSurface * entity)
{
	if (entity == NULL)
		return;
	Surfaces[n_user] = *entity;
	Surfaces.find(n_user)->second.Set_n_user_both(n_user);
}
void 
cxxStorageBin::Set_Surface(int n_user, cxxSurface & entity)
{
	Surfaces[n_user] = entity;
	Surfaces.find(n_user)->second.Set_n_user_both(n_user);
}
void 
cxxStorageBin::Remove_Surface(int n_user)
{
	Surfaces.erase(n_user);
}

cxxMix *
cxxStorageBin::Get_Mix(int n_user)
{
	if (this->Mixes.find(n_user) != this->Mixes.end())
	{
		return (&(this->Mixes.find(n_user)->second));
	}
	return (NULL);
}
void 
cxxStorageBin::Set_Mix(int n_user, cxxMix * entity)
{
	if (entity == NULL)
		return;
	Mixes[n_user] = *entity;
	Mixes.find(n_user)->second.Set_n_user_both(n_user);
}
void 
cxxStorageBin::Set_Mix(int n_user, cxxMix & entity)
{
	Mixes[n_user] = entity;
	Mixes.find(n_user)->second.Set_n_user_both(n_user);
}
void 
cxxStorageBin::Remove_Mix(int n_user)
{
	Mixes.erase(n_user);
}

cxxReaction *
cxxStorageBin::Get_Reaction(int n_user)
{
	if (this->Reactions.find(n_user) != this->Reactions.end())
	{
		return (&(this->Reactions.find(n_user)->second));
	}
	return (NULL);
}
void 
cxxStorageBin::Set_Reaction(int n_user, cxxReaction * entity)
{
	if (entity == NULL)
		return;
	Reactions[n_user] = *entity;
	Reactions.find(n_user)->second.Set_n_user_both(n_user);
}
void 
cxxStorageBin::Set_Reaction(int n_user, cxxReaction & entity)
{
	Reactions[n_user] = entity;
	Reactions.find(n_user)->second.Set_n_user_both(n_user);
}
void 
cxxStorageBin::Remove_Reaction(int n_user)
{
	Reactions.erase(n_user);
}

cxxTemperature *
cxxStorageBin::Get_Temperature(int n_user)
{
	if (this->Temperatures.find(n_user) != this->Temperatures.end())
	{
		return (&(this->Temperatures.find(n_user)->second));
	}
	return (NULL);
}

void 
cxxStorageBin::Set_Temperature(int n_user, cxxTemperature * entity)
{
	if (entity == NULL)
		return;
	Temperatures[n_user] = *entity;
	Temperatures.find(n_user)->second.Set_n_user_both(n_user);
}
void 
cxxStorageBin::Set_Temperature(int n_user, cxxTemperature & entity)
{
	Temperatures[n_user] = entity;
	Temperatures.find(n_user)->second.Set_n_user_both(n_user);
}
void 
cxxStorageBin::Remove_Temperature(int n_user)
{
	Temperatures.erase(n_user);
}

cxxPressure *
cxxStorageBin::Get_Pressure(int n_user)
{
	return Utilities::Rxn_find(this->Pressures, n_user);
}

void
cxxStorageBin::Set_Pressure(int n_user, cxxPressure * entity)
{
	if (entity == NULL)
		return;
	Pressures[n_user] = *entity;
	Pressures.find(n_user)->second.Set_n_user_both(n_user);
}
void
cxxStorageBin::Set_Pressure(int n_user, cxxPressure & entity)
{
	Pressures[n_user] = entity;
	Pressures.find(n_user)->second.Set_n_user_both(n_user);
}
void 

cxxStorageBin::Remove_Pressure(int n_user)
{
	Pressures.erase(n_user);
}

std::map < int, cxxSolution > &
cxxStorageBin::Get_Solutions()
{
	return this->Solutions;
}
std::map < int, cxxExchange > &
cxxStorageBin::Get_Exchangers()
{
	return this->Exchangers;
}
std::map < int, cxxGasPhase > &
cxxStorageBin::Get_GasPhases() 
{
	return this->GasPhases;
}
std::map < int, cxxKinetics > &
cxxStorageBin::Get_Kinetics()
{
	return this->Kinetics;
}
std::map < int, cxxPPassemblage > &
cxxStorageBin::Get_PPassemblages()
{
	return this->PPassemblages;
}
std::map < int, cxxSSassemblage > &
cxxStorageBin::Get_SSassemblages()
{
	return this->SSassemblages;
}
std::map < int, cxxSurface > &
cxxStorageBin::Get_Surfaces()
{
	return this->Surfaces;
}
std::map < int, cxxMix > &
cxxStorageBin::Get_Mixes()
{
	return this->Mixes;
}
std::map < int, cxxReaction > &
cxxStorageBin::Get_Reactions()
{
	return this->Reactions;
}
std::map < int, cxxTemperature > &
cxxStorageBin::Get_Temperatures()
{
	return this->Temperatures;
}
std::map < int, cxxPressure > &
cxxStorageBin::Get_Pressures()
{
	return this->Pressures;
}
#ifdef SKIP
void
cxxStorageBin::dump_xml(std::ostream & s_oss, unsigned int indent) const const
{
	unsigned int i;
	s_oss.precision(DBL_DIG - 1);
	std::string indent0(""), indent1(""), indent2("");
	for (i = 0; i < indent; ++i)
		indent0.append(Utilities::INDENT);
	for (i = 0; i < indent + 1; ++i)
		indent1.append(Utilities::INDENT);
	for (i = 0; i < indent + 2; ++i)
		indent2.append(Utilities::INDENT);

	// StorageBin element and attributes
	s_oss << indent0;
	s_oss << "<mix " << "\n";

	s_oss << indent1;
	s_oss << "pitzer_mix_gammas=\"" << this->
		pitzer_mix_gammas << "\"" << "\n";

	// components
	s_oss << indent1;
	s_oss << "<component " << "\n";
	for (std::list < cxxStorageBinComp >::const_iterator it =
		 mixComps.begin(); it != mixComps.end(); ++it)
	{
		it->dump_xml(s_oss, indent + 2);
	}

	return;
}
#endif

void
cxxStorageBin::dump_raw(std::ostream & s_oss, unsigned int indent) const
{
	// Dump all data
	s_oss.precision(DBL_DIG - 1);

	// Solutions
	Utilities::Rxn_dump_raw(Solutions, s_oss, indent);

	// Exchange
	Utilities::Rxn_dump_raw(Exchangers, s_oss, indent);

	// Gas Phases
	Utilities::Rxn_dump_raw(GasPhases, s_oss, indent);

	// Kinetics
	Utilities::Rxn_dump_raw(Kinetics, s_oss, indent);

	// PPassemblage
	Utilities::Rxn_dump_raw(PPassemblages, s_oss, indent);

	// SSassemblage
	Utilities::Rxn_dump_raw(SSassemblages, s_oss, indent);

	// Surface
	Utilities::Rxn_dump_raw(Surfaces, s_oss, indent);

	// Mix
	Utilities::Rxn_dump_raw(Mixes, s_oss, indent);

	// Reactions
	Utilities::Rxn_dump_raw(Reactions, s_oss, indent);

	// Temperature
	Utilities::Rxn_dump_raw(Temperatures, s_oss, indent);
}

void
cxxStorageBin::dump_raw(std::ostream & s_oss, int n, unsigned int indent, int *n_out)
{
	// Dump one user number, optionally change number from n to n_out
	int n_user_local = (n_out != NULL) ? *n_out : n;
	s_oss.precision(DBL_DIG - 1);

	// Solutions
	if (this->Get_Solution(n) != NULL)
	{
		this->Get_Solution(n)->dump_raw(s_oss, indent, &n_user_local);
	}

	// Exchange
	if (this->Get_Exchange(n) != NULL)
	{
		this->Get_Exchange(n)->dump_raw(s_oss, indent, &n_user_local);
	}

	// Gas Phases
	if (this->Get_GasPhase(n) != NULL)
	{
		this->Get_GasPhase(n)->dump_raw(s_oss, indent, &n_user_local);
	}

	// Kinetics
	if (this->Get_Kinetics(n) != NULL)
	{
		this->Get_Kinetics(n)->dump_raw(s_oss, indent, &n_user_local);
	}

	// PPassemblage
	if (this->Get_PPassemblage(n) != NULL)
	{
		this->Get_PPassemblage(n)->dump_raw(s_oss, indent, &n_user_local);
	}

	// SSassemblage
	if (this->Get_SSassemblage(n) != NULL)
	{
		this->Get_SSassemblage(n)->dump_raw(s_oss, indent, &n_user_local);
	}

	// Surface
	if (this->Get_Surface(n) != NULL)
	{
		this->Get_Surface(n)->dump_raw(s_oss, indent, &n_user_local);
	}

	// Mix
	if (this->Get_Mix(n) != NULL)
	{
		this->Get_Mix(n)->dump_raw(s_oss, indent, &n_user_local);
	}

	// Reaction
	if (this->Get_Reaction(n) != NULL)
	{
		this->Get_Reaction(n)->dump_raw(s_oss, indent, &n_user_local);
	}

	// Temperature
	if (this->Get_Temperature(n) != NULL)
	{
		this->Get_Temperature(n)->dump_raw(s_oss, indent, &n_user_local);
	}
}

void
cxxStorageBin::read_raw(CParser & parser)
{
	PHRQ_io::LINE_TYPE i;
	while ((i =
			parser.check_line("StorageBin read_raw", false, true, true,
							  true)) != PHRQ_io::LT_KEYWORD)
	{
		if (i == PHRQ_io::LT_EOF)
			return;				// PHRQ_io::LT_EOF;
	}

	for (;;)
	{
		switch (parser.next_keyword())
		{
		case Keywords::KEY_END:
		case Keywords::KEY_NONE:
			goto END_OF_SIMULATION_INPUT;
			break;
		case Keywords::KEY_SOLUTION_RAW:
			{
				cxxSolution entity(this->Get_io());
				entity.read_raw(parser);
				Solutions[entity.Get_n_user()] = entity;
			}
			break;
		case Keywords::KEY_EXCHANGE_RAW:
			{
				cxxExchange entity(this->Get_io());
				entity.read_raw(parser);
				Exchangers[entity.Get_n_user()] = entity;
			}
			break;
		case Keywords::KEY_GAS_PHASE_RAW:
			{
				cxxGasPhase entity(this->Get_io());
				entity.read_raw(parser);
				GasPhases[entity.Get_n_user()] = entity;
			}
			break;
		case Keywords::KEY_KINETICS_RAW:
			{
				cxxKinetics entity(this->Get_io());
				entity.read_raw(parser);
				Kinetics[entity.Get_n_user()] = entity;
			}
			break;

		case Keywords::KEY_EQUILIBRIUM_PHASES_RAW:
			{
				cxxPPassemblage entity(this->Get_io());
				entity.read_raw(parser);
				PPassemblages[entity.Get_n_user()] = entity;
			}
			break;

		case Keywords::KEY_SOLID_SOLUTIONS_RAW:
			{
				cxxSSassemblage entity;
				entity.read_raw(parser);
				SSassemblages[entity.Get_n_user()] = entity;
			}
			break;

		case Keywords::KEY_SURFACE_RAW:
			{
				cxxSurface entity(this->Get_io());
				entity.read_raw(parser);
				Surfaces[entity.Get_n_user()] = entity;
			}
			break;

		case Keywords::KEY_REACTION_TEMPERATURE_RAW:
			{
				cxxTemperature entity(this->Get_io());
				entity.read_raw(parser);
				Temperatures[entity.Get_n_user()] = entity;
			}
			break;

		case Keywords::KEY_REACTION_RAW:
			{
				cxxReaction entity;
				entity.read_raw(parser, true);
				Reactions[entity.Get_n_user()] = entity;
			}
			break;
		case Keywords::KEY_MIX_RAW:
			{
				cxxMix entity;
				entity.read_raw(parser);
				Mixes[entity.Get_n_user()] = entity;
			}
			break;
		default:
			{
				for (;;)
				{
					PHRQ_io::LINE_TYPE lt;
					lt = parser.check_line("read_raw", false, true, true, false);
					if (lt == PHRQ_io::LT_KEYWORD)
						break;
					if (lt == PHRQ_io::LT_EOF)
						goto END_OF_SIMULATION_INPUT;
				}
			}
			break;
		}
	}

  END_OF_SIMULATION_INPUT:
	return;						//PHRQ_io::LT_OK;
}

int
cxxStorageBin::read_raw_keyword(CParser & parser)
{
	PHRQ_io::LINE_TYPE i;
	int entity_number = -999;

	switch (parser.next_keyword())
	{
	case Keywords::KEY_NONE:
	case Keywords::KEY_END:
		while ((i =
				parser.check_line("StorageBin read_raw_keyword", false, true,
								  true, true)) != PHRQ_io::LT_KEYWORD)
		{
			if (i == PHRQ_io::LT_EOF)
				break;			// PHRQ_io::LT_EOF;
		}
		break;

	case Keywords::KEY_SOLUTION_RAW:
		{
			cxxSolution entity(this->Get_io());
			entity.read_raw(parser);
			Solutions[entity.Get_n_user()] = entity;
			entity_number = entity.Get_n_user();
		}
		break;

	case Keywords::KEY_EXCHANGE_RAW:
		{
			cxxExchange entity(this->Get_io());
			entity.read_raw(parser);
			Exchangers[entity.Get_n_user()] = entity;
			entity_number = entity.Get_n_user();
		}
		break;

	case Keywords::KEY_GAS_PHASE_RAW:
		{
			cxxGasPhase entity(this->Get_io());
			entity.read_raw(parser);
			GasPhases[entity.Get_n_user()] = entity;
			entity_number = entity.Get_n_user();
		}
		break;

	case Keywords::KEY_KINETICS_RAW:
		{
			cxxKinetics entity(this->Get_io());
			entity.read_raw(parser);
			Kinetics[entity.Get_n_user()] = entity;
			entity_number = entity.Get_n_user();
		}
		break;

	case Keywords::KEY_EQUILIBRIUM_PHASES_RAW:
		{
			cxxPPassemblage entity(this->Get_io());
			entity.read_raw(parser);
			PPassemblages[entity.Get_n_user()] = entity;
			entity_number = entity.Get_n_user();
		}
		break;

	case Keywords::KEY_SOLID_SOLUTIONS_RAW:
		{
			cxxSSassemblage entity;
			entity.read_raw(parser);
			SSassemblages[entity.Get_n_user()] = entity;
			entity_number = entity.Get_n_user();
		}
		break;

	case Keywords::KEY_SURFACE_RAW:
		{
			cxxSurface entity(this->Get_io());
			entity.read_raw(parser);
			Surfaces[entity.Get_n_user()] = entity;
			entity_number = entity.Get_n_user();
		}
		break;

	case Keywords::KEY_REACTION_TEMPERATURE_RAW:
		{
			cxxTemperature entity(this->Get_io());
			entity.read_raw(parser);
			Temperatures[entity.Get_n_user()] = entity;
			entity_number = entity.Get_n_user();
		}
		break;

	case Keywords::KEY_REACTION_RAW:
		{
			cxxReaction entity;
			entity.read_raw(parser, true);
			Reactions[entity.Get_n_user()] = entity;
			entity_number = entity.Get_n_user();
		}
		break;

	case Keywords::KEY_MIX_RAW:
		{
			cxxMix entity;
			entity.read_raw(parser);
			Mixes[entity.Get_n_user()] = entity;
			entity_number = entity.Get_n_user();
		}
		break;

	default:
		break;
	}
	return (entity_number);		//PHRQ_io::LT_OK;
}

void
cxxStorageBin::Remove(int n)
{
	// Solution
	this->Solutions.erase(n);

	// Exchanger
	this->Exchangers.erase(n);

	// GasPhase
	this->GasPhases.erase(n);

	// Kinetics
	this->Kinetics.erase(n);

	// PPassemblage
	this->PPassemblages.erase(n);

	// SSassemblage
	this->SSassemblages.erase(n);

	// Surface
	this->Surfaces.erase(n);

	// Mixes
	this->Mixes.erase(n);

	// Reactions
	this->Reactions.erase(n);

	// Temperature
	this->Temperatures.erase(n);

	// Pressure
	this->Pressures.erase(n);
}
void
cxxStorageBin::Clear(void) 
{
	// Delete all data

	// Solutions
	this->Solutions.clear();

	// Exchange
	this->Exchangers.clear();

	// Gas Phases
	this->GasPhases.clear();

	// Kinetics
	this->Kinetics.clear();

	// PPassemblage
	this->PPassemblages.clear();

	// SSassemblage
	this->SSassemblages.clear();

	// Surface
	this->Surfaces.clear();

	// Mix
	this->Mixes.clear();

	// Reactions
	this->Reactions.clear();

	// Temperature
	this->Temperatures.clear();

	// Pressure
	this->Pressures.clear();
}
#ifdef SKIP
cxxSolution *
cxxStorageBin::mix_cxxSolutions(cxxMix & mixmap)
{
/*
 *   mixes solutions based on cxxMix structure, returns new solution
 *   return solution must be freed by calling method
 */
	LDBLE intensive, extensive;
	cxxSolution *cxxsoln_ptr, *cxxsoln_ptr1;
/*
 *   Zero out global solution data
 */
	cxxsoln_ptr = new cxxSolution(0.0);
/*
 *   Determine sum of mixing fractions
 */
	extensive = 0.0;

	std::map < int, LDBLE >*mixcomps = mixmap.comps();

	std::map < int, LDBLE >::const_iterator it;
	for (it = mixcomps->begin(); it != mixcomps->end(); it++)
	{
		extensive += it->second;
	}
/*
 *   Add solutions 
 */
	for (it = mixcomps->begin(); it != mixcomps->end(); it++)
	{
		cxxsoln_ptr1 = &((this->Solutions.find(it->first))->second);
		if (cxxsoln_ptr1 == NULL)
		{
			error_string = sformatf(
					"Solution %d not found in mix_cxxSolutions.", it->first);
			error_msg(error_string, CONTINUE);
			phreeqc_ptr-> input_error++;
			return (NULL);
		}
		intensive = it->second / extensive;
		cxxsoln_ptr->add(*cxxsoln_ptr1, intensive, it->second);
	}
	return (cxxsoln_ptr);
}
#endif

#ifdef SKIP_OR_MOVE_TO_STRUCTURES
struct system *
cxxStorageBin::cxxStorageBin2system(Phreeqc * phreeqc_ptr, int n)
		//
		// make a system from storagebin
		//
{
	struct system *system_ptr =
		(struct system *) phreeqc_ptr-> PHRQ_malloc(sizeof(struct system));
	if (system_ptr == NULL)
		phreeqc_ptr-> malloc_error();

	// Solutions

	if (this->getSolution(n) != NULL)
	{
		//system_ptr->solution = (this->getSolution(n))->cxxSolution2solution(phreeqc_ptr);
		system_ptr->solution = phreeqc_ptr-> cxxSolution2solution(this->getSolution(n));
	}
	else
	{
		system_ptr->solution = NULL;
	}

	// Exchangers
	if (this->getExchange(n) != NULL)
	{
		//system_ptr->exchange = (this->getExchange(n))->cxxExchange2exchange(phreeqc_ptr);
		system_ptr->exchange = phreeqc_ptr-> cxxExchange2exchange(this->getExchange(n));
	}
	else
	{
		system_ptr->exchange = NULL;
	}

	// GasPhases
	if (this->getGasPhase(n) != NULL)
	{
		//system_ptr->gas_phase = (this->getGasPhase(n))->cxxGasPhase2gas_phase(phreeqc_ptr);
		system_ptr->gas_phase = phreeqc_ptr-> cxxGasPhase2gas_phase(this->getGasPhase(n));
	}
	else
	{
		system_ptr->gas_phase = NULL;
	}

	// Kinetics
	if (this->getKinetics(n) != NULL)
	{
		//system_ptr->kinetics = (this->getKinetics(n))->cxxKinetics2kinetics(phreeqc_ptr);
		system_ptr->kinetics = phreeqc_ptr-> cxxKinetics2kinetics(this->getKinetics(n));
		
	}
	else
	{
		system_ptr->kinetics = NULL;
	}

	// PPassemblages
	if (this->getPPassemblage(n) != NULL)
	{
		//system_ptr->pp_assemblage =
		//	(this->getPPassemblage(n))->cxxPPassemblage2pp_assemblage(phreeqc_ptr);
		system_ptr->pp_assemblage =
			phreeqc_ptr-> cxxPPassemblage2pp_assemblage(this->getPPassemblage(n));
	}
	else
	{
		system_ptr->pp_assemblage = NULL;
	}

	// SSassemblages
	if (this->getSSassemblage(n) != NULL)
	{
		//system_ptr->ss_assemblage =
		//	(this->getSSassemblage(n))->cxxSSassemblage2ss_assemblage(phreeqc_ptr);
		system_ptr->ss_assemblage =
			phreeqc_ptr-> cxxSSassemblage2ss_assemblage((this->getSSassemblage(n)));
	}
	else
	{
		system_ptr->ss_assemblage = NULL;
	}

	// Surfaces
	if (this->getSurface(n) != NULL)
	{
		//system_ptr->surface = (this->getSurface(n))->cxxSurface2surface(phreeqc_ptr);
		system_ptr->surface = phreeqc_ptr-> cxxSurface2surface((this->getSurface(n)));
	}
	else
	{
		system_ptr->surface = NULL;
	}
	return system_ptr;
}
#endif

#ifdef SKIP
cxxExchange *
cxxStorageBin::mix_cxxExchange(cxxMix & mixmap)
{
/*
 *   mixes exchangers based on cxxMix structure, returns new exchanger
 *   return exchanger must be freed by calling method
 */
	cxxExchange *new_exch_ptr, *old_exch_ptr;
/*
 *   Zero out global solution data
 */
	new_exch_ptr = new cxxExchange();

	std::map < int, LDBLE >::const_iterator it_mix;
	std::map < int, LDBLE >*mixcomps = mixmap.comps();

	// Pitzer_exchange_gammas
	it_mix = mixcomps->begin();
	old_exch_ptr = &((this->Exchangers.find(it_mix->first))->second);
	if (old_exch_ptr == NULL)
	{
		error_string = sformatf( "Exchange %d not found in mix_cxxExchange.",
				it_mix->first);
		error_msg(error_string, CONTINUE);
		phreeqc_ptr-> input_error++;
		return (NULL);
	}
	new_exch_ptr->set_pitzer_exchange_gammas(old_exch_ptr->
											 get_pitzer_exchange_gammas());
/*
 *   Make list of ExchComps
 */
	std::vector < cxxExchComp > ec_vector;
	std::vector < LDBLE >f_vector;
	//
	// make list of all exchange components and their mix fractions
	//
	for (it_mix = mixcomps->begin(); it_mix != mixcomps->end(); it_mix++)
	{
		old_exch_ptr = &((this->Exchangers.find(it_mix->first))->second);
		if (old_exch_ptr == NULL)
		{
			error_string = sformatf( "Exchange %d not found in mix_cxxExchange.",
					it_mix->first);
			error_msg(error_string, CONTINUE);
			phreeqc_ptr-> input_error++;
			return (NULL);
		}
		//  Add exchange components to vector ec_vector
		std::list < cxxExchComp >::const_iterator it_ec;
		std::list < cxxExchComp > &eclist = old_exch_ptr->get_exchComps();
		for (it_ec = eclist.begin(); it_ec != eclist.end(); it_ec++)
		{
			f_vector.push_back(it_mix->second);
			//cxxExchComp ec = *it_ec;
			//ec_vector.push_back(ec);
			ec_vector.push_back(*it_ec);
		}
	}
	//
	// Process list to make mixture
	//
	char *current_formula = ec_vector.begin()->get_formula();
	while (current_formula != NULL)
	{

		std::vector < cxxExchComp > ec_subvector;
		std::vector < LDBLE >f_subvector;
		std::vector < cxxExchComp >::iterator it_ec = ec_vector.begin();
		std::vector < LDBLE >::iterator it_f = f_vector.begin();
		current_formula = NULL;
		for (; it_ec != ec_vector.end(); it_ec++)
		{
			if (*it_f != 0)
			{
				if (current_formula == NULL)
					current_formula = it_ec->get_formula();
				if (it_ec->get_formula() == current_formula)
				{
					ec_subvector.push_back(*it_ec);
					f_subvector.push_back(*it_f);
					*it_f = 0;
					//ec_vector.erase(it_ec);
					//f_vector.erase(it_f);
				}
			}
			it_f++;
		}
		//
		//  mix ec_subvector to make
		// one exchange component
		//
		if (current_formula != NULL)
		{
			cxxExchComp new_comp(ec_subvector, f_subvector);
			new_exch_ptr->get_exchComps().push_back(new_comp);
		}
	}
	/*
	   std::ostringstream oss;
	   new_exch_ptr->dump_raw(oss, 0);
	   std::cerr << oss.str();
	 */
	return (new_exch_ptr);
}
#endif

cxxSystem &
cxxStorageBin::Get_System(void)
{
	return this->system;
}

void
cxxStorageBin::Set_System(cxxUse *use_ptr)
{
	// Initialize
	this->system.Initialize();
	// Solution
	if (use_ptr->Get_solution_ptr() != NULL)
	{
		std::map < int, cxxSolution >::iterator it =
			this->Solutions.find(use_ptr->Get_n_solution_user());
		if (it != this->Solutions.end())
		{
			this->system.Set_Solution(&(it->second));
		}
	}
	// Exchange
	if (use_ptr->Get_exchange_ptr() != NULL)
	{
		std::map < int, cxxExchange >::iterator it =
			this->Exchangers.find(use_ptr->Get_n_exchange_user());
		if (it != this->Exchangers.end())
		{
			this->system.Set_Exchange(&(it->second));
		}
	}
	// gas_phase
	if (use_ptr->Get_gas_phase_ptr() != NULL)
	{
		std::map < int, cxxGasPhase >::iterator it =
			this->GasPhases.find(use_ptr->Get_n_gas_phase_user());
		if (it != this->GasPhases.end())
		{
			this->system.Set_GasPhase(&(it->second));
		}
	}
	// kinetics
	if (use_ptr->Get_kinetics_ptr() != NULL)
	{
		std::map < int, cxxKinetics >::iterator it =
			this->Kinetics.find(use_ptr->Get_n_kinetics_user());
		if (it != this->Kinetics.end())
		{
			this->system.Set_Kinetics(&(it->second));
		}
	}
	// pp_assemblage
	if (use_ptr->Get_pp_assemblage_ptr() != NULL)
	{
		std::map < int, cxxPPassemblage >::iterator it =
			this->PPassemblages.find(use_ptr->Get_n_pp_assemblage_user());
		if (it != this->PPassemblages.end())
		{
			this->system.Set_PPassemblage(&(it->second));
		}
	}
	// ss_assemblage
	if (use_ptr->Get_ss_assemblage_ptr() != NULL)
	{
		std::map < int, cxxSSassemblage >::iterator it =
			this->SSassemblages.find(use_ptr->Get_n_ss_assemblage_user());
		if (it != this->SSassemblages.end())
		{
			this->system.Set_SSassemblage(&(it->second));
		}
	}
	// surface
	if (use_ptr->Get_surface_ptr() != NULL)
	{
		std::map < int, cxxSurface >::iterator it =
			this->Surfaces.find(use_ptr->Get_n_surface_user());
		if (it != this->Surfaces.end())
		{
			this->system.Set_Surface(&(it->second));
		}
	}
	// mix
	if (use_ptr->Get_mix_ptr() != NULL)
	{
		std::map < int, cxxMix >::iterator it =
			this->Mixes.find(use_ptr->Get_n_mix_user());
		if (it != this->Mixes.end())
		{
			this->system.Set_Mix(&(it->second));
		}
	}
	// reaction
	if (use_ptr->Get_reaction_ptr() != NULL)
	{
		std::map < int, cxxReaction >::iterator it =
			this->Reactions.find(use_ptr->Get_n_reaction_user());
		if (it != this->Reactions.end())
		{
			this->system.Set_Reaction(&(it->second));
		}
	}
	// reaction temperature
	if (use_ptr->Get_temperature_ptr() != NULL)
	{
		std::map < int, cxxTemperature >::iterator it =
			this->Temperatures.find(use_ptr->Get_n_temperature_user());
		if (it != this->Temperatures.end())
		{
			this->system.Set_Temperature(&(it->second));
		}
	}
	// reaction pressure
	if (use_ptr->Get_pressure_ptr() != NULL)
	{
		cxxPressure * p = Utilities::Rxn_find(this->Pressures, use_ptr->Get_n_pressure_user());
		if (p != NULL)
		{
			this->system.Set_Pressure(p);
		}
	}
}
void
cxxStorageBin::Set_System(int i)
{
	// Initialize
	this->system.Initialize();
	// Solution
	{
		std::map < int, cxxSolution >::iterator it = this->Solutions.find(i);
		if (it != this->Solutions.end())
		{
			this->system.Set_Solution(&(it->second));
		}
	}

	// Exchange
	{
		std::map < int, cxxExchange >::iterator it = this->Exchangers.find(i);
		if (it != this->Exchangers.end())
		{
			this->system.Set_Exchange(&(it->second));
		}
	}

	// gas_phase
	{
		std::map < int, cxxGasPhase >::iterator it = this->GasPhases.find(i);
		if (it != this->GasPhases.end())
		{
			this->system.Set_GasPhase(&(it->second));
		}
	}
	// kinetics
	{
		std::map < int, cxxKinetics >::iterator it = this->Kinetics.find(i);
		if (it != this->Kinetics.end())
		{
			this->system.Set_Kinetics(&(it->second));
		}
	}
	// pp_assemblage
	{
		std::map < int, cxxPPassemblage >::iterator it = this->PPassemblages.find(i);
		if (it != this->PPassemblages.end())
		{
			this->system.Set_PPassemblage(&(it->second));
		}
	}
	// ss_assemblage
	{
		std::map < int, cxxSSassemblage >::iterator it = this->SSassemblages.find(i);
		if (it != this->SSassemblages.end())
		{
			this->system.Set_SSassemblage(&(it->second));
		}
	}
	// surface
	{
		std::map < int, cxxSurface >::iterator it = this->Surfaces.find(i);
		if (it != this->Surfaces.end())
		{
			this->system.Set_Surface(&(it->second));
		}
	}
	// mix
	{
		std::map < int, cxxMix >::iterator it =	this->Mixes.find(i);
		if (it != this->Mixes.end())
		{
			this->system.Set_Mix(&(it->second));
		}
	}
	// reaction
	{
		std::map < int, cxxReaction >::iterator it = this->Reactions.find(i);
		if (it != this->Reactions.end())
		{
			this->system.Set_Reaction(&(it->second));
		}
	}
	// reaction temperature
	{
		std::map < int, cxxTemperature >::iterator it = this->Temperatures.find(i);
		if (it != this->Temperatures.end())
		{
			this->system.Set_Temperature(&(it->second));
		}
	}

	// reaction pressure
	{
		this->system.Set_Pressure(Utilities::Rxn_find(this->Pressures, i));
	}
}
