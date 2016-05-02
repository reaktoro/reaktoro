#include "Serializer.h"
#include "Phreeqc.h"
#include "Utils.h"
#include "Solution.h"
#include "Exchange.h"
#include "Temperature.h"
#include "GasPhase.h"
#include "cxxKinetics.h"
#include "PPassemblage.h"
#include "SSassemblage.h"
#include "Surface.h"
Serializer::Serializer(PHRQ_io *io)
	: PHRQ_base(io)
{
}
Serializer::~Serializer(void)
{
}
bool Serializer::Serialize(Phreeqc &phreeqc_ref, int start, int end, bool include_t, bool include_p, PHRQ_io *io)
{
	for (int i = start; i <= end; i++)
	{
		cxxSolution *soln_ptr = Utilities::Rxn_find(phreeqc_ref.Get_Rxn_solution_map(), i);
		if (soln_ptr)
		{
			ints.push_back((int) PT_SOLUTION);
			soln_ptr->Serialize(this->dictionary, this->ints, this->doubles);
		}
		// Exchangers
		{
			cxxExchange *entity_ptr = Utilities::Rxn_find(phreeqc_ref.Get_Rxn_exchange_map(), i);
			if (entity_ptr)
			{
				ints.push_back((int) PT_EXCHANGE);
				entity_ptr->Serialize(this->dictionary, this->ints, this->doubles);
			}
		}
		// GasPhases
		{
			cxxGasPhase *entity_ptr = Utilities::Rxn_find(phreeqc_ref.Get_Rxn_gas_phase_map(), i);
			if (entity_ptr)
			{
				ints.push_back((int) PT_GASPHASE);
				entity_ptr->Serialize(this->dictionary, this->ints, this->doubles);
			}
		}
		// Kinetics
		{
			cxxKinetics *entity_ptr = Utilities::Rxn_find(phreeqc_ref.Get_Rxn_kinetics_map(), i);
			if (entity_ptr)
			{
				ints.push_back((int) PT_KINETICS);
				entity_ptr->Serialize(this->dictionary, this->ints, this->doubles);
			}
		}
		// PPassemblages
		{
			cxxPPassemblage *entity_ptr = Utilities::Rxn_find(phreeqc_ref.Get_Rxn_pp_assemblage_map(), i);
			if (entity_ptr)
			{
				ints.push_back((int) PT_PPASSEMBLAGE);
				entity_ptr->Serialize(this->dictionary, this->ints, this->doubles);
			}
		}
		// SSassemblages
		{
			cxxSSassemblage *entity_ptr = Utilities::Rxn_find(phreeqc_ref.Get_Rxn_ss_assemblage_map(), i);
			if (entity_ptr)
			{
				ints.push_back((int) PT_SSASSEMBLAGE);
				entity_ptr->Serialize(this->dictionary, this->ints, this->doubles);
			}
		}
		// Surfaces
		{
			cxxSurface *entity_ptr = Utilities::Rxn_find(phreeqc_ref.Get_Rxn_surface_map(), i);
			if (entity_ptr)
			{
				ints.push_back((int) PT_SURFACES);
				entity_ptr->Serialize(this->dictionary, this->ints, this->doubles);
			}
		}
		// Temperature
		if (include_t)
		{		
			cxxTemperature *entity_ptr = Utilities::Rxn_find(phreeqc_ref.Get_Rxn_temperature_map(), i);
			if (entity_ptr)
			{
				ints.push_back((int) PT_TEMPERATURE);
				entity_ptr->Serialize(this->dictionary, this->ints, this->doubles);
			}
		}
		// Pressure
		if (include_p)
		{		
			cxxPressure *entity_ptr = Utilities::Rxn_find(phreeqc_ref.Get_Rxn_pressure_map(), i);
			if (entity_ptr)
			{
				ints.push_back((int) PT_PRESSURE);
				entity_ptr->Serialize(this->dictionary, this->ints, this->doubles);
			}
		}			
	}	
	return true;
}

bool 
Serializer::Deserialize(Phreeqc &phreeqc_ref, Dictionary &dictionary, std::vector<int> &ints, std::vector<double> &doubles)
{
	int ii = 0;
	int dd = 0;
	while (ii < (int) ints.size())
	{
		PACK_TYPE type = (PACK_TYPE) ints[ii++];
		switch (type)
		{
		case PT_SOLUTION:	
			{
				cxxSolution soln;
				soln.Deserialize(dictionary, ints, doubles, ii, dd);
				int n_user = soln.Get_n_user();
				//std::cerr << "unpacked solution " << n_user << std::endl;
				phreeqc_ref.Get_Rxn_solution_map()[n_user] = soln;
			}
			break;
		case PT_EXCHANGE:	
			{
				cxxExchange entity;
				entity.Deserialize(dictionary, ints, doubles, ii, dd);
				int n_user = entity.Get_n_user();
				phreeqc_ref.Get_Rxn_exchange_map()[n_user] = entity;
			}
			break;
		case PT_GASPHASE:
			{
				cxxGasPhase entity;
				entity.Deserialize(dictionary, ints, doubles, ii, dd);
				int n_user = entity.Get_n_user();
				phreeqc_ref.Get_Rxn_gas_phase_map()[n_user] = entity;
			}
			break;
		case PT_KINETICS:	
			{
				cxxKinetics entity;
				entity.Deserialize(dictionary, ints, doubles, ii, dd);
				int n_user = entity.Get_n_user();
				phreeqc_ref.Get_Rxn_kinetics_map()[n_user] = entity;
			}
			break;
		case PT_PPASSEMBLAGE:	
			{
				cxxPPassemblage entity;
				entity.Deserialize(dictionary, ints, doubles, ii, dd);
				int n_user = entity.Get_n_user();
				//std::cerr << "unpacked pp assemblage " << n_user << std::endl;
				phreeqc_ref.Get_Rxn_pp_assemblage_map()[n_user] = entity;
			}
			break;
		case PT_SSASSEMBLAGE:	
			{
				cxxSSassemblage entity;
				entity.Deserialize(dictionary, ints, doubles, ii, dd);
				int n_user = entity.Get_n_user();
				phreeqc_ref.Get_Rxn_ss_assemblage_map()[n_user] = entity;
			}
			break;
		case PT_SURFACES:	
			{
				cxxSurface entity;
				entity.Deserialize(dictionary, ints, doubles, ii, dd);
				int n_user = entity.Get_n_user();
				phreeqc_ref.Get_Rxn_surface_map()[n_user] = entity;
			}
			break;
		case PT_TEMPERATURE:
			{
				cxxTemperature entity;
				entity.Deserialize(dictionary, ints, doubles, ii, dd);
				int n_user = entity.Get_n_user();
				phreeqc_ref.Get_Rxn_temperature_map()[n_user] = entity;
			}
			break;
		case PT_PRESSURE:	
			{
				cxxPressure entity;
				entity.Deserialize(dictionary, ints, doubles, ii, dd);
				int n_user = entity.Get_n_user();
				phreeqc_ref.Get_Rxn_pressure_map()[n_user] = entity;
			}
			break;
		default:
#if !defined(R_SO)
			std::cerr << "Unknown pack type in deserialize " << type << std::endl;
			exit(4);
#endif
			break;
		}
	}
	return true;
}
