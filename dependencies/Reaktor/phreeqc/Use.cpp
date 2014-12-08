#include <stdio.h>
#include "Use.h"

cxxUse::cxxUse()
{
	this->init();

}


cxxUse::~cxxUse(void)
{
}
void cxxUse::
init(void)
{
	solution_in = false;
	n_solution_user = -999;
	solution_ptr = NULL;

	pp_assemblage_in = false;
	n_pp_assemblage_user = -999;
	pp_assemblage_ptr = NULL;

	mix_in = false;
	n_mix_user = -999;
	mix_ptr = NULL;
	n_mix_user_orig = -999;

	reaction_in = false;
	n_reaction_user = -999;
	reaction_ptr = NULL;

	exchange_in = false;
	n_exchange_user = -999;
	exchange_ptr = NULL;

	kinetics_in = false;
	n_kinetics_user = -999;
	kinetics_ptr = NULL;

	surface_in = false;
	n_surface_user = -999;
	surface_ptr = NULL;

	pressure_in = false;
	n_pressure_user = -999;
	pressure_ptr = NULL;

	temperature_in = false;
	n_temperature_user = -999;
	temperature_ptr = NULL;

	inverse_in = false;
	n_inverse_user = -999;
	inverse_ptr = NULL;

	gas_phase_in = false;
	n_gas_phase_user = -999;
	gas_phase_ptr = NULL;

	ss_assemblage_in = false;
	n_ss_assemblage_user = -999;
	ss_assemblage_ptr = NULL;

	trans_in = false;
	advect_in = false;
}


