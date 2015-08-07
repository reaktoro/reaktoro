#include <algorithm>			// std::replace 

#include "Phreeqc.h"
#include "System.h"
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

cxxSystem::cxxSystem(PHRQ_io *io)
:
PHRQ_base(io)
{
	this->solution = NULL;
	this->exchange = NULL;
	this->ppassemblage = NULL;
	this->gasphase = NULL;
	this->ssassemblage = NULL;
	this->kinetics = NULL;
	this->surface = NULL;
	this->mix = NULL;
	this->reaction = NULL;
	this->temperature = NULL;
	this->pressure = NULL;
} 
cxxSystem::~cxxSystem(void) 
{
} 
void
cxxSystem::Initialize(void) 
{
	this->solution = NULL;
	this->exchange = NULL;
	this->ppassemblage = NULL;
	this->gasphase = NULL;
	this->ssassemblage = NULL;
	this->kinetics = NULL;
	this->surface = NULL;
	this->mix = NULL;
	this->reaction = NULL;
	this->temperature = NULL;
	this->pressure = NULL;
} 
void
cxxSystem::totalize(Phreeqc * phreeqc_ptr) 
{
	//initialize
	this->totals.clear();
	//add solution
	if (this->solution != NULL)
	{
		char token[MAX_LENGTH];
		strcpy(token, "O");
		this->totals[token] = this->solution->Get_total_o();
		strcpy(token, "H");
		this->totals[token] = this->solution->Get_total_h();
		strcpy(token, "Charge");
		this->totals[token] = this->solution->Get_cb();
		this->totals.add_extensive(this->solution->Get_totals(), 1.0);
	}
	if (this->exchange != NULL)
	{
		this->exchange->totalize();
		this->totals.add_extensive(this->exchange->Get_totals(), 1.0);
	}
	if (this->ppassemblage != NULL)
	{
		this->ppassemblage->totalize(phreeqc_ptr);
		this->totals.add_extensive(this->ppassemblage->Get_assemblage_totals(), 1.0);
	}
	if (this->gasphase != NULL)
	{
		this->gasphase->totalize(phreeqc_ptr);
		this->totals.add_extensive(this->gasphase->Get_totals(), 1.0);
	}
	if (this->ssassemblage != NULL)
	{
		this->ssassemblage->totalize(phreeqc_ptr);
		this->totals.add_extensive(this->ssassemblage->Get_totals(), 1.0);
	}
	if (this->surface != NULL)
	{
		this->surface->totalize();
		this->totals.add_extensive(this->surface->Get_totals(), 1.0);
	}
	return;
}
