#if !defined(SOLUTION_H_INCLUDED)
#define SOLUTION_H_INCLUDED

#include <cassert>				// assert
#include <map>					// std::map
#include <string>				// std::string
#include <vector>				// std::vector
#include <iostream>
#include "NumKeyword.h"
#include "SolutionIsotope.h"
#include "NameDouble.h"
#include "PHRQ_base.h"
#include "PHRQ_io.h"
#include "ISolution.h"
class cxxMix;

class cxxSolution:public cxxNumKeyword
{

  public:
	cxxSolution(PHRQ_io *io=NULL);
	cxxSolution(const cxxSolution &old_sol);
	const cxxSolution & operator = (const cxxSolution &rhs);
	cxxSolution(std::map < int, cxxSolution > &solution_map,
				  cxxMix & mx, int n_user, PHRQ_io *io=NULL);
	virtual ~cxxSolution();

	bool Get_new_def() const  {return this->new_def;}
	void Set_new_def(bool p)  {this->new_def = p;}
	LDBLE Get_patm() const    {return this->patm;}
	void Set_patm(LDBLE p)    {this->patm = p;}
	LDBLE Get_tc() const      {return this->tc;}
	void Set_tc(LDBLE l_tc)   {this->tc = l_tc;}
	LDBLE Get_ph() const      {return this->ph;}
	void Set_ph(LDBLE pH)     {this->ph = pH;}
	LDBLE Get_pe() const      {return this->pe;}
	void Set_pe(LDBLE l_pe)   {this->pe = l_pe;}
	LDBLE Get_mu() const      {return this->mu;}
	void Set_mu(LDBLE l_mu)   {this->mu = l_mu;}
	LDBLE Get_ah2o() const                              {return this->ah2o;}
	void Set_ah2o(LDBLE l_ah2o)                         {this->ah2o = l_ah2o;}
	LDBLE Get_total_h() const                           {return this->total_h;}
	void Set_total_h(LDBLE l_total_h)                   {this->total_h = l_total_h;}
	LDBLE Get_total_o() const                           {return this->total_o;}
	void Set_total_o(LDBLE l_total_o)                   {this->total_o = l_total_o;}
	LDBLE Get_cb() const                                {return this->cb;}
	void Set_cb(LDBLE l_cb)                             {this->cb = l_cb;}
	LDBLE Get_density() const                           {return this->density;}
	void Set_density(LDBLE l_density)                   {this->density = l_density;}
	LDBLE Get_mass_water() const                        {return this->mass_water;}
	void Set_mass_water(LDBLE l_mass_water)             {this->mass_water = l_mass_water;}
	LDBLE Get_total_alkalinity() const                  {return this->total_alkalinity;}
	void Set_total_alkalinity(LDBLE l_total_alkalinity) {this->total_alkalinity = l_total_alkalinity;}
	LDBLE Get_soln_vol() const                          {return this->soln_vol;}
	void Set_soln_vol(LDBLE t)                          {this->soln_vol = t;}
	cxxNameDouble & Get_totals(void)                    {return this->totals;}
	const cxxNameDouble & Get_totals(void)const         {return this->totals;}
	void Set_totals(cxxNameDouble & nd)
	{
		this->totals = nd;
		this->totals.type = cxxNameDouble::ND_ELT_MOLES;
	}
	cxxNameDouble & Get_master_activity(void)           {return this->master_activity;}
	cxxNameDouble & Get_species_gamma(void)             {return this->species_gamma;}
	std::map<int, double> & Get_species_map(void)       {return this->species_map;}
	std::map < std::string, cxxSolutionIsotope > & Get_isotopes(void)             {return this->isotopes;}	
	const std::map < std::string, cxxSolutionIsotope > & Get_isotopes(void)const  {return this->isotopes;}	
	void Set_isotopes(const std::map < std::string, cxxSolutionIsotope > &iso ) {this->isotopes = iso;}	
	cxxISolution *Get_initial_data()                    {return this->initial_data;}
	const cxxISolution *Get_initial_data()const         {return this->initial_data;}
	void Set_initial_data(const cxxISolution * id)
	{
		if (this->initial_data != NULL)
			delete this->initial_data;
		this->initial_data = new cxxISolution(*id);
	}
	void Create_initial_data()
	{
		if (this->initial_data != NULL)
			delete initial_data;
		initial_data = new cxxISolution;
	}
	void Destroy_initial_data()
	{
		if (this->initial_data != NULL)
			delete initial_data;
		initial_data = NULL;
	}

	void clear_totals()					{this->totals.clear();}
	void clear_master_activity()    	{this->master_activity.clear();}

	void zero();
	void add(const cxxSolution & addee, LDBLE extensive);
	void Add_isotopes(const std::map < std::string, cxxSolutionIsotope > &old, LDBLE intensive, LDBLE extensive);
	void multiply(LDBLE extensive);
	void Multiply_isotopes(LDBLE extensive);
	// not checked
	void dump_xml(std::ostream & os, unsigned int indent = 0) const;
	void dump_raw(std::ostream & s_oss, unsigned int indent, int *n_out=NULL) const;

	LDBLE Get_master_activity(char *string) const;
	void Set_master_activity(char *string, LDBLE value);	
	LDBLE Get_total(const char *string) const;
	LDBLE Get_total_element(const char *string) const;
	void Set_total(char *string, LDBLE value);
	void read_raw(CParser & parser, bool check = true);
	//void modify_activities(const cxxSolution & original);
	//void Simplify_totals();
	void Update(const cxxNameDouble &nd);
	void Update(LDBLE h_tot, LDBLE o_tot, LDBLE charge, const cxxNameDouble &nd);
	void Update_activities(const cxxNameDouble &original_tot);

  protected:
	bool new_def;
	LDBLE patm;
	LDBLE tc;
	LDBLE ph;
	LDBLE pe;
	LDBLE mu;
	LDBLE ah2o;
	LDBLE total_h;
	LDBLE total_o;
	LDBLE cb;
	LDBLE mass_water;
	LDBLE density;
	LDBLE soln_vol;
	LDBLE total_alkalinity;
	cxxNameDouble totals;
	cxxNameDouble master_activity;
	cxxNameDouble species_gamma;
	//cxxSolutionIsotopeList isotopes;
	std::map < std::string, cxxSolutionIsotope > isotopes;
	cxxISolution *initial_data;
	const static std::vector < std::string > vopts;
	std::map<int, double> species_map;
};

#endif // !defined(SOLUTION_H_INCLUDED)
