#if !defined(SURFACECHARGE_H_INCLUDED)
#define SURFACECHARGE_H_INCLUDED

#include <cassert>				// assert
#include <map>					// std::map
#include <string>				// std::string
#include <list>					// std::list
#include <vector>				// std::vector

#include "NameDouble.h"
#include "PHRQ_base.h"
class cxxSpeciesDL
{
public:
	cxxSpeciesDL()
	{
		g_moles = dg_g_moles = dx_moles = dh2o_moles = drelated_moles = 0;
	}
	LDBLE Get_g_moles(void) const {return g_moles;}
	void Set_g_moles(LDBLE t) {g_moles = t;}
	LDBLE Get_dg_g_moles(void) const {return dg_g_moles;}
	void Set_dg_g_moles(LDBLE t) {dg_g_moles = t;}
	LDBLE Get_dx_moles(void) const {return dx_moles;}
	void Set_dx_moles(LDBLE t) {dx_moles = t;}
	LDBLE Get_dh2o_moles(void) const {return dh2o_moles;}
	void Set_dh2o_moles(LDBLE t) {dh2o_moles = t;}
	LDBLE Get_drelated_moles(void) const {return drelated_moles;}
	void Set_drelated_moles(LDBLE t) {drelated_moles = t;}

	LDBLE *Get_g_moles_address(void) {return &g_moles;}
	LDBLE *Get_dg_g_moles_address(void) {return &dg_g_moles;}
	LDBLE *Get_dx_moles_address(void) {return &dx_moles;}
	LDBLE *Get_dh2o_moles_address(void) {return &dh2o_moles;}
	LDBLE *Get_drelated_moles_address(void) {return &drelated_moles;}


protected:
	LDBLE g_moles;
	LDBLE dg_g_moles;			/* g_moles*dgterm */
	LDBLE dx_moles;
	LDBLE dh2o_moles;			/* moles*g*Ws/Waq */
	LDBLE drelated_moles;		/* for related phase */
};
class cxxSurfDL
{
public:
	cxxSurfDL()
	{
		g = dg = psi_to_z = 0;
	}
	LDBLE Get_g(void) const {return g;}
	void Set_g(LDBLE t) {g = t;}
	LDBLE Get_dg(void) const {return dg;}
	void Set_dg(LDBLE t) {dg = t;}
	LDBLE Get_psi_to_z(void) const {return psi_to_z;}
	void Set_psi_to_z(LDBLE t) {psi_to_z = t;}
protected:
	LDBLE g;
	LDBLE dg;
	LDBLE psi_to_z;
};
class cxxSurfaceCharge: public PHRQ_base
{

public:

	cxxSurfaceCharge(PHRQ_io *io=NULL);
	cxxSurfaceCharge(struct surface_charge *, PHRQ_io *io=NULL);
	virtual ~cxxSurfaceCharge();

	struct master *Get_psi_master();
	void dump_xml(std::ostream & os, unsigned int indent = 0) const;
	void dump_raw(std::ostream & s_oss, unsigned int indent) const;
	void read_raw(CParser & parser, bool check = true);
	void add(const cxxSurfaceCharge & comp, LDBLE extensive);
	void multiply(LDBLE extensive);

	const std::string &Get_name() const	{return this->name;}
	void Set_name(const char * f) {this->name = f ? f : "";}
	LDBLE Get_specific_area() const {return this->specific_area;}
	void Set_specific_area(LDBLE d) {this->specific_area = d;}
	LDBLE Get_grams() const {return this->grams;}
	void Set_grams(LDBLE d) {this->grams = d;}
	LDBLE Get_charge_balance() const {return this->charge_balance;}
	void Set_charge_balance(LDBLE d) {this->charge_balance = d;}
	LDBLE Get_mass_water() const {return this->mass_water;}
	void Set_mass_water(LDBLE d) {this->mass_water = d;}
	LDBLE Get_la_psi() const {return this->la_psi;}
	void Set_la_psi(LDBLE d) {this->la_psi = d;}
	LDBLE Get_capacitance0() const {return this->capacitance[0];}
	void Set_capacitance0(LDBLE d) {this->capacitance[0] = d;}
	LDBLE Get_capacitance1() const {return this->capacitance[1];}
	void Set_capacitance1(LDBLE d) {this->capacitance[1] = d;}
	LDBLE Get_sigma0() const {return this->sigma0;}
	void Set_sigma0(LDBLE d) {this->sigma0 = d;}
	LDBLE Get_sigma1() const {return this->sigma1;}
	void Set_sigma1(LDBLE d) {this->sigma1 = d;}
	LDBLE Get_sigma2() const {return this->sigma2;}
	void Set_sigma2(LDBLE d) {this->sigma2 = d;}
	LDBLE Get_sigmaddl() const {return this->sigmaddl;}
	void Set_sigmaddl(LDBLE d) {this->sigmaddl = d;}
	const cxxNameDouble & Get_diffuse_layer_totals(void) const {return this->diffuse_layer_totals;}
	void Set_diffuse_layer_totals(cxxNameDouble & nd) {this->diffuse_layer_totals = nd;}
	std::map<LDBLE, cxxSurfDL> &Get_g_map(void) {return g_map;}
	void Set_g_map(std::map<LDBLE, cxxSurfDL> & t) {g_map = t;}

protected:
	std::string name;
	LDBLE specific_area;
	LDBLE grams;
	LDBLE charge_balance;
	LDBLE mass_water;
	LDBLE la_psi;
	LDBLE capacitance[2];
	cxxNameDouble diffuse_layer_totals;
	// workspace variables
	LDBLE sigma0, sigma1, sigma2, sigmaddl;
	std::map<LDBLE, cxxSurfDL> g_map;
	const static std::vector < std::string > vopts;
};

#endif // !defined(SURFACECHARGE_H_INCLUDED)
