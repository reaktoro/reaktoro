#if !defined(SURFACE_H_INCLUDED)
#define SURFACE_H_INCLUDED

#include <cassert>				// assert
#include <map>					// std::map
#include <string>				// std::string
#include <list>					// std::list
#include <vector>				// std::vector

#include "NumKeyword.h"
#include "SurfaceComp.h"
#include "SurfaceCharge.h"
class cxxMix;

class cxxSurface:public cxxNumKeyword
{

public:
	enum SURFACE_TYPE
	{ UNKNOWN_DL, NO_EDL, DDL, CD_MUSIC, CCM };
	enum DIFFUSE_LAYER_TYPE
	{ NO_DL, BORKOVEK_DL, DONNAN_DL };
	enum SITES_UNITS
	{ SITES_ABSOLUTE, SITES_DENSITY };

	cxxSurface(PHRQ_io *io=NULL);
	cxxSurface(std::map < int, cxxSurface > &entity_map, cxxMix & mx,
				 int n_user, PHRQ_io *io=NULL);
	~cxxSurface();

	//void dump_xml(std::ostream & os, unsigned int indent = 0) const;
	void dump_raw(std::ostream & s_oss, unsigned int indent, int *n_out=NULL) const;
	void read_raw(CParser & parser, bool check = true);
	bool Get_related_phases(void) const;
	bool Get_related_rate(void) const;
	void totalize();
	cxxSurfaceComp *Find_comp(std::string str);
	cxxSurfaceCharge *Find_charge(const std::string str);
	const cxxSurfaceCharge *Find_charge(const std::string str)const;
	void Sort_comps();

	void add(const cxxSurface & addee, LDBLE extensive);
	void multiply(LDBLE extensive);

	std::vector < cxxSurfaceComp > & Get_surface_comps() {return this->surface_comps;}
	const std::vector < cxxSurfaceComp > & Get_surface_comps()const {return this->surface_comps;}
	void Set_surface_comps(std::vector < cxxSurfaceComp > &sc) {this->surface_comps = sc;}
	std::vector < cxxSurfaceCharge > & Get_surface_charges() {return this->surface_charges;}
	void Set_surface_charges(std::vector < cxxSurfaceCharge > &sc) {this->surface_charges = sc;}
	bool Get_new_def(void) {return new_def;}
	void Set_new_def(bool tf) {new_def = tf;}
	SURFACE_TYPE Get_type(void) const {return this->type;}
	void Set_type(SURFACE_TYPE t) {this->type = t;}
	DIFFUSE_LAYER_TYPE Get_dl_type(void) const {return dl_type;}
	void Set_dl_type(DIFFUSE_LAYER_TYPE t) {dl_type = t;}
	SITES_UNITS Get_sites_units(void) const {return sites_units;}
	void Set_sites_units(SITES_UNITS u) {sites_units = u;}
	bool Get_only_counter_ions(void) const {return only_counter_ions;}
	void Set_only_counter_ions(bool tf) {only_counter_ions = tf;}
	LDBLE Get_thickness(void) const {return thickness;}
	void Set_thickness(LDBLE t) {thickness = t;}
	LDBLE Get_debye_lengths(void) const {return debye_lengths;}
	void Set_debye_lengths(LDBLE t) {debye_lengths = t;}
	LDBLE Get_DDL_viscosity(void) const {return DDL_viscosity;}
	void Set_DDL_viscosity(LDBLE t) {DDL_viscosity = t;}
	LDBLE Get_DDL_limit(void) const {return DDL_limit;}
	void Set_DDL_limit(LDBLE t) {DDL_limit = t;}
	bool Get_transport(void) const {return transport;}
	void Set_transport(bool tf) {transport = tf;}
	cxxNameDouble & Get_totals() {return this->totals;}
	void Get_totals(cxxNameDouble &nd) {this->totals = nd;}
	bool Get_solution_equilibria(void)const {return solution_equilibria;}
	void Set_solution_equilibria(bool tf) {solution_equilibria = tf;}
	int Get_n_solution(void)const {return n_solution;}
	void Set_n_solution(int i) {n_solution = i;}
	void Serialize(Dictionary & dictionary, std::vector < int >&ints, std::vector < double >&doubles);
	void Deserialize(Dictionary & dictionary, std::vector < int >&ints, std::vector < double >&doubles, int &ii, int &dd);
	
protected:
	std::vector < cxxSurfaceComp > surface_comps;
	std::vector < cxxSurfaceCharge > surface_charges;
	bool new_def;
	SURFACE_TYPE type;
	DIFFUSE_LAYER_TYPE dl_type;
	SITES_UNITS sites_units;
	bool only_counter_ions;
	LDBLE thickness;
	LDBLE debye_lengths;
	LDBLE DDL_viscosity;
	LDBLE DDL_limit;
	bool transport;
	cxxNameDouble totals;
	bool solution_equilibria;
	int n_solution;
	const static std::vector < std::string > vopts;
};

#endif // !defined(SURFACE_H_INCLUDED)
