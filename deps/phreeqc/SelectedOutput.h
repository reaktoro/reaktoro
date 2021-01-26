#if !defined(SELECTEDOUTPUT_H_INCLUDED)
#define SELECTEDOUTPUT_H_INCLUDED
#include <string>				// std::string
#include <vector>
#include <map>
#include "NumKeyword.h"

class SelectedOutput:public cxxNumKeyword
{
public:
	SelectedOutput(int n=1, PHRQ_io *io=NULL);
	~SelectedOutput(void);

	void Reset(bool tf);

	// vector getters
	inline std::vector< std::pair< std::string, void * > > & Get_totals(void)           {return this->totals;}
	inline std::vector< std::pair< std::string, void * > > & Get_molalities(void)       {return this->molalities;}
	inline std::vector< std::pair< std::string, void * > > & Get_activities(void)       {return this->activities;}
	inline std::vector< std::pair< std::string, void * > > & Get_pure_phases(void)      {return this->pure_phases;}
	inline std::vector< std::pair< std::string, void * > > & Get_si(void)               {return this->si;}
	inline std::vector< std::pair< std::string, void * > > & Get_gases(void)            {return this->gases;}
	inline std::vector< std::pair< std::string, void * > > & Get_s_s(void)              {return this->s_s;}
	inline std::vector< std::pair< std::string, void * > > & Get_kinetics(void)         {return this->kinetics;}
	inline std::vector< std::pair< std::string, void * > > & Get_isotopes(void)         {return this->isotopes;}
	inline std::vector< std::pair< std::string, void * > > & Get_calculate_values(void) {return this->calculate_values;}

	// const vector getters
	inline const std::vector< std::pair< std::string, void * > > & Get_totals(void)const           {return this->totals;}
	inline const std::vector< std::pair< std::string, void * > > & Get_molalities(void)const       {return this->molalities;}
	inline const std::vector< std::pair< std::string, void * > > & Get_activities(void)const       {return this->activities;}
	inline const std::vector< std::pair< std::string, void * > > & Get_pure_phases(void)const      {return this->pure_phases;}
	inline const std::vector< std::pair< std::string, void * > > & Get_si(void)const               {return this->si;}
	inline const std::vector< std::pair< std::string, void * > > & Get_gases(void)const            {return this->gases;}
	inline const std::vector< std::pair< std::string, void * > > & Get_s_s(void)const              {return this->s_s;}
	inline const std::vector< std::pair< std::string, void * > > & Get_kinetics(void)const         {return this->kinetics;}
	inline const std::vector< std::pair< std::string, void * > > & Get_isotopes(void)const         {return this->isotopes;}
	inline const std::vector< std::pair< std::string, void * > > & Get_calculate_values(void)const {return this->calculate_values;}

	// file_name getters/setters
	void Set_file_name(int i);
	inline void Set_file_name(std::string s)                          {this->file_name = s;}
	inline std::string & Get_file_name(void)                          {return this->file_name;}
	inline const std::string & Get_file_name(void)const               {return this->file_name;}

	// punch_ostream getters/setters
	inline std::ostream* Get_punch_ostream(void)                      {return this->punch_ostream;}
	inline const std::ostream* Get_punch_ostream(void)const           {return this->punch_ostream;}
	inline void Set_punch_ostream(std::ostream * os)                  {this->punch_ostream = os;}

	// state var getters
	inline bool Get_active(void)const                                 {return this->active;}
	inline bool Get_new_def(void)const                                {return this->new_def;}
	inline bool Get_user_punch_new_def(void)const                     {return this->user_punch_new_def;}
	inline bool Get_have_punch_name(void)const                        {return this->have_punch_name;}

	// state var setters
	inline void Set_active(bool tf)                                   {this->active = tf;}
	inline void Set_new_def(bool tf)                                  {this->new_def = tf;}
	inline void Set_user_punch_new_def(bool tf)                       {this->user_punch_new_def = tf;}
	inline void Set_have_punch_name(bool tf)                          {this->have_punch_name = tf;}

	// as_is getters
	inline bool Get_user_punch(void)const                             {return this->user_punch;}
	inline bool Get_high_precision(void)const                         {return this->high_precision;}
	inline bool Get_inverse(void)const                                {return this->inverse;}

	inline bool Get_sim(void)const                                    {return this->sim;}
	inline bool Get_state(void)const                                  {return this->state;}
	inline bool Get_soln(void)const                                   {return this->soln;}
	inline bool Get_dist(void)const                                   {return this->dist;}
	inline bool Get_time(void)const                                   {return this->time;}

	inline bool Get_step(void)const                                   {return this->step;}
	inline bool Get_ph(void)const                                     {return this->ph;}
	inline bool Get_pe(void)const                                     {return this->pe;}
	inline bool Get_rxn(void)const                                    {return this->rxn;}
	inline bool Get_temp(void)const                                   {return this->temp;}

	inline bool Get_alk(void)const                                    {return this->alk;}
	inline bool Get_mu(void)const                                     {return this->mu;}
	inline bool Get_water(void)const                                  {return this->water;}
	inline bool Get_charge_balance(void)const                         {return this->charge_balance;}
	inline bool Get_percent_error(void)const                          {return this->percent_error;}
	inline bool Get_new_line(void)const                               {return this->new_line; }

	// as-is setters
	inline void Set_user_punch(bool tf)                               {this->user_punch = tf;              this->set_user_punch = true;}
	inline void Set_high_precision(bool tf)                           {this->high_precision = tf;          this->set_high_precision = true;}
	inline void Set_inverse(bool tf)                                  {this->inverse = tf;                 this->set_inverse = true;}

	inline void Set_sim(bool tf)                                      {this->sim = tf;                     this->set_sim = true;}
	inline void Set_state(bool tf)                                    {this->state = tf;                   this->set_state = true;}
	inline void Set_soln(bool tf)                                     {this->soln = tf;                    this->set_soln = true;}
	inline void Set_dist(bool tf)                                     {this->dist = tf;                    this->set_dist = true;}
	inline void Set_time(bool tf)                                     {this->time = tf;                    this->set_time = true;}

	inline void Set_step(bool tf)                                     {this->step = tf;                    this->set_step = true;}
	inline void Set_ph(bool tf)                                       {this->ph = tf;                      this->set_ph = true;}
	inline void Set_pe(bool tf)                                       {this->pe = tf;                      this->set_pe = true;}
	inline void Set_rxn(bool tf)                                      {this->rxn = tf;                     this->set_rxn = true;}
	inline void Set_temp(bool tf)                                     {this->temp = tf;                    this->set_temp = true;}

	inline void Set_alk(bool tf)                                      {this->alk = tf;                     this->set_alk = true;}
	inline void Set_mu(bool tf)                                       {this->mu = tf;                      this->set_mu = true;}
	inline void Set_water(bool tf)                                    {this->water = tf;                   this->set_water = true;}
	inline void Set_charge_balance(bool tf)                           {this->charge_balance = tf;          this->set_charge_balance = true;}
	inline void Set_percent_error(bool tf)                            {this->percent_error = tf;           this->set_percent_error = true;}
	inline void Set_new_line(bool tf)                                 {this->new_line = tf;                this->set_new_line = true;}

	// set flag getters
	inline bool was_set_user_punch()const                             {return this->set_user_punch;}
	inline bool was_set_high_precision()const                         {return this->set_high_precision;}
	inline bool was_set_inverse()const                                {return this->set_inverse;}

	inline bool was_set_sim()const                                    {return this->set_sim;}
	inline bool was_set_state()const                                  {return this->set_state;}
	inline bool was_set_soln()const                                   {return this->set_soln;}
	inline bool was_set_dist()const                                   {return this->set_dist;}
	inline bool was_set_time()const                                   {return this->set_time;}

	inline bool was_set_step()const                                   {return this->set_step;}
	inline bool was_set_ph()const                                     {return this->set_ph;}
	inline bool was_set_pe()const                                     {return this->set_pe;}
	inline bool was_set_rxn()const                                    {return this->set_rxn;}
	inline bool was_set_temp()const                                   {return this->set_temp;}

	inline bool was_set_alk()const                                    {return this->set_alk;}
	inline bool was_set_mu()const                                     {return this->set_mu;}
	inline bool was_set_water()const                                  {return this->set_water;}
	inline bool was_set_charge_balance()const                         {return this->set_charge_balance;}
	inline bool was_set_percent_error()const                          {return this->set_percent_error;}
	inline bool was_set_new_line()const                               {return this->set_new_line;}

protected:

	// vectors
	std::vector< std::pair< std::string, void * > > totals;
	std::vector< std::pair< std::string, void * > > molalities;
	std::vector< std::pair< std::string, void * > > activities;
	std::vector< std::pair< std::string, void * > > pure_phases;
	std::vector< std::pair< std::string, void * > > si;
	std::vector< std::pair< std::string, void * > > gases;
	std::vector< std::pair< std::string, void * > > s_s;
	std::vector< std::pair< std::string, void * > > kinetics;
	std::vector< std::pair< std::string, void * > > isotopes;
	std::vector< std::pair< std::string, void * > > calculate_values;

	// file_name
	std::string file_name;

	// punch_ostream
	std::ostream * punch_ostream;

	// state vars
	bool active;
	bool new_def;
	bool user_punch_new_def;
	bool have_punch_name;

	// as-is vars
	bool user_punch;
	bool high_precision;
	bool inverse;

	bool sim;
	bool state;
	bool soln;
	bool dist;
	bool time;

	bool step;
	bool ph;
	bool pe;
	bool rxn;
	bool temp;

	bool alk;
	bool mu;
	bool water;
	bool charge_balance;
	bool percent_error;
	bool new_line;

	// as-is set flags
	bool set_user_punch;
	bool set_high_precision;
	bool set_inverse;

	bool set_sim;
	bool set_state;
	bool set_soln;
	bool set_dist;
	bool set_time;

	bool set_step;
	bool set_ph;
	bool set_pe;
	bool set_rxn;
	bool set_temp;

	bool set_alk;
	bool set_mu;
	bool set_water;
	bool set_charge_balance;
	bool set_percent_error;
	bool set_new_line;
};
#endif // !defined(SELECTEDOUTPUT_H_INCLUDED)
