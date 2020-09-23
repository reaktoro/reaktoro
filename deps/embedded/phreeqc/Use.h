#if !defined(USE_H_INCLUDED)
#define USE_H_INCLUDED
class cxxPPassemblage;
class cxxMix;
class cxxReaction;
class cxxExchange;
class cxxGasPhase;
class cxxPressure;
class cxxTemperature;
class cxxSSassemblage;
class cxxKinetics;
class cxxSurface;
class cxxSolution;

class cxxUse
{
public:
	cxxUse();
	virtual ~cxxUse(void);
	void init(void);

	bool Get_solution_in(void) const      {return this->solution_in;}
	bool Get_pp_assemblage_in(void) const {return this->pp_assemblage_in;}
	bool Get_mix_in(void) const           {return this->mix_in;}
	bool Get_reaction_in(void) const      {return this->reaction_in;}
	bool Get_exchange_in(void) const      {return this->exchange_in;}
	bool Get_kinetics_in(void) const      {return this->kinetics_in;}
	bool Get_surface_in(void) const       {return this->surface_in;}
	bool Get_pressure_in(void) const      {return this->pressure_in;}
	bool Get_temperature_in(void) const   {return this->temperature_in;}
	bool Get_gas_phase_in(void) const     {return this->gas_phase_in;}
	bool Get_inverse_in(void) const       {return this->inverse_in;}
	bool Get_ss_assemblage_in(void) const {return this->ss_assemblage_in;}
	bool Get_advect_in(void) const        {return this->advect_in;}
	bool Get_trans_in(void) const         {return this->trans_in;}

	void Set_solution_in(bool tf)         {solution_in = tf;}
	void Set_pp_assemblage_in(bool tf)    {pp_assemblage_in = tf;}
	void Set_mix_in(bool tf)              {mix_in = tf;}
	void Set_reaction_in(bool tf)         {reaction_in = tf;}
	void Set_exchange_in(bool tf)         {exchange_in = tf;}
	void Set_kinetics_in(bool tf)         {kinetics_in = tf;}
	void Set_surface_in(bool tf)          {surface_in = tf;}
	void Set_pressure_in(bool tf)         {pressure_in = tf;}
	void Set_temperature_in(bool tf)      {temperature_in = tf;}
	void Set_gas_phase_in(bool tf)        {gas_phase_in = tf;}
	void Set_inverse_in(bool tf)          {inverse_in = tf;}
	void Set_ss_assemblage_in(bool tf)    {ss_assemblage_in = tf;}
	void Set_advect_in(bool tf)           {advect_in = tf;}
	void Set_trans_in(bool tf)            {trans_in = tf;}

	int Get_n_solution_user(void) const       {return n_solution_user;}
	int Get_n_pp_assemblage_user(void) const  {return n_pp_assemblage_user;}
	int Get_n_mix_user(void) const            {return n_mix_user;}
	int Get_n_mix_user_orig(void) const       {return n_mix_user_orig;}
	int Get_n_reaction_user(void) const       {return n_reaction_user;}
	int Get_n_exchange_user(void) const       {return n_exchange_user;}
	int Get_n_kinetics_user(void) const       {return n_kinetics_user;}
	int Get_n_surface_user(void) const        {return n_surface_user;}
	int Get_n_pressure_user(void) const       {return n_pressure_user;}
	int Get_n_temperature_user(void) const    {return n_temperature_user;}
	int Get_n_gas_phase_user(void) const      {return n_gas_phase_user;}
	int Get_n_inverse_user(void) const        {return n_inverse_user;}
	int Get_n_ss_assemblage_user(void) const  {return n_ss_assemblage_user;}

	void Set_n_solution_user(int i)       {n_solution_user = i;}
	void Set_n_pp_assemblage_user(int i)  {n_pp_assemblage_user = i;}
	void Set_n_mix_user(int i)            {n_mix_user = i;}
	void Set_n_mix_user_orig(int i)       {n_mix_user_orig = i;}
	void Set_n_reaction_user(int i)       {n_reaction_user = i;}
	void Set_n_exchange_user(int i)       {n_exchange_user = i;}
	void Set_n_kinetics_user(int i)       {n_kinetics_user = i;}
	void Set_n_surface_user(int i)        {n_surface_user = i;}
	void Set_n_pressure_user(int i)       {n_pressure_user = i;}
	void Set_n_temperature_user(int i)    {n_temperature_user = i;}
	void Set_n_gas_phase_user(int i)      {n_gas_phase_user = i;}
	void Set_n_inverse_user(int i)        {n_inverse_user = i;}
	void Set_n_ss_assemblage_user(int i)  {n_ss_assemblage_user = i;}

	cxxSolution * Get_solution_ptr(void) const               {return this->solution_ptr;}
	cxxPPassemblage * Get_pp_assemblage_ptr(void) const      {return this->pp_assemblage_ptr;}
	cxxMix * Get_mix_ptr(void) const                         {return this->mix_ptr;}
	cxxReaction * Get_reaction_ptr(void) const               {return this->reaction_ptr;}
	cxxExchange * Get_exchange_ptr(void) const               {return this->exchange_ptr;}
	cxxKinetics * Get_kinetics_ptr(void) const               {return this->kinetics_ptr;}
	cxxSurface * Get_surface_ptr(void) const                 {return this->surface_ptr;}
	cxxPressure * Get_pressure_ptr(void) const               {return this->pressure_ptr;}
	cxxTemperature * Get_temperature_ptr(void) const         {return this->temperature_ptr;}
	cxxGasPhase * Get_gas_phase_ptr(void) const              {return this->gas_phase_ptr;}
	struct inverse * Get_inverse_ptr(void) const             {return this->inverse_ptr;}
	cxxSSassemblage * Get_ss_assemblage_ptr(void)            {return this->ss_assemblage_ptr;}

	void Set_solution_ptr(cxxSolution * p)                   {this->solution_ptr = p;}
	void Set_pp_assemblage_ptr(cxxPPassemblage * p)          {this->pp_assemblage_ptr = p;}
	void Set_mix_ptr(cxxMix * p)                             {this->mix_ptr = p;}
	void Set_reaction_ptr(cxxReaction * p)                   {this->reaction_ptr = p;}
	void Set_exchange_ptr(cxxExchange * p)                   {this->exchange_ptr = p;}
	void Set_kinetics_ptr(cxxKinetics * p)                   {this->kinetics_ptr = p;}
	void Set_surface_ptr(cxxSurface * p)                     {this->surface_ptr = p;}
	void Set_pressure_ptr(cxxPressure * p)                   {this->pressure_ptr = p;}
	void Set_temperature_ptr(cxxTemperature * p)             {this->temperature_ptr = p;}
	void Set_gas_phase_ptr(cxxGasPhase * p)                  {this->gas_phase_ptr = p;}
	void Set_inverse_ptr(struct inverse * p)                 {this->inverse_ptr = p;}
	void Set_ss_assemblage_ptr(cxxSSassemblage * p)     {this->ss_assemblage_ptr = p;}

protected:
	bool solution_in;
	int n_solution_user;
	cxxSolution *solution_ptr;

	bool pp_assemblage_in;
	int n_pp_assemblage_user;
	cxxPPassemblage *pp_assemblage_ptr;

	bool mix_in;
	int n_mix_user;
	cxxMix * mix_ptr;
	int n_mix_user_orig;

	bool reaction_in;
	int n_reaction_user;
	cxxReaction * reaction_ptr;

	bool exchange_in;
	int n_exchange_user;
	cxxExchange * exchange_ptr;

	bool kinetics_in;
	int n_kinetics_user;
	cxxKinetics *kinetics_ptr;

	bool surface_in;
	int n_surface_user;
	cxxSurface *surface_ptr;

	bool pressure_in;
	int n_pressure_user;
	cxxPressure *pressure_ptr;

	bool temperature_in;
	int n_temperature_user;
	cxxTemperature *temperature_ptr;

	bool inverse_in;
	int n_inverse_user;
	struct inverse *inverse_ptr;

	bool gas_phase_in;
	int n_gas_phase_user;
	cxxGasPhase * gas_phase_ptr;

	bool ss_assemblage_in;
	int n_ss_assemblage_user;
	cxxSSassemblage *ss_assemblage_ptr;

	bool trans_in;
	bool advect_in;
};
#endif // !defined(USE_H_INCLUDED)
