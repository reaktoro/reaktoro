#if !defined(RUNNER_H_INCLUDED)
#define RUNNER_H_INCLUDED
#include <set>					// std::set
#include <string>				// std::string

#include "phrqtype.h"
#include "StorageBinList.h"
#include "PHRQ_base.h"
class CParser;

class runner: public PHRQ_base
{
public:
	runner(PHRQ_io *io=NULL);
	runner(CParser & parser, PHRQ_io *io=NULL);
	virtual ~runner(void);
	bool Read(CParser & parser);
	StorageBinListItem & Get_cells(void) { return(this->cells); };
	LDBLE Get_time_step() { return(this->time_step); };
	LDBLE Get_start_time() { return(this->start_time); };
	void Set_time_step(LDBLE ts) { this->time_step = ts; };
	void Set_start_time(LDBLE st) { this->start_time = st; };
	bool Get_run_cells() { return(this->run_cells); };
	void Set_run_cells(bool tf) { this->run_cells = tf; };

protected:
	LDBLE time_step;
	LDBLE start_time;
	StorageBinListItem cells;
	bool run_cells;
	const static std::vector < std::string > vopts;
};
#endif // !defined(RUNNER_H_INCLUDED)
