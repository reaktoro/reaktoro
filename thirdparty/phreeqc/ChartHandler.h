#if !defined(CHARTHANDLER_H_INCLUDED)
#define CHARTHANDLER_H_INCLUDED
#if defined MULTICHART
#include <vector>
#include <map>
#include <string>
#include "Parser.h"
#include "ChartObject.h"
#include "PHRQ_base.h"

class ChartHandler: public PHRQ_base
{

public:
	ChartHandler(PHRQ_io *io = NULL);
	virtual ~ChartHandler();

	size_t Get_chart_count()const
	{
		return this->chart_map.size();
	}
	ChartObject * Get_current_chart()
	{
		return this->current_chart;
	}
	const ChartObject * Get_current_chart()const
	{
		return this->current_chart;
	}
	bool Get_timer()
	{
		return timer;
	}
	int Get_active_charts() {return this->active_charts;}
	void Increment_active_charts()
	{
		System::Threading::Interlocked::Increment(this->active_charts);
	}
	void Decrement_active_charts()
	{
		System::Threading::Interlocked::Decrement(this->active_charts);
	}
	bool Read(Phreeqc * phreeqc_ptr, CParser &parser);
	void Punch_user_graph(Phreeqc * phreeqc_ptr);
	bool End_timer();
	bool dump(std::ostream & oss, unsigned int indent);
protected:
	std::map<int, ChartObject *> chart_map;
	int current_chart_n_user;
	ChartObject * current_chart;
	bool u_g_defined;
	bool timer;
	int active_charts;

public:

};
#endif // MULTICHART
#endif // !defined(CHARTHANDLER_H_INCLUDED)
