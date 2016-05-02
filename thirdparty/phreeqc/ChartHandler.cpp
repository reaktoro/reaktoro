// ChartHandler.cpp: implementation of the ChartHandler class.
//
//////////////////////////////////////////////////////////////////////
#if defined MULTICHART
#ifdef _DEBUG
#pragma warning(disable : 4786)	// disable truncation warning (Only used by debugger)
#endif
#include "ChartHandler.h"
#include "phreeqc.h"
#include <iostream>


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ChartHandler::ChartHandler(PHRQ_io *io)
:
PHRQ_base(io)
	//
	// default constructor for ChartHandler
	//
{
	current_chart = NULL;
	current_chart_n_user = -1000;
	u_g_defined = false;
	timer = true;
	active_charts = 0;
}

ChartHandler::~ChartHandler()
{
	std::map<int, ChartObject *>::iterator it;
	for (it = this->chart_map.begin(); it != chart_map.end(); it++)
	{
		delete it->second;
	}
}
void
ChartHandler::Punch_user_graph(Phreeqc * phreeqc_ptr)
{
	std::map<int, ChartObject *>::iterator it = this->chart_map.begin();
	for ( ; it != chart_map.end(); it++)
	{
		if (it->second->Get_active())
		{
#if defined(__cplusplus_cli)
			while (0 != System::Threading::Interlocked::CompareExchange(it->second->usingResource, 4, 0))
			{
				System::Threading::Thread::Sleep(5);
			}
#endif
			try
			{
				this->current_chart = it->second;
				phreeqc_ptr-> punch_user_graph();
			}
			catch (...)
			{
#if defined(__cplusplus_cli)
				int n = System::Threading::Interlocked::Exchange(it->second->usingResource, 0);
				assert(n == 4);
#endif
				throw;
			}
#if defined(__cplusplus_cli)
			System::Threading::Interlocked::Exchange(it->second->usingResource, 0);
#endif
		}
	}
}

bool
ChartHandler::Read(Phreeqc * phreeqc_ptr, CParser &parser)
{
	int n_user;
	std::string token;

	// reads line, next character is after keyword
	parser.check_line("ChartHandler", true, false, true, false);

	std::istringstream iss(parser.line());
	// keyword
	iss >> token;
	// number
	if (!(iss >> n_user))
	{
		n_user = 1;
	}

	// makes new ChartObject if necessary
	std::map<int, ChartObject *>::iterator it = this->chart_map.find(n_user);
	if (it == this->chart_map.end())
	{
		chart_map[n_user] = new ChartObject(this->Get_io());
		it = this->chart_map.find(n_user);
		it->second->Set_phreeqc(phreeqc_ptr);
	}

	// Read/update ChartObject
#if defined(__cplusplus_cli)
	while (0 != System::Threading::Interlocked::CompareExchange(it->second->usingResource, 5, 0))
	{
		System::Threading::Thread::Sleep(5);
	}
#endif
	try
	{
		{
			it->second->Read(parser);
			current_chart_n_user = n_user;
			current_chart = it->second;
			u_g_defined = true;
		}

		// if detached, set timer_end and free
		if (it->second->Get_detach() && it->second->Get_form_started())
		{
			it->second->Set_end_timer(true);
			it->second->Rate_free();
		}
	}
	catch(...)
	{
#if defined(__cplusplus_cli)
		// Release lock
		int n = System::Threading::Interlocked::Exchange(it->second->usingResource, 0);
		assert(n == 5);
		throw;
#endif
	}
#if defined(__cplusplus_cli)
	// Release lock
	int n = System::Threading::Interlocked::Exchange(it->second->usingResource, 0);
	assert(n == 5);
#endif

	// if detached, wait for thread to acknowledge and then erase chart
	if (it->second->Get_detach())
	{
		while (it->second->Get_form_started() && it->second->Get_done() != true) 
		{
#if defined(__cplusplus_cli)
			System::Threading::Thread::Sleep(5);
#endif
		}
		delete it->second;
		this->chart_map.erase(it);
	}	
	return true;
}
bool
ChartHandler::End_timer()
{
	
	size_t max_tries = 6000; // 1 h, but not used
	std::map<int, ChartObject *>::iterator it = this->chart_map.begin();
	if (chart_map.size() > 0) 
	{
		screen_msg("Detaching charts...");
		if (io != NULL)
		{
			io->error_flush();
		}
	}
	size_t i(0), i2(0);
	for  ( ; it != chart_map.end(); it++)
	{
		i = 0;
		it->second->Rate_free();
		if (it->second->Get_form_started())
		{
#if defined(__cplusplus_cli)
			while (0 != System::Threading::Interlocked::CompareExchange(it->second->usingResource, 6, 0))
			{
				//if (i > max_tries) break;
				i++;
				System::Threading::Thread::Sleep(60);
			}
#endif
			it->second->Set_end_timer(true);
			//it->second->Set_phreeqc(NULL);
#if defined(__cplusplus_cli)
			int n = System::Threading::Interlocked::Exchange(it->second->usingResource, 0);
			assert(n == 6);
#endif

			i2 = 0;
			while (it->second->Get_done() != true) 
			{
				//if (i2 > max_tries) break;
				i2++;
#if defined(__cplusplus_cli)
				System::Threading::Thread::Sleep(60);
#endif
			}
			//if (i >= max_tries || i2 >= max_tries)
			//{
			//	error_msg("\nChart did not respond.", CONTINUE);
			//}
		}
	}
	if (chart_map.size() > 0)
	{
		screen_msg("\rCharts detached.         \n");
		if (io != NULL)
		{
			io->error_flush();
		}
	}
	this->timer = false;

	return true;
}
bool
ChartHandler::dump(std::ostream & oss, unsigned int indent)
{
	std::map<int, ChartObject *>::iterator it = this->chart_map.begin();
	for  ( ; it != chart_map.end(); it++)
	{
		size_t i = 0;
		it->second->dump(oss, indent);
	}
	return true;
}
#endif

