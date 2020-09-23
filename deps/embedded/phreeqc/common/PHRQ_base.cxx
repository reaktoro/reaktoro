#include "PHRQ_base.h"
#include <iostream>
#include "PHRQ_io.h"
PHRQ_base::
PHRQ_base(void)
{
	this->io = NULL;
	base_error_count = 0;
}

PHRQ_base::
PHRQ_base(PHRQ_io * p_io)
{
	this->io = p_io;
	base_error_count = 0;
}

PHRQ_base::
~PHRQ_base()
{
}

void PHRQ_base::
error_msg(const std::string & stdstr, int stop)
{
	this->base_error_count++;
	std::ostringstream msg;
	msg << "ERROR: " << stdstr << "\n";
	if (this->io)
	{
		this->io->output_msg(msg.str().c_str());
		this->io->log_msg(msg.str().c_str());
		this->io->error_msg("\n");
		this->io->error_msg(msg.str().c_str(), stop!=0);
	}
	else
	{
#if !defined(R_SO)
		std::cerr << msg.str().c_str();
		std::cout << msg.str().c_str();
#endif
	}
	if (stop != 0)
	{
		throw PhreeqcStop();
	}
}

void PHRQ_base::
warning_msg(const std::string & stdstr)
{
	if (this->io)
	{
		this->io->warning_msg(stdstr.c_str());
	}
	else
	{
#if !defined(R_SO)
		std::cerr << stdstr << "\n";
		std::cout << stdstr << "\n";
#endif
	}
}

void PHRQ_base::
output_msg(const std::string & stdstr)
{
	if (this->io)
	{
		this->io->output_msg(stdstr.c_str());
	}
	else
	{
#if !defined(R_SO)
		std::cout << stdstr << "\n";
#endif
	}
}

void PHRQ_base::
screen_msg(const std::string & stdstr)
{
	if (this->io)
	{
		this->io->screen_msg(stdstr.c_str());
	}
	else
	{
#if !defined(R_SO)
		std::cerr << stdstr << "\n";
#endif
	}
}

void PHRQ_base::
echo_msg(const std::string & stdstr)
{
	if (this->io)
	{
		this->io->echo_msg(stdstr.c_str());
	}
	else
	{
#if !defined(R_SO)
		std::cout << stdstr << "\n";
#endif
	}
}