#include "UserPunch.h"
#include "Phreeqc.h"
UserPunch::UserPunch(int n, PHRQ_io *io)
:	cxxNumKeyword(io)
{
	this->PhreeqcPtr    = NULL;
	this->rate          = NULL;
}


UserPunch::~UserPunch(void)
{
	if (this->rate != NULL)
	{
		this->PhreeqcPtr->rate_free(this->rate);
	}
	this->PhreeqcPtr->free_check_null(this->rate);
	this->rate = NULL;
}
