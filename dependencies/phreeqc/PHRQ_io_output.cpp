#include <assert.h>
#include "Phreeqc.h"
#include "phqalloc.h"

/* ---------------------------------------------------------------------- */
int Phreeqc::
warning_msg(const char *err_str)
/* ---------------------------------------------------------------------- */
{
	if (state == TRANSPORT && transport_warnings == FALSE)
		return (OK);
	if (state == ADVECTION && advection_warnings == FALSE)
		return (OK);
	count_warnings++;
	if (pr.warnings >= 0)
	{
		if (count_warnings > pr.warnings)
			return (OK);
	}
	if (phrq_io)
	{
		if (status_on)
		{
			phrq_io->screen_msg("\n");
		}
		std::ostringstream msg;
		msg << "WARNING: " << err_str;
		phrq_io->warning_msg(msg.str().c_str());
		status_on = false;
	}
	
	return OK;
}

/* ---------------------------------------------------------------------- */
void Phreeqc::
echo_msg(const char *str)
/* ---------------------------------------------------------------------- */
{
	if (pr.echo_input == TRUE)
	{
		if (phrq_io) phrq_io->echo_msg(str);
	}
}

/* ---------------------------------------------------------------------- */
void Phreeqc::
set_forward_output_to_log(int value)
/* ---------------------------------------------------------------------- */
{
	forward_output_to_log = value;
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
get_forward_output_to_log(void)
/* ---------------------------------------------------------------------- */
{
	return forward_output_to_log;
}

void Phreeqc::
fpunchf_heading(const char *name)
{
	if (pr.punch == TRUE && current_selected_output != NULL)
	{
		punch_msg(name);
	}
}
void Phreeqc::
fpunchf(const char *name, const char *format, double d)
{
	try
	{
		if (phrq_io) phrq_io->fpunchf(name, format, d);
	}
	catch(std::bad_alloc)
	{
		malloc_error();
	}
}
void Phreeqc::
fpunchf(const char *name, const char *format, char * s)
{
	try
	{
		if (phrq_io) phrq_io->fpunchf(name, format, s);
	}
	catch(std::bad_alloc)
	{
		malloc_error();
	}
}
void Phreeqc::
fpunchf(const char *name, const char *format, int d)
{
	try
	{
		if (phrq_io) phrq_io->fpunchf(name, format, d);
	}
	catch(std::bad_alloc)
	{
		malloc_error();
	}
}

void Phreeqc::
fpunchf_user(int user_index, const char *format, double d)
{
	const char *name;
	
	if (current_user_punch == NULL)
		return;
	// check headings
	//if (user_index < user_punch_count_headings)
	int user_punch_count_headings = (int) current_user_punch->Get_headings().size();
	if (user_index < user_punch_count_headings)
	{
		//name = user_punch_headings[user_index];
		name = current_user_punch->Get_headings()[user_index].c_str();
	}
	else
	{
		if (fpunchf_user_s_warning == 0)
		{
			error_string = sformatf(
					"USER_PUNCH: Headings count does not match number of calls to PUNCH.\n");
			warning_msg(error_string);
			fpunchf_user_s_warning = 1;
		}
		sprintf(fpunchf_user_buffer, "no_heading_%d",
				(user_index - user_punch_count_headings) + 1);
		name = fpunchf_user_buffer;
	}
	try
	{
		if (phrq_io) phrq_io->fpunchf(name, format, (double) d);
	}
	catch(std::bad_alloc)
	{
		malloc_error();
	}
}

void Phreeqc::
fpunchf_user(int user_index, const char *format, char * d)
{
	const char *name;
	
	if (current_user_punch == NULL)
		return;
	int user_punch_count_headings = (int) current_user_punch->Get_headings().size();
	// check headings
	if (user_index < user_punch_count_headings)
	{
		//name = user_punch_headings[user_index];
		name = current_user_punch->Get_headings()[user_index].c_str();
	}
	else
	{
		if (fpunchf_user_s_warning == 0)
		{
			error_string = sformatf(
					"USER_PUNCH: Headings count does not match number of calls to PUNCH.\n");
			warning_msg(error_string);
			fpunchf_user_s_warning = 1;
		}
		sprintf(fpunchf_user_buffer, "no_heading_%d",
				(user_index - user_punch_count_headings) + 1);
		name = fpunchf_user_buffer;
	}
	try
	{
		if (phrq_io) phrq_io->fpunchf(name, format, d);
	}
	catch(std::bad_alloc)
	{
		malloc_error();
	}
}

int Phreeqc::
fpunchf_end_row(const char *format)
{
	if (phrq_io) 
	{
		phrq_io->fpunchf_end_row(format);
	}
	return OK;
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
screen_msg(const char *err_str)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) phrq_io->screen_msg(err_str);
}
// ---------------------------------------------------------------------- */
// dump file methods
// ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
bool Phreeqc::
dump_open(const char *file_name)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io)
		return this->phrq_io->dump_open(file_name);
	return false;
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
dump_flush(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->dump_flush();
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
dump_close(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->dump_close();
}

/* ---------------------------------------------------------------------- */
void Phreeqc::
dump_msg(const char * str)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->dump_msg(str);
}
// ---------------------------------------------------------------------- */
// error file methods
// ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
bool Phreeqc::
error_open(const char *file_name)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io)
		return this->phrq_io->error_open(file_name);
	return false;
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
error_flush(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->error_flush();
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
error_close(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->error_close();
}

/* ---------------------------------------------------------------------- */
void Phreeqc::
error_msg(const char *err_str, bool stop)
/* ---------------------------------------------------------------------- */
{
	if (get_input_errors() <= 0)
		input_error = 1;
	if (phrq_io)
	{
		std::ostringstream msg;
		msg << "ERROR: " << err_str << "\n";

		phrq_io->output_msg(msg.str().c_str());
		phrq_io->log_msg(msg.str().c_str());

		if (status_on)
		{
			phrq_io->screen_msg("\n");
		}
		status_on = false;
		phrq_io->error_msg(msg.str().c_str(), stop);
	}

	if (stop)
	{
		throw PhreeqcStop();
	}
}
// ---------------------------------------------------------------------- */
// log file methods
// ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
bool Phreeqc::
log_open(const char *file_name)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io)
		return this->phrq_io->log_open(file_name);
	return false;
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
log_flush(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->log_flush();
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
log_close(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->log_close();
}

/* ---------------------------------------------------------------------- */
void Phreeqc::
log_msg(const char * str)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->log_msg(str);
}
// ---------------------------------------------------------------------- */
// output_temp file methods
// ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
bool Phreeqc::
output_open(const char *file_name)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) 
		return this->phrq_io->output_open(file_name);
	return false;
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
output_flush(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->output_flush();
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
output_close(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->output_close();
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
output_msg(const char * str)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io)
	{
		if (get_forward_output_to_log())
		{
			phrq_io->log_msg(str);
		}
		else
		{
			phrq_io->output_msg(str);
		}
	}
}
// ---------------------------------------------------------------------- */
// punch file methods
// ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
bool Phreeqc::
punch_open(const char *file_name, int n_user)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io)
		return this->phrq_io->punch_open(file_name, std::ios_base::out, n_user);
	return false;
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
punch_flush(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->punch_flush();
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
punch_close(void)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->punch_close();
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
punch_msg(const char * str)
/* ---------------------------------------------------------------------- */
{
	if (phrq_io) this->phrq_io->punch_msg(str);
}
