#include <assert.h>
#include "Utils.h"
#include "Phreeqc.h"
#include <istream>
#include <fstream>
#include "phqalloc.h"

/* ---------------------------------------------------------------------- */
void Phreeqc::
set_reading_database(int reading_database)
/* ---------------------------------------------------------------------- */
{
	reading_db = reading_database;
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
reading_database(void)
/* ---------------------------------------------------------------------- */
{
	return reading_db;
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
check_line(const char *string, int allow_empty, int allow_eof,
		   int allow_keyword, int print)
/* ---------------------------------------------------------------------- */
{
	if (reading_database())
		print = FALSE;
	return check_line_impl(string, allow_empty, allow_eof, allow_keyword,
						   print);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
check_line_impl(const char *string, int allow_empty, int allow_eof,
				int allow_keyword, int print)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function gets a new line and checks for empty, eof, and keywords.
 *
 *   Arguments:
 *      string        Input, character string used in printing error message
 *      allow_empty   Input, True or false, if a blank line is accepable
 *                       if false, another line is read
 *      allow_eof     Input, True or false, if EOF is acceptable
 *      allow_keyword Input, True or false, if a keyword is acceptable
 *
 *   Returns:
 *      EMPTY         if empty line read and allow_empty == true
 *      KEYWORD       if line begins with keyword
 *      EOF           if eof and allow_eof == true
 *      OK            otherwise
 *      OPTION        if line begins with -[alpha]
 *
 *   Terminates       if EOF and allow_eof == false.
 */
	int i;


/* Get line */
	do
	{
		i = get_line();
		if ((print == TRUE && i != EOF) || i == KEYWORD)
		{
			echo_msg(sformatf( "\t%s\n", line_save));
		}
	}
	while (i == EMPTY && allow_empty == FALSE);
/* Check eof */
	if (i == EOF && allow_eof == FALSE)
	{
		error_string = sformatf(
				"Unexpected eof while reading %s\nExecution terminated.\n",
				string);
		error_msg(error_string, STOP);
	}
/* Check keyword */
	if (i == KEYWORD && allow_keyword == FALSE)
	{
		error_string = sformatf(
				"Expected data for %s, but got a keyword ending data block.",
				string);
		error_msg(error_string, CONTINUE);
		input_error++;
	}
	check_line_return = i;
	return (i);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
get_line(void)
/* ---------------------------------------------------------------------- */
{
	PHRQ_io::LINE_TYPE j = phrq_io->get_line();
	// check_key sets next_keyword
	next_keyword = phrq_io->Get_m_next_keyword();

	// copy parser line to line and line_save
	// make sure there is enough space
	size_t l1 = strlen(phrq_io->Get_m_line().c_str()) + 1;
	size_t l2 = strlen(phrq_io->Get_m_line_save().c_str()) + 1;
	size_t l = (l1 > l2) ? l1 : l2;
	if (l >= (size_t) max_line)
	{
		max_line = (int) l * 2;
		line_save =	(char *) PHRQ_realloc(line_save,
			(size_t) max_line * sizeof(char));
		if (line_save == NULL)
			malloc_error();
		line = (char *) PHRQ_realloc(line, (size_t) max_line * sizeof(char));
		if (line == NULL)
			malloc_error();
	}
	strcpy(line, phrq_io->Get_m_line().c_str());
	strcpy(line_save, phrq_io->Get_m_line_save().c_str());
	return j;
}
