#include <assert.h>
#include "PHRQ_io.h"
#include "Parser.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <set>

#include <string.h>
#include <stdio.h>
#include <stdarg.h>

PHRQ_io::
PHRQ_io(void)
{
	output_ostream = NULL;
	log_ostream = NULL;
	punch_ostream = NULL;
#ifdef ERROR_OSTREAM
	error_ostream = NULL;
#else
	error_file = NULL;
#endif
	dump_ostream = NULL;
	io_error_count = 0;

	output_on = true;
	log_on = false;
	punch_on = true;
	error_on = true;
	dump_on = true;
	echo_on = true;
	screen_on = true;
	echo_destination = ECHO_OUTPUT;

	m_next_keyword = Keywords::KEY_NONE;
	accumulate = false;
	m_line_type = PHRQ_io::LT_EMPTY; 
}

PHRQ_io::
~PHRQ_io()
{
}
// ---------------------------------------------------------------------- */
// output ostream methods
// ---------------------------------------------------------------------- */

bool PHRQ_io::
ofstream_open(std::ostream **os, const char *file_name, std::ios_base::openmode mode)
{
	std::ofstream *ofs = new std::ofstream(file_name, mode);
	if (ofs && ofs->is_open())
	{
		*os = ofs;
		return true;
	}
	delete ofs;
	return false;
}

/* ---------------------------------------------------------------------- */
bool PHRQ_io::
output_open(const char *file_name, std::ios_base::openmode mode)
/* ---------------------------------------------------------------------- */
{
	return ofstream_open(&output_ostream, file_name, mode);
}
/* ---------------------------------------------------------------------- */
void PHRQ_io::
output_flush(void)
/* ---------------------------------------------------------------------- */
{
	if (output_ostream)
	{
		output_ostream->flush();
	}
}
/* ---------------------------------------------------------------------- */
void PHRQ_io::
output_close(void)
/* ---------------------------------------------------------------------- */
{
	safe_close(&output_ostream);
}
/* ---------------------------------------------------------------------- */
void PHRQ_io::
output_msg(const char * str)
/* ---------------------------------------------------------------------- */
{
	if (output_ostream != NULL && output_on)
	{
		(*output_ostream) << str;
	}
	//output_flush();
}
// ---------------------------------------------------------------------- */
// log ostream methods
// ---------------------------------------------------------------------- */
bool PHRQ_io::
log_open(const char *file_name, std::ios_base::openmode mode)
/* ---------------------------------------------------------------------- */
{
	return ofstream_open(&log_ostream, file_name, mode);
}
/* ---------------------------------------------------------------------- */
void PHRQ_io::
log_flush(void)
/* ---------------------------------------------------------------------- */
{
	if (log_ostream)
	{
		log_ostream->flush();
	}
}
/* ---------------------------------------------------------------------- */
void PHRQ_io::
log_close(void)
/* ---------------------------------------------------------------------- */
{
	safe_close(&log_ostream);
}
/* ---------------------------------------------------------------------- */
void PHRQ_io::
log_msg(const char * str)
/* ---------------------------------------------------------------------- */
{
	if (log_ostream != NULL && log_on)
	{
		(*log_ostream) << str;
	}
}
// ---------------------------------------------------------------------- */
// punch ostream methods
// ---------------------------------------------------------------------- */
bool PHRQ_io::
punch_open(const char *file_name, std::ios_base::openmode mode, int n_user)
/* ---------------------------------------------------------------------- */
{
	return ofstream_open(&punch_ostream, file_name, mode);
}
/* ---------------------------------------------------------------------- */
void PHRQ_io::
punch_flush(void)
/* ---------------------------------------------------------------------- */
{
	if (punch_ostream)
	{
		punch_ostream->flush();
	}
}
/* ---------------------------------------------------------------------- */
void PHRQ_io::
punch_close(void)
/* ---------------------------------------------------------------------- */
{
	safe_close(&punch_ostream);
}
/* ---------------------------------------------------------------------- */
void PHRQ_io::
punch_msg(const char * str)
/* ---------------------------------------------------------------------- */
{
	if (punch_ostream != NULL && punch_on)
	{
		(*punch_ostream) << str;
	}
}
// ---------------------------------------------------------------------- */
// error file methods
// ---------------------------------------------------------------------- */
#ifdef ERROR_OSTREAM
/* ---------------------------------------------------------------------- */
bool PHRQ_io::
error_open(const char *file_name, std::ios_base::openmode mode)
/* ---------------------------------------------------------------------- */
{
	if (file_name != NULL)
	{
		if (!ofstream_open(&error_ostream, file_name, mode))
		{
#if !defined(R_SO)
			error_ostream = &std::cerr;
#endif
			return false;
		}
	}
	else
	{
#if !defined(R_SO)
		error_ostream = &std::cerr;
#endif
	}
	return true;
}
/* ---------------------------------------------------------------------- */
void PHRQ_io::
error_flush(void)
/* ---------------------------------------------------------------------- */
{
	if (error_ostream)
	{
		error_ostream->flush();
	}
}
/* ---------------------------------------------------------------------- */
void PHRQ_io::
error_close(void)
/* ---------------------------------------------------------------------- */
{
	safe_close(&error_ostream);
}
/* ---------------------------------------------------------------------- */
void PHRQ_io::
error_msg(const char *err_str, bool stop)
/* ---------------------------------------------------------------------- */
{

	io_error_count++;
	if (error_ostream != NULL && error_on)
	{
		//(*error_ostream) << err_str;
		screen_msg(err_str);
		error_flush();
	}
	if (stop)
	{
		if (error_ostream != NULL && error_on)
		{
			//(*error_ostream) << "Stopping.\n";
			screen_msg("Stopping.\n");
			error_ostream->flush();
		}
		output_msg("Stopping.\n");
		log_msg("Stopping.\n");

		throw PhreeqcStop();
	}
}
/* ---------------------------------------------------------------------- */
void PHRQ_io::
warning_msg(const char *err_str)
/* ---------------------------------------------------------------------- */
{
	if (error_ostream != NULL && error_on)
	{
		//(*error_ostream) << err_str << "\n";
		std::string err_stdstr(err_str);
		err_stdstr.append("\n");
		screen_msg(err_stdstr.c_str());
		error_ostream->flush();
	}
	std::ostringstream warn_str;
	warn_str << err_str << "\n";
	log_msg(warn_str.str().c_str());
	log_flush();
	output_msg(warn_str.str().c_str());
	output_flush();
}
/* ---------------------------------------------------------------------- */
void PHRQ_io::
screen_msg(const char * str)
/* ---------------------------------------------------------------------- */
{
	if (error_ostream != NULL && screen_on)
	{
		(*error_ostream) << str;
	}
}
#else
/* ---------------------------------------------------------------------- */
bool PHRQ_io::
error_open(const char *file_name, const char * mode)
/* ---------------------------------------------------------------------- */
{
	if (file_name != NULL)
	{
		if ((error_file = fopen(file_name, mode)) == NULL)
		{
			error_file = stderr;
			return false;
		}
	}
	else
	{
		error_file = stderr;
	}
	return true;
}
/* ---------------------------------------------------------------------- */
void PHRQ_io::
error_flush(void)
/* ---------------------------------------------------------------------- */
{
	if (error_file)
	{
		fflush(error_file);
	}
}
/* ---------------------------------------------------------------------- */
void PHRQ_io::
error_close(void)
/* ---------------------------------------------------------------------- */
{
	safe_close(&error_file);
}
/* ---------------------------------------------------------------------- */
void PHRQ_io::
error_msg(const char *err_str, bool stop)
/* ---------------------------------------------------------------------- */
{

	io_error_count++;
	if (error_file != NULL && error_on)
	{
		//(*error_file) << err_str;
		screen_msg(err_str);
		error_flush();
	}
	if (stop)
	{
		if (error_file != NULL && error_on)
		{
			//(*error_file) << "Stopping.\n";
			screen_msg("Stopping.\n");
			fflush(error_file);
		}
		output_msg("Stopping.\n");
		log_msg("Stopping.\n");
	}
}
/* ---------------------------------------------------------------------- */
void PHRQ_io::
warning_msg(const char *err_str)
/* ---------------------------------------------------------------------- */
{
	if (error_file != NULL && error_on)
	{
		//(*error_file) << err_str << "\n";
		std::string err_stdstr(err_str);
		err_stdstr.append("\n");
		screen_msg(err_stdstr.c_str());
		error_flush();
	}
	std::ostringstream warn_str;
	warn_str << err_str << "\n";
	log_msg(warn_str.str().c_str());
	log_flush();
	output_msg(warn_str.str().c_str());
	output_flush();
}
/* ---------------------------------------------------------------------- */
void PHRQ_io::
screen_msg(const char * str)
/* ---------------------------------------------------------------------- */
{
	std::string stdstr(str);
	if (error_file != NULL && screen_on)
	{
		fprintf(error_file, "%s", str);
	}
}
#endif
// ---------------------------------------------------------------------- */
// dump ostream methods
// ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
bool PHRQ_io::
dump_open(const char *file_name, std::ios_base::openmode mode)
/* ---------------------------------------------------------------------- */
{
	return ofstream_open(&dump_ostream, file_name, mode);
}
/* ---------------------------------------------------------------------- */
void PHRQ_io::
dump_flush(void)
/* ---------------------------------------------------------------------- */
{
	if (dump_ostream)
	{
		dump_ostream->flush();
	}
}
/* ---------------------------------------------------------------------- */
void PHRQ_io::
dump_close(void)
/* ---------------------------------------------------------------------- */
{
	safe_close(&dump_ostream);
}
/* ---------------------------------------------------------------------- */
void PHRQ_io::
dump_msg(const char * str)
/* ---------------------------------------------------------------------- */
{
	if (dump_ostream != NULL && dump_on)
	{
		(*dump_ostream) << str;
	}
}
int PHRQ_io::
getc(void)
{
	if (std::istream* is = get_istream())
	{
		int n = is->get();
		if (n == 13 && is->peek() == 10)
		{
			n = is->get();
		}
		return n;
	}
	return EOF;
}

/* ---------------------------------------------------------------------- */
void PHRQ_io::
fpunchf(const char *name, const char *format, double d)
/* ---------------------------------------------------------------------- */
{
	if (punch_ostream != NULL && punch_on)
	{
		fpunchf_helper(punch_ostream, format, d);
	}
}
/* ---------------------------------------------------------------------- */
void PHRQ_io::
fpunchf(const char *name, const char *format, char * s)
/* ---------------------------------------------------------------------- */
{
	if (punch_ostream != NULL && punch_on)
	{
		fpunchf_helper(punch_ostream, format, s);
	}
}
/* ---------------------------------------------------------------------- */
void PHRQ_io::
fpunchf(const char *name, const char *format, int d)
/* ---------------------------------------------------------------------- */
{
	if (punch_ostream != NULL && punch_on)
	{
		fpunchf_helper(punch_ostream, format, d);
	}
}
/* ---------------------------------------------------------------------- */
void PHRQ_io::
fpunchf_helper(std::ostream *os, const char *format, ...)
/* ---------------------------------------------------------------------- */
{
	if (os)
	{
		const size_t STACK_MAX = 2048;
		char stack_buffer[STACK_MAX];

		va_list args;
		va_start(args, format);
		int j = ::vsnprintf(stack_buffer, STACK_MAX, format, args);
		bool success = (j >= 0 && j < (int) STACK_MAX);
		va_end(args);

		if (success)
		{
			(*os) << stack_buffer;
		}
		else
		{
			size_t alloc_buffer_size = STACK_MAX * 2;
			char *alloc_buffer = new char[alloc_buffer_size];
			do 
			{
				va_list args;
				va_start(args, format);
				j = ::vsnprintf(alloc_buffer, alloc_buffer_size, format, args);
				success = (j >= 0 && j < (int) alloc_buffer_size);
				va_end(args);
				if (!success)
				{
					delete[] alloc_buffer;
					alloc_buffer_size *= 2;
					alloc_buffer = new char[alloc_buffer_size];
				}
			}
			while (!success);

			(*os) << alloc_buffer;
			delete[] alloc_buffer;
		}
	}
}
/* ---------------------------------------------------------------------- */
void PHRQ_io::
fpunchf_helper(std::string *str, const char *format, ...)
/* ---------------------------------------------------------------------- */
{
	if (str)
	{
		const size_t STACK_MAX = 2048;
		char stack_buffer[STACK_MAX];

		va_list args;
		va_start(args, format);
		int j = ::vsnprintf(stack_buffer, STACK_MAX, format, args);
		bool success = (j >= 0 && j < (int) STACK_MAX);
		va_end(args);

		if (success)
		{
			(*str) += stack_buffer;
		}
		else
		{
			size_t alloc_buffer_size = STACK_MAX * 2;
			char *alloc_buffer = new char[alloc_buffer_size];
			do 
			{
				va_list args;
				va_start(args, format);
				j = ::vsnprintf(alloc_buffer, alloc_buffer_size, format, args);
				success = (j >= 0 && j < (int) alloc_buffer_size);
				va_end(args);
				if (!success)
				{
					delete[] alloc_buffer;
					alloc_buffer_size *= 2;
					alloc_buffer = new char[alloc_buffer_size];
				}
			}
			while (!success);

			(*str) += alloc_buffer;
			delete[] alloc_buffer;
		}
	}
}

/* ---------------------------------------------------------------------- */
void PHRQ_io::fpunchf_end_row(const char *format)
/* ---------------------------------------------------------------------- */
{
	//NOOP for Phreeqc
}
/* ---------------------------------------------------------------------- */
void PHRQ_io::
close_ostreams(void)
/* ---------------------------------------------------------------------- */
{
	std::set<std::ostream *> streams;

	streams.insert(output_ostream);
	streams.insert(log_ostream);
//	streams.insert(punch_ostream);   // Should be deleted in ~SelectedOutput
#ifdef ERROR_OSTREAM
	streams.insert(error_ostream);
#else
	safe_close(&error_file);
#endif
	streams.insert(dump_ostream);

	std::set<std::ostream *>::iterator it = streams.begin();
	for (; it != streams.end(); it++)
	{
		std::ostream * x = *it;
		safe_close(&x);
	}

	output_ostream = NULL;
	log_ostream = NULL;
	punch_ostream = NULL;
#ifdef ERROR_OSTREAM
	error_ostream = NULL;
#else
	error_file = NULL;
#endif
	dump_ostream = NULL;
}
//safe_close is static method
/* ---------------------------------------------------------------------- */
void PHRQ_io::
safe_close(std::ostream **stream_ptr)
/* ---------------------------------------------------------------------- */
{
	if (
#if !defined(R_SO)
		*stream_ptr != &std::cerr &&
		*stream_ptr != &std::cout &&
		*stream_ptr != &std::clog &&
#endif
		*stream_ptr != NULL)
	{
		delete *stream_ptr;
		*stream_ptr = NULL;
	}
}
void PHRQ_io::
safe_close(FILE **file_ptr)
/* ---------------------------------------------------------------------- */
{
	if (
#if !defined(R_SO)
		*file_ptr != stderr &&
		*file_ptr != stdout &&
		*file_ptr != stdin &&
#endif
		*file_ptr != NULL)
	{
		fclose(*file_ptr);
		*file_ptr = NULL;
	}
}
/* ---------------------------------------------------------------------- */
void PHRQ_io::
echo_msg(const char * str)
/* ---------------------------------------------------------------------- */
{
	if (echo_on)
	{
		switch (this->echo_destination)
		{
		case ECHO_LOG:
			log_msg(str);
			break;
		case ECHO_OUTPUT:
			output_msg(str);
			break;
		}
	}
}

std::istream * PHRQ_io::
get_istream()
{
	if (istream_list.size() > 0)
	{
		return istream_list.front();
	}
	else
	{
		return NULL;
	}
}
void PHRQ_io::
push_istream(std::istream * cookie, bool auto_delete)
{
	istream_list.push_front(cookie);
	delete_istream_list.push_front(auto_delete);
}
void PHRQ_io::
clear_istream(void)
{
	while (istream_list.size() > 0)
	{
		pop_istream();
	}
}
void PHRQ_io::
pop_istream()
{
	if (istream_list.size() > 0)
	{
		if (delete_istream_list.front())
		{
			delete istream_list.front();
		}
		istream_list.pop_front();
		delete_istream_list.pop_front();
	}
}
/* ---------------------------------------------------------------------- */
PHRQ_io::LINE_TYPE PHRQ_io::
get_line(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Read a line from input file put in "line".
 *   Copy of input line is stored in "line_save".
 *   Characters after # are discarded in line but retained in "line_save"
 *
 *   Arguments:
 *      fp is file name
 *   Returns:
 *      EMPTY,
 *      EOF,
 *      KEYWORD,
 *      OK,
 *      OPTION
 */
	std::string stdtoken;
	bool continue_loop = true;;

	PHRQ_io::LINE_TYPE return_value;
	// loop for include files
	for (;;)
	{
		if (this->get_istream() == NULL)
		{
			break;
		}
		return_value = LT_EMPTY;
		while (return_value == LT_EMPTY)
		{
			/*
			*   Eliminate all characters after # sign as a comment
			*/
			/*
			*   Get line, check for eof
			*/
			continue_loop = false;

			if (get_logical_line() == LT_EOF)
			{
				//pop next file
				this->pop_istream();
				continue_loop = true;
				break;
			}
			/*
			*   Get long lines
			*/
			bool empty = true;
			m_line = m_line_save.substr(0, m_line_save.find_first_of('#'));
			for (unsigned int i = 0; i < m_line.size(); ++i)
			{
				if (!::isspace(m_line[i]))
				{
					empty = false;
					break;
				}
			}

			if (this->accumulate)
			{
				this->accumulated.append(m_line_save);
				this->accumulated.append("\n");
			}
			//
			// New line character encountered
			//
			return_value = (empty ? LT_EMPTY : LT_OK);

		}
		if (continue_loop) continue;
		//
		// Determine return_value
		//
		if (return_value == LT_OK)
		{
			if (check_key(m_line.begin(), m_line.end()))
			{
				return_value = LT_KEYWORD;
			}
			else
			{
				std::string::iterator beg = m_line.begin();
				std::string::iterator end = m_line.end();
				std::string token;
				CParser::copy_token(token, beg, end);

				if (token.size() > 1 && token[0] == '-' &&::isalpha(token[1]))
				{
					return_value = LT_OPTION;
				}
			}
		}

		// add new include file to stack
		std::string::iterator beg = m_line.begin();
		std::string::iterator end = m_line.end();
		CParser::copy_token(stdtoken, beg, end);
		std::transform(stdtoken.begin(), stdtoken.end(), stdtoken.begin(), ::tolower);
		if ((strstr(stdtoken.c_str(),"include$") == stdtoken.c_str()) ||
			(strstr(stdtoken.c_str(),"include_file") == stdtoken.c_str()))
		{
			std::string file_name;
			file_name.assign(beg, end);
			file_name = trim(file_name);

			if (file_name.size() > 0)
			{
				std::ifstream *next_stream = new std::ifstream(file_name.c_str(), std::ios_base::in);
				if (!next_stream->is_open())
				{
					std::ostringstream errstr;
					errstr << "\n***********  Could not open include file " << file_name 
						   <<".\n             Please, write the full path to this file. ***********\n\n";
					delete next_stream;
#if defined(PHREEQCI_GUI)
					warning_msg(errstr.str().c_str());
					continue;
#else
					output_msg(errstr.str().c_str());
					error_msg(errstr.str().c_str(), OT_STOP);
#endif
				}
				else
				{
					this->push_istream(next_stream);
				}
				continue;
			}
		}
		return return_value;
	}
	m_next_keyword = Keywords::KEY_END;
	return LT_EOF;
}

/**
        Reads input stream until end of line, ";", or eof
        stores characters in line_save

        returns:
                EOF on empty line on end of file or
                OK otherwise
*/
PHRQ_io::LINE_TYPE PHRQ_io::
get_logical_line(void)
{
	int j;
	unsigned int pos;
	char c;

	m_line_save.erase(m_line_save.begin(), m_line_save.end());	// m_line_save.clear();
	while ((j = getc()) != EOF)
	{
		c = (char) j;
		if (c == '#')
		{
			// ignore all chars after # until newline
			do
			{
				c = (char) j;
				if (c == '\n')
				{
					break;
				}
				m_line_save += c;
			}
			while ((j = getc()) != EOF);
		}
		if (c == ';')
			break;
		if (c == '\n')
		{
			break;
		}
		if (c == '\\')
		{
			pos = (int) m_line_save.size();
			m_line_save += c;
			while ((j = getc()) != EOF)
			{
				c = (char) j;
				if (c == '\\')
				{
					pos = (int) m_line_save.size();
					m_line_save += c;
					continue;
				}
				if (c == '\n')
				{
					// remove '\\'
					m_line_save = m_line_save.substr(0,pos);
					break;
				}
				m_line_save += c;
				if (!::isspace(j))
					break;
			}
		}
		else
		{
			m_line_save += c;
		}
	}
	if (j == std::char_traits < char >::eof() && m_line_save.size() == 0)
	{
		return (LT_EOF);
	}
	return (LT_OK);
}

bool PHRQ_io::
check_key(std::string::iterator begin, std::string::iterator end)
{
	std::string lowercase;
	CParser::copy_token(lowercase, begin, end);
	std::transform(lowercase.begin(), lowercase.end(), lowercase.begin(),
				   tolower);

	m_next_keyword = Keywords::Keyword_search(lowercase);
	if (m_next_keyword == Keywords::KEY_NONE)
	{
		return false;
	}
	return true;
}

