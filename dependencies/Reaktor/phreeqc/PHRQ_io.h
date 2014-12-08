#ifndef _PHRQIO_H
#define _PHRQIO_H

#if defined(_WINDLL)
#define IPQ_DLL_EXPORT __declspec(dllexport)
#else
#define IPQ_DLL_EXPORT
#endif

#include <iostream>
#include <exception>
#include <list>
#include "Keywords.h"
#include <time.h>

#define ERROR_OSTREAM

class PhreeqcStop : public std::exception
{
};

class IPQ_DLL_EXPORT PHRQ_io
{
public:
	enum LINE_TYPE
	{
		LT_EOF = -1,
		LT_OK = 1,
		LT_EMPTY = 2,
		LT_KEYWORD = 3,
		LT_OPTION = 8
	};

	enum ONERROR_TYPE
	{
		OT_CONTINUE = 0,
		OT_STOP = 1
	};
	// constructor/destructor
	PHRQ_io(void);
	virtual ~ PHRQ_io();

	// methods
	static void safe_close(std::ostream **stream_ptr);
	static void safe_close(FILE **file_ptr);
	void close_ostreams(void);
	void Set_io_error_count(int i)				{this->io_error_count = i;};
	int Get_io_error_count(void)				{return this->io_error_count;};


	// istreams
	std::istream *get_istream();
	void pop_istream();
	void push_istream(std::istream * cookie, bool auto_delete = true);
	void clear_istream(void);

	// helper
	bool ofstream_open(std::ostream **os, const char *file_name, std::ios_base::openmode mode = std::ios_base::out);

	// output_ostream
	virtual bool output_open(const char *file_name, std::ios_base::openmode mode = std::ios_base::out);
	void output_flush(void);
	void output_close(void);
	virtual void output_msg(const char * str);
	void Set_output_ostream(std::ostream * out)		{this->output_ostream = out;};
	std::ostream *Get_output_ostream(void)			{return this->output_ostream;};
	void Set_output_on(bool tf)						{this->output_on = tf;};
	bool Get_output_on(void)						{return this->output_on;};

	// log_ostream
	virtual bool log_open(const char *file_name, std::ios_base::openmode mode = std::ios_base::out);
	void log_flush(void);
	void log_close(void);
	virtual void log_msg(const char * str);
	void Set_log_ostream(std::ostream * out)		{this->log_ostream = out;}
	std::ostream *Get_log_ostream(void)				{return this->log_ostream;}
	void Set_log_on(bool tf)						{this->log_on = tf;}
	bool Get_log_on(void)							{return this->log_on;}

	// punch_ostream
	virtual bool punch_open(const char *file_name, std::ios_base::openmode mode = std::ios_base::out, int n_user = 1);
	void punch_flush(void);
	void punch_close(void);
	virtual void punch_msg(const char * str);
	void Set_punch_ostream(std::ostream * out)		{this->punch_ostream = out;}
	std::ostream *Get_punch_ostream(void)			{return this->punch_ostream;}
	void Set_punch_on(bool tf)						{this->punch_on = tf;}
	bool Get_punch_on(void)							{return this->punch_on;}
	
	// error_ostream
#ifdef ERROR_OSTREAM
	virtual bool error_open(const char *file_name, std::ios_base::openmode mode = std::ios_base::out);
	void error_flush(void);
	void error_close(void);
	virtual void error_msg(const char * str, bool stop=false);
	void Set_error_ostream(std::ostream * out)		{this->error_ostream = out;}
	std::ostream *Get_error_ostream(void)			{return this->error_ostream;}
	void Set_error_on(bool tf)						{this->error_on = tf;}
	bool Get_error_on(void)							{return this->error_on;}
	virtual void warning_msg(const char *err_str);
#else
	virtual bool error_open(const char *file_name, const char * mode = "w");
	void error_flush(void);
	void error_close(void);
	virtual void error_msg(const char * str, bool stop=false);
	void Set_error_file(FILE * out)				{this->error_file = out;}
	FILE *Get_error_file(void)					{return this->error_file;}
	void Set_error_on(bool tf)						{this->error_on = tf;}
	bool Get_error_on(void)							{return this->error_on;}
	virtual void warning_msg(const char *err_str);
#endif

	// dump_ostream
	virtual bool dump_open(const char *file_name, std::ios_base::openmode mode = std::ios_base::out);
	void dump_flush(void);
	void dump_close(void);
	virtual void dump_msg(const char * str);
	void Set_dump_ostream(std::ostream * out)		{this->dump_ostream = out;};
	std::ostream *Get_dump_ostream(void)			{return this->dump_ostream;};
	void Set_dump_on(bool tf)						{this->dump_on = tf;};
	bool Get_dump_on(void)							{return this->dump_on;};

	// fpunchf
	virtual void fpunchf(const char *name, const char *format, double d);
	virtual void fpunchf(const char *name, const char *format, char * d);
	virtual void fpunchf(const char *name, const char *format, int d);
	virtual void fpunchf_end_row(const char *format);
	static void fpunchf_helper(std::ostream *os, const char *format, ...);
	static void fpunchf_helper(std::string *str, const char *format, ...);

	virtual void screen_msg(const char * str);
	void Set_screen_on(bool tf)						{this->screen_on = tf;};
	bool Get_screen_on(void)						{return this->screen_on;};

	// input methods
	virtual int getc(void);
	virtual LINE_TYPE get_line(void);
	virtual LINE_TYPE get_logical_line(void);
	bool check_key(std::string::iterator begin, std::string::iterator end);
	std::string & Get_m_line()       {return m_line;}
	std::string & Get_m_line_save()  {return m_line_save;}
	std::string & Get_accumulated()	 {return accumulated;}
	LINE_TYPE Get_m_line_type()      {return m_line_type;};
	void Set_accumulate(bool tf) 
	{ 
		if (tf)
		{
			accumulated.clear();
		}
		this->accumulate = tf; 
	}
	Keywords::KEYWORDS Get_m_next_keyword() const {return m_next_keyword;}

	// echo 
	enum ECHO_OPTION
	{
		ECHO_LOG,
		ECHO_OUTPUT
	};
	virtual void echo_msg(const char * str);
	void Set_echo_on(bool tf)					{this->echo_on = tf;};
	bool Get_echo_on(void)						{return this->echo_on;};
	void Set_echo_destination(ECHO_OPTION eo)   {this->echo_destination = eo;};
	ECHO_OPTION Get_echo_destination(void)      {return this->echo_destination;};

	// data
protected:
	std::ostream *output_ostream;	
	std::ostream *log_ostream;		
	std::ostream *punch_ostream;	
#ifdef ERROR_OSTREAM
	std::ostream *error_ostream;
#else
	FILE * error_file;
#endif
	std::ostream *dump_ostream;
	int io_error_count;

	bool output_on;
	bool log_on;
	bool punch_on;
	bool error_on;
	bool dump_on;
	bool echo_on;
	bool screen_on;
	ECHO_OPTION echo_destination;

#if defined(_MSC_VER)
/* disable warning C4251: 'identifier' : class 'type' needs to have dll-interface to be used by clients of class 'type2' */
#pragma warning(disable:4251)
#endif

	std::list <std::istream *> istream_list;
	std::list <bool> delete_istream_list;

	std::string m_line;
	std::string m_line_save;
	std::string accumulated;

#if defined(_MSC_VER)
/* reset warning C4251 */
#pragma warning(default:4251)
#endif

	// input data members
	Keywords::KEYWORDS m_next_keyword;
	bool accumulate;
	LINE_TYPE m_line_type;
};

#endif /* _PHRQIO_H */
