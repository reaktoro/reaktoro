#if !defined(PARSER_H_INCLUDED)
#define PARSER_H_INCLUDED
#if defined(WIN32)
#include <windows.h>
#endif
#include <string>				// std::string
#include <map>					// std::map
#include <vector>				// std::vector
#include <sstream>				// std::istringstream std::ostringstream
#include <ostream>				// std::ostream
#include <istream>				// std::istream
#include <cctype>               // std::isspace
#include <algorithm>            // std::find_if
#include <functional>           // std::ptr_fun std::not1

#include "PHRQ_base.h"
#include "Keywords.h"
#include "PHRQ_io.h"

#ifdef _DEBUG
#define isspace(a) isspace((a) < -1 ? (a) + 256 : (a)) 
#define isupper(a) isupper((a) < -1 ? (a) + 256 : (a))
#define islower(a) islower((a) < -1 ? (a) + 256 : (a)) 
#define isdigit(a) isdigit((a) < -1 ? (a) + 256 : (a)) 
#define isalpha(a) isalpha((a) < -1 ? (a) + 256 : (a))
#endif

class CParser: public PHRQ_base
{
  public:
	CParser(PHRQ_io *io=NULL);
	CParser(std::istream & input, PHRQ_io *io=NULL);
	//CParser(std::istream & input, std::ostream & output, PHRQ_io *io=NULL);
	//CParser(std::istream & input, std::ostream & output,
	//		std::ostream & error, PHRQ_io *io=NULL);

	virtual ~ CParser();
	enum TOKEN_TYPE
	{
		TT_EMPTY = 2,
		TT_UPPER = 4,
		TT_LOWER = 5,
		TT_DIGIT = 6,
		TT_UNKNOWN = 7
	};

	enum FIND_TYPE
	{
		FT_OK = 0,
		FT_ERROR = 1
	};

	enum OPT_TYPE
	{
		OPT_DEFAULT = -4,
		OPT_ERROR = -3,
		OPT_KEYWORD = -2,
		OPT_EOF = -1
	};
	enum ECHO_OPTION
	{
		EO_NONE = 0,
		EO_ALL = 1,
		EO_KEYWORDS = 2,
		EO_NOKEYWORDS = 3
	};

	enum STATUS_TYPE
	{
		PARSER_ERROR = 0,
		PARSER_OK = 1
	};


		/**
                Function gets a new line and checks for empty, eof, and keywords.

                Arguments:
                        string        Input, character string used in printing error message
                        allow_empty   Input, True or false, if a blank line is accepable
                                      if false, another line is read
                        allow_eof     Input, True or false, if EOF is acceptable
                        allow_keyword Input, True or false, if a keyword is acceptable

                Returns:
                        LT_EMPTY         if empty line read and allow_empty == true
                        LT_KEYWORD       if line begins with keyword
                        LT_EOF           if eof and allow_eof == true
                        LT_OK            otherwise
                        LT_OPTION        if line begins with -[alpha]

           Terminates       if EOF and allow_eof == false.
        */
	PHRQ_io::LINE_TYPE check_line(const std::string & str, bool allow_empty,
						 bool allow_eof, bool allow_keyword, bool print);

		/**
                Read a line from input file put in "line".
                Copy of input line is stored in "line_save".
                Characters after # are discarded in line but retained in "line_save"

                Arguments:
                        None
                Returns:
                        LT_EMPTY,
                        LT_EOF,
                        LT_KEYWORD,
                        LT_OK,
                        LT_OPTION
        */
	PHRQ_io::LINE_TYPE get_line();
	PHRQ_io::LINE_TYPE get_line_phrq_io();

	// bool check_key(const std::string::iterator ptr);
	bool check_key(std::string::iterator begin, std::string::iterator end);

	STATUS_TYPE check_units(std::string & tot_units, bool alkalinity,
							bool check_compatibility,
							const std::string & default_units, bool print);


	//KEY_TYPE next_keyword() const
	Keywords::KEYWORDS next_keyword() const
	{
		return m_next_keyword;
	}
	int get_option(const std::vector < std::string > &opt_list,
				   std::string::iterator & next_char);
	int get_option(const std::vector < std::string > &opt_list,
				   std::istream::pos_type & next_pos);
	int getOptionFromLastLine(const std::vector < std::string > &opt_list,
							  std::string::iterator & next_char, bool flag_error);
	int getOptionFromLastLine(const std::vector < std::string > &opt_list,
							  std::istream::pos_type & next_pos, bool flag_error);


	std::string & line()
	{
		return m_line;
	}
	std::string & line_save()
	{
		return m_line_save;
	}
	std::string & get_accumulated()
	{
		return accumulated;
	}
	void set_accumulate(bool tf) 
	{ 
		if (tf)
		{
			accumulated.clear();
		}
		this->accumulate = tf; 
	}
	std::istringstream & get_iss()
	{
		return m_line_iss;
	}
	int incr_input_error();
	int get_input_error()
	{
		return m_input_error;
	}

		/**
                Copies from begin to token until first space is encountered.

                Arguments:
                        token    output, the token
                        begin    input, begin iterator
                        end      input, end iterator

                Returns:
                        TT_EMPTY
                        TT_UPPER
                        TT_LOWER
                        TT_DIGIT
                        TT_UNKNOWN
        */
	static CParser::TOKEN_TYPE copy_token(std::string & token,
								 std::string::iterator & begin,
								 std::string::iterator & end);
	static CParser::TOKEN_TYPE copy_title(std::string & token,
										std::string::iterator & begin,
										std::string::iterator & end);
	static CParser::TOKEN_TYPE token_type(const std::string & token);
	static CParser::TOKEN_TYPE copy_token(std::string & token, std::istream & is);
	CParser::TOKEN_TYPE copy_token(std::string & token, std::istream::pos_type & pos);
	bool get_true_false(std::istream::pos_type & pos, bool def);
	CParser::TOKEN_TYPE get_rest_of_line(std::string &token);
	static CParser::TOKEN_TYPE parse_delimited(std::string & source, std::string & result, const std::string& t);
	CParser::TOKEN_TYPE peek_token();
	PHRQ_io::LINE_TYPE get_m_line_type(void) const {return this->m_line_type;}

		/**
                Function reads an element name out of the equation string.
                An element name is composed of a capital letter followed by any number
                of lower case characters.

                Arguments:
                        begin    input, points to position in the equation to begin
                                        output, points to next character of equation after
                                        element name.
                        end      input, points to last position in the equation
                        element  input pointer to place to return element character string
        */
	STATUS_TYPE get_elt(std::string::iterator & begin,
						const std::string::iterator end,
						std::string & element);


		/**
                Compares a string value to match beginning letters of a list of options

                Arguments:
                        item    entry: pointer to string to compare
                        n       exit:  item in list that was matched
                        list    entry: pointer to list of character values, assumed to
                                be lower case
                        count_list entry: number of character values in list

                Returns:
                        OK      item matched
                        ERROR   item not matched
                        n       -1      item not matched
                                i       position of match in list
        */
	static FIND_TYPE find_option(const std::string & item, int *n,
								 const std::vector < std::string > &list,
								 bool exact);

	void set_echo_file(ECHO_OPTION opt)
	{
		echo_file = opt;
	}
	ECHO_OPTION get_echo_file()
	{
		return this->echo_file;
	};

	void set_echo_stream(ECHO_OPTION opt)
	{
		echo_stream = opt;
	}
	ECHO_OPTION get_echo_stream()
	{
		return this->echo_stream;
	};

	STATUS_TYPE parse_couple(std::string & token);

	template <class T>
	STATUS_TYPE addPair(std::map < std::string, T >&totals,
						std::istream::pos_type & pos);

  protected:
	PHRQ_io::LINE_TYPE get_logical_line();

  protected:
	std::istream & m_input_stream;
	//std::ostream & m_output_stream;
	//std::ostream & m_error_stream;
	int m_input_error;
	//KEY_TYPE m_next_keyword;
	Keywords::KEYWORDS m_next_keyword;
	std::string m_line;
	std::string m_line_save;
	std::istringstream m_line_iss;
	PHRQ_io::LINE_TYPE m_line_type;
	ECHO_OPTION echo_stream;
	ECHO_OPTION echo_file;
	std::string accumulated;
	bool accumulate;
	bool phrq_io_only;

};

// Global functions
static inline std::string &trim_left(std::string &s)
{ 
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace)))); 
	return s; 
} 
static inline std::string &trim_right(std::string &s)
{ 
	s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end()); 
	return s; 
} 
static inline std::string &trim(std::string &s)
{ 
	return trim_left(trim_right(s)); 
} 

#endif // PARSER_H_INCLUDED
