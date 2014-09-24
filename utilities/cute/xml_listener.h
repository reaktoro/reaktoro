/*********************************************************************************
 * This file is part of CUTE.
 *
 * CUTE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CUTE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with CUTE.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Copyright 2007-2013 Peter Sommerlad
 *
 *********************************************************************************/

#ifndef XML_LISTENER_H_
#define XML_LISTENER_H_
#include "cute_listener.h"
#include "cute_xml_file.h" // for convenience
#include <ostream>
namespace cute {
	template <typename Listener=null_listener>
	class xml_listener:public Listener
	{
	protected:
		std::string mask_xml_chars(std::string in){
			std::string::size_type pos=0;
			while((pos=in.find_first_of("<&\"",pos))!=std::string::npos){
				switch(in[pos]){
				case '&': in.replace(pos,1,"&amp;"); pos +=5; break;
				case '<': in.replace(pos,1,"&lt;"); pos += 4; break;
				case '"': in.replace(pos,1,"&quot;"); pos+=6; break;
				default: throw "oops";break;
				}
			}
			return in;
		}
		std::ostream &out;
		std::string current_suite;
	public:
		xml_listener(std::ostream &os):out(os) {
			out << "<testsuites>\n";
		}
		~xml_listener(){
			out << "</testsuites>\n"<< std::flush;
		}

		void begin(suite const &t,char const *info, size_t n_of_tests){
			current_suite=mask_xml_chars(info);
			out << "\t<testsuite name=\"" << current_suite << "\" tests=\"" << n_of_tests << "\">\n";
			Listener::begin(t,info, n_of_tests);
		}
		void end(suite const &t, char const *info){
			out << "\t</testsuite>\n";
			current_suite.clear();
			Listener::end(t,info);
		}
		void start(test const &t){
			out << "\t\t<testcase classname=\""<<current_suite <<"\" name=\""<< mask_xml_chars(t.name())<<"\"";
			Listener::start(t);
		}
		void success(test const &t, char const *msg){
			out << "/>\n";
			Listener::success(t,msg);
		}
		void failure(test const &t,test_failure const &e){
			out <<  ">\n\t\t\t<failure message=\"" << mask_xml_chars(e.filename) << ":" << e.lineno << " "
				<< mask_xml_chars(e.reason) << "\">\n"<<mask_xml_chars(e.reason)<<"\n\t\t\t</failure>\n\t\t</testcase>\n";
			Listener::failure(t,e);
		}
		void error(test const &t, char const *what){
			out << ">\n\t\t\t<error message=\"" << mask_xml_chars(t.name()) << " " << mask_xml_chars(what)
				<< "\" type=\"unexpected exception\">\n"<<mask_xml_chars(what)
				<<"\n\t\t\t</error>\n\t\t</testcase>\n";
			Listener::error(t,what);
		}
	};
}

#endif /*IDE_LISTENER_H_*/
