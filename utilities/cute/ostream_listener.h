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
 * Copyright 2007-2009 Peter Sommerlad
 *
 *********************************************************************************/

#ifndef OSTREAM_LISTENER_H_
#define OSTREAM_LISTENER_H_
#include "cute_listener.h"
#include <iostream>
namespace cute {
	// a "root" listener displaying output, use it as an example on how to build your own, e.g., for XML output
	template <typename Listener=null_listener>
	struct ostream_listener:Listener
	{
		std::ostream &out;
	public:
		ostream_listener(std::ostream &os=std::cerr):out(os) {}
		void begin(suite const &t,char const *info, size_t n_of_tests=t.size()){
			out << "beginning: " << info<<std::endl;
			Listener::begin(t,info, n_of_tests);
		}
		void end(suite const &t, char const *info){
			out << "ending: " << info<<std::endl;
			Listener::end(t,info);
		}
		void start(test const &t){
			out << "starting: " <<t.name()<< std::endl;
			Listener::start(t);
		}
		void success(test const &t, char const *msg){
			out <<  t.name() <<" " << msg<< std::endl;
			Listener::success(t,msg);
		}
		void failure(test const &t,test_failure const &e){
			out << e.filename << ":" << e.lineno << ": testcase failed: " <<e.reason << " in " << t.name()<< std::endl;
			Listener::failure(t,e);
		}
		void error(test const &t, char const *what){
			out << what << " in " << t.name() << std::endl;
			Listener::error(t,what);
		}
	};
}
#endif /*OSTREAM_LISTENER_H_*/
