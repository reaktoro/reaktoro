/*
 * cute_xml_file.h
 *
 *  Created on: 07.06.2013
 *      Author: sop
 */

#ifndef CUTE_XML_FILE_H_
#define CUTE_XML_FILE_H_

#include <fstream>
#include <string>
namespace cute {
struct xml_file_opener {
	std::string filename;
	std::ofstream out;
	xml_file_opener(int argc, char const *const* argv)
	:filename(argc>0&&argv[0]?basename(argv[0]):"testresult.xml")
	,out(filename.c_str()){}
	std::string basename(std::string path){
#if defined( _MSC_VER ) || defined(__MINGW32__)
		char const sep='\\';
#else
		char const sep='/';
#endif
		std::string::size_type pos=path.find_last_of(sep,path.size()-1);
		if (pos != std::string::npos) path.erase(0,pos+1);
		path+=".xml";
		return path;
	}
};
}

#endif /* CUTE_XML_FILE_H_ */
