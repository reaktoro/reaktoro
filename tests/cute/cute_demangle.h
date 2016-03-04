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
 * Copyright 2009 Peter Sommerlad
 *
 *********************************************************************************/

#ifndef CUTE_DEMANGLE_H_
#define CUTE_DEMANGLE_H_
#include <string>
// needs adaptation for different compilers
// dependency to demangle is a given,
// otherwise we have to use macros everywhere
#ifdef __GNUG__ // also for clang...
#include <cxxabi.h> // __cxa_demangle
#include <cstdlib> // ::free() 
namespace cute {
	extern inline std::string demangle(char const *name);

namespace cute_impl_demangle {
inline std::string plain_demangle(char const *name){
	if (!name) return "unknown";
	char const *toBeFreed = abi::__cxa_demangle(name,0,0,0);
	std::string result(toBeFreed?toBeFreed:name);
	::free(const_cast<char*>(toBeFreed));
	return result;
}
#ifdef _LIBCPP_NAMESPACE
inline void patch_library_namespace(std::string &mightcontaininlinenamespace) {
// libc++ uses inline namespace std::_LIBCPP_NAMESPACE:: for its classes. This breaks the tests relying on meta information. re-normalize the names back to std::
	std::string::size_type pos=std::string::npos;
#define XNS(X) #X
#define NS(X) XNS(X)
#define TOREPLACE "::" NS(_LIBCPP_NAMESPACE)
	std::string const nothing;
	while (std::string::npos != (pos= mightcontaininlinenamespace.find(TOREPLACE)))
			mightcontaininlinenamespace.erase(pos,sizeof(TOREPLACE)-1);
#undef NS
#undef XNS
#undef TOREPLACE
}
inline void patchresultforstring(std::string& result) {
	static const std::string stringid=plain_demangle(typeid(std::string).name());
	std::string::size_type pos=std::string::npos;
	while(std::string::npos != (pos=result.find(stringid))){
		if (!result.compare(pos+stringid.size(),2," >",2)) result.erase(pos+stringid.size(),1); // makes templates look nice
		result.replace(pos,stringid.size(),"std::string");
	}
	patch_library_namespace(result);
}
#endif

}
inline std::string demangle(char const *name){
	if (!name) return "unknown";
	std::string result(cute_impl_demangle::plain_demangle(name));
#ifdef _LIBCPP_NAMESPACE
	cute_impl_demangle::patchresultforstring(result);
#endif
	return result;
}
}
#else
namespace cute {
#ifdef _MSC_VER
namespace cute_demangle_impl {

inline void removeMSKeyword(std::string &name,std::string const &kw){
	std::string::size_type pos=std::string::npos;
	while (std::string::npos != (pos= name.find(kw)))
			name.erase(pos,kw.size());

}
inline void patchresultforstring(std::string& result) {
	static const std::string stringid=(typeid(std::string).name());
	std::string::size_type pos=std::string::npos;
	while(std::string::npos != (pos=result.find(stringid))){
		if (!result.compare(pos+stringid.size(),2," >",2)) result.erase(pos+stringid.size(),1); // makes templates look nice
		result.replace(pos,stringid.size(),"std::string");
	}
}

inline void patchMSMangling(std::string &name){
	patchresultforstring(name);
	removeMSKeyword(name,"class ");
	removeMSKeyword(name,"struct ");
	for (std::string::iterator i=name.begin(); i != name.end(); ++i){
		if (*i==','){  i = name.insert(i+1,' ');}
	}
	std::string::size_type pos=0;
	while(std::string::npos !=(pos=name.find(" ,",pos))){
		name.erase(pos,1);
		++pos;
	}
}
}
inline std::string demangle(char const *name){
	std::string result(name?name:"unknown");
	cute_demangle_impl::patchMSMangling(result);
	return result;
}

#else
// this default works reasonably with MSVC71 and 8, hopefully for others as well
inline std::string demangle(char const *name){
	return std::string(name?name:"unknown");
}
#endif
}
#endif

#endif /* CUTE_DEMANGLE_H_ */
