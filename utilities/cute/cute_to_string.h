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

#ifndef CUTE_TO_STRING_H_
#define CUTE_TO_STRING_H_
#include <string>
#include <algorithm>
namespace cute {
namespace cute_to_string {
		static inline std::string backslashQuoteTabNewline(std::string const &input){
			std::string result;
			result.reserve(input.size());
			for (std::string::size_type i=0; i < input.length() ; ++i){
				switch(input[i]) {
					case '\n': result += "\\n"; break;
					case '\t': result += "\\t"; break;
					case '\\': result += "\\\\"; break;
					case '\r': result += "\\r"; break;
					default: result += input[i];
				}
			}
			return result;
		}
		// common overloads of interface that work without an ostream
		static inline std::string to_string(char const *const &s){
			return s;
		}
		static inline std::string to_string(std::string const &s){
			return s;
		}
	}
}
#ifndef DONT_USE_IOSTREAM
#include <ostream>
#include <sstream>
#include <typeinfo>
#ifdef _MSC_VER
#include <map>
#include <set>
#endif
#include "cute_demangle.h"
namespace cute {
namespace cute_to_string {
		template <typename T>
		std::ostream &to_stream(std::ostream &os,T const &t); // recursion needs forward

		// the following code was stolen and adapted from Boost Exception library.
		// it avoids compile errors, if a type used with ASSERT_EQUALS doesn't provide an output shift operator
		namespace to_string_detail {
			template <class T,class CharT,class Traits>
			char operator<<( std::basic_ostream<CharT,Traits> &, T const & );
			template <class T,class CharT,class Traits>
			struct is_output_streamable_impl {
				static std::basic_ostream<CharT,Traits> & f();
				static T const & g();
				enum e { value = (sizeof(char) != sizeof(f()<<g())) }; // assumes sizeof(char)!=sizeof(ostream&)
			};
			// specialization for pointer types to map char * to operator<<(std::ostream&,char const *)
			template <class T,class CharT,class Traits>
			struct is_output_streamable_impl<T*,CharT,Traits> {
				static std::basic_ostream<CharT,Traits> & f();
				static T const * g();
				enum e { value = (sizeof(char) != sizeof(f()<<g())) }; // assumes sizeof(char)!=sizeof(ostream&)
			};
			template <class CONT>
			struct has_begin_end_const_member {
				template <typename T, T, T> struct type_check;
				template <typename C> static typename C::const_iterator test(
						type_check<typename C::const_iterator (C::*)()const,&C::begin, &C::end>*);
				template <typename C> static char test(...);
				enum e { value = (sizeof(char) != sizeof(test<CONT>(0)))
				};
			}; // this doesn't work with VS library for set/map due to implementation in parent
			   // no useful workaround possible
			   //-> internal compiler errors for VS2010/12/12CTP Nov 12
		}
		template <class T, class CharT=char, class Traits=std::char_traits<CharT> >
		struct is_output_streamable {
			enum e { value=to_string_detail::is_output_streamable_impl<T,CharT,Traits>::value };
		};
		// detect standard container conforming begin() end() iterator accessors.
		// might employ begin/end traits from c++0x for loop in the future. --> select_container
		template <typename T>
		struct printItWithDelimiter
		{
			std::ostream &os;
			bool first; // allow use of for_each algorithm
			printItWithDelimiter(std::ostream &os):os(os),first(true){}
			void operator()(T const &t){
				if (!first) os<<',';
				else first=false;
				os << '\n'; // use newlines so that CUTE's plug-in result viewer gives nice diffs
				cute_to_string::to_stream<T>(os,t);
			}
		};
		//try print_pair with specialization of template function instead:
		// the generic version prints about missing operator<< that is the last resort
		template <typename T>
		std::ostream &print_pair(std::ostream &os,T const &t){
			return os << "no operator<<(ostream&, " <<cute::demangle(typeid(T).name())<<')';
		}
		//the std::pair overload is useful for std::map etc. however,
		template <typename K, typename V>
		std::ostream &print_pair(std::ostream &os,std::pair<K,V> const &p){
			os << '[' ;
			cute_to_string::to_stream(os,p.first);
			os << " -> ";
			cute_to_string::to_stream(os,p.second);
			os << ']';
			return os;
		}
		// overload for Arrays
		template <typename T, size_t N>
		std::ostream &print_pair(std::ostream &os,T const (&t)[N]){
			printItWithDelimiter<T> printer(os);
			os << cute::demangle(typeid(T).name()) <<'['<<N<<']'<< '{';
			std::for_each(t,t+N,printer);
			return os << '}';
		}

		template <typename T, bool select>
		struct select_container {
			std::ostream &os;
			select_container(std::ostream &os):os(os){}
			std::ostream& operator()(T const &t){
				printItWithDelimiter<typename T::value_type> printer(os);
				os << cute::demangle(typeid(T).name()) << '{';
				std::for_each(t.begin(),t.end(),printer);
				return os << '}';
			}
		};

		template <typename T>
		struct select_container<T,false> {
			std::ostream &os;
			select_container(std::ostream &os):os(os){}
			std::ostream & operator()(T const &t){
				//  look for std::pair. a future with tuple might be useful as well, but not now.
				return print_pair(os,t); // here a simple template function overload works.
			}
		};

		template <typename T, bool select>
		struct select_built_in_shift_if {
			std::ostream &os;
			select_built_in_shift_if(std::ostream &ros):os(ros){}
			std::ostream& operator()(T const &t){
				return os << t ; // default uses operator<<(std::ostream&,T const&) if available
			}
		};

		template <typename T>
		struct select_built_in_shift_if<T,false> {
			std::ostream &os;
			select_built_in_shift_if(std::ostream &ros):os(ros){}
			std::ostream & operator()(T const &t){
				// if no operator<< is found, try if it is a container or std::pair
				return select_container<T,bool(to_string_detail::has_begin_end_const_member<T>::value) >(os)(t);
			}
		};
		template <typename T>
		std::ostream &to_stream(std::ostream &os,T const &t){
			select_built_in_shift_if<T,cute_to_string::is_output_streamable<T>::value > out(os);
			return out(t);
		}
#ifdef _MSC_VER
		// special overloads because VC can not detect begin/end
		inline std::ostream& to_stream(std::ostream &os,std::string const &s){
			return os<<s;
		} // needed to compensate for following overload, hope nothing else matches
		template <template<typename,typename,typename> class S,
		          typename K, typename CMP, typename ALLOC>
		std::ostream &to_stream(std::ostream &os,S<K,CMP,ALLOC> const &t){
			printItWithDelimiter<typename S<K,CMP,ALLOC>::value_type> printer(os);
			os << cute::demangle(typeid(S<K,CMP,ALLOC>).name()) << '{';
			std::for_each(t.begin(),t.end(),printer);
			return os << '}';
		}
		template <template<typename,typename,typename,typename> class M,typename K, typename V, typename CMP, typename ALLOC>
		std::ostream &to_stream(std::ostream &os,M<K,V,CMP,ALLOC> const &t){
			printItWithDelimiter<typename M<K,V,CMP,ALLOC>::value_type> printer(os);
			os << cute::demangle(typeid(M<K,V,CMP,ALLOC>).name()) << '{';
			std::for_each(t.begin(),t.end(),printer);
			return os << '}';
		}
#endif

		// this is the interface:
		template <typename T>
		std::string to_string(T const &t) {
			std::ostringstream os;
			to_stream(os,t);
			return os.str();
		}
	}
}
#else
#include "cute_determine_traits.h"
#include <limits>
// traits
namespace cute{
namespace cute_to_string {
		template <typename T>
		void adjust_long(T const &,std::string &to_adjust){ // assumes T is an integral type
			if (sizeof(T) <= sizeof(int)) return; // don't mark int sized integrals with L
			if (sizeof(T)>=sizeof(long)) to_adjust+='L';
			if (sizeof(T)> sizeof(long)) to_adjust+='L'; // if there is support for (unsigned) long long
		}
		template <typename T>
		std::string to_string_embedded_int_signed(T const &t, impl_place_for_traits::true_type ){
			std::string convert; // t is an integral value
			T x=t;
			bool negative=t<0;
			bool minint=false;
			if (x == std::numeric_limits<T>::min()){ // can not easily convert it, assuming 2s complement
				minint=true;
				x +=1;
			}
			if (x < 0) x = -x;
			if (x == 0) convert += '0';
			while (x > 0) {
				convert += "0123456789"[x%10];
				x /= 10;
			}
			if (minint) ++ convert[0]; // adjust last digit
			if (negative) convert += '-';
			reverse(convert.begin(),convert.end());
			cute::cute_to_string::adjust_long(t,convert);
			return convert;
		}
		template <typename T>
		std::string hexit(T const &t){ // must be an unsigned type
			std::string hexed;
			if (t == 0) hexed+='0';
			for (T x=t;x>0;x /= 16){
				hexed += "0123456789ABCDEF"[x%16];
			}
			reverse(hexed.begin(),hexed.end());
			return hexed;
		}
		template <typename T>
		std::string to_string_embedded_int_signed(T const &t, impl_place_for_traits::false_type ){
			// manual hex conversion to avoid ostream dependency for unsigned values
			std::string hexed="0x"+cute::cute_to_string::hexit(t);
			cute::cute_to_string::adjust_long(t,hexed);
			return hexed;
		}
		template <typename T>
		std::string to_string_embedded_int(T const &t, impl_place_for_traits::true_type ){
			return to_string_embedded_int_signed(t,impl_place_for_traits::is_signed<T>());
		}
		template <typename T>
		std::string to_string_embedded_int(T const &t, impl_place_for_traits::false_type ){
			return "no to_string";
		}
		// convenience for pointers.... useful?
		template <typename T>
		std::string to_string(T * const&t) {
			std::string result;
			if (sizeof(T *) <= sizeof(unsigned long))
				result = cute::cute_to_string::hexit(reinterpret_cast<unsigned long>(t));
			else
#if defined(USE_STD11) /* should allow for all compilers supporting ULL*/
			result = "p"+cute::cute_to_string::hexit(reinterpret_cast<unsigned long long>(t));
#else
			return "no to_string";
#endif
			result.insert(0u,sizeof(T*)*2-result.size(),'0');
			result.insert(0,1,'p');
			return result;
		}

		// this is the interface:
		template <typename T>
		std::string to_string(T const &t) {
			return to_string_embedded_int(t, impl_place_for_traits::is_integral<T>());
		}
	}
}
#endif
#endif /* CUTE_TO_STRING_H_ */
