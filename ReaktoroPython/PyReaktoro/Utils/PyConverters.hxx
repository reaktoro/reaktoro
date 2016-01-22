// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

// Boost includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
namespace py = boost::python;

// C++ includes
#include <iostream>
#include <vector>

namespace std {

template<typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& vec)
{
    out << "[";

    for(size_t i = 0; i < vec.size(); ++i)
        out << vec[i] << ((i+1 != vec.size()) ? ", " : "");

    out << "]";

    return out;
}

template<typename T, typename X>
std::ostream& operator<<(std::ostream& out, const std::map<T, X>& map)
{
    out << "{";

    unsigned size = map.size();
    unsigned i = 0;
    for(const auto& pair : map)
    {
        out << pair.first << ": " << pair.second << ((i+1 != size) ? ", " : "");
        ++i;
    }

    out << "}";

    return out;
}

} // namespace std

namespace Reaktoro {

template<typename T>
struct std_vector_to_python_list
{
    static PyObject* convert(std::vector<T> const& v)
    {
        py::list l;
        for(const T& val : v)
            l.append(val);
        return py::incref(l.ptr());
    }
};

template<typename T>
struct std_vector_from_python_list
{
    std_vector_from_python_list()
    {
        py::converter::registry::push_back(&convertible, &construct, py::type_id<std::vector<T>>());
    }

    static void* convertible(PyObject* obj_ptr)
    {
        if(!PyList_Check(obj_ptr))
            return 0;
        return obj_ptr;
    }

    static void construct(PyObject* obj_ptr, py::converter::rvalue_from_python_stage1_data* data)
    {
        py::extract<py::list> x(obj_ptr);
        if(!x.check())
            py::throw_error_already_set();

        py::list l = x();

        void *storage =
            ((py::converter::rvalue_from_python_storage<std::vector<T>>*) data)->
                storage.bytes;

        new (storage) std::vector<T>();

        std::vector<T>& v = *reinterpret_cast<std::vector<T>*>(storage);

        for(int idx = 0; idx < len(l); ++idx)
        {
            py::extract<T> ext(l[idx]);
            if(!ext.check())
            {
                v.~vector();
                py::throw_error_already_set();
            }

            v.push_back(ext());
        }

        data->convertible = storage;
    }
};

template<typename T>
void init_converter_std_vector_from_python_list()
{
    std_vector_from_python_list<T>();
}

template<typename T>
void export_std_vector(const char* type)
{
    init_converter_std_vector_from_python_list<T>();

    py::class_<std::vector<T>>(type)
        .def(py::vector_indexing_suite<std::vector<T>>());
}

template<typename T>
void export_std_vector_with_str(const char* type)
{
    init_converter_std_vector_from_python_list<T>();

    py::class_<std::vector<T>>(type)
        .def(py::vector_indexing_suite<std::vector<T>>())
        .def(py::self_ns::str(py::self_ns::self));
}

template<typename Key, typename Value>
void export_std_map(const char* type)
{
    py::class_<std::map<Key, Value>>(type)
        .def(py::map_indexing_suite<std::map<Key, Value>>());
}

template<typename Key, typename Value>
void export_std_map_with_str(const char* type)
{
    py::class_<std::map<Key, Value>>(type)
        .def(py::map_indexing_suite<std::map<Key, Value>>())
        .def(py::self_ns::str(py::self_ns::self));
}

} // namespace Reaktoro
