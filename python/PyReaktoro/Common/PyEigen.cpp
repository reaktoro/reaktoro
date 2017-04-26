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

#include "PyEigen.hpp"

// Numpy includes
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/ndarrayobject.h>

// Boost includes
#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <boost/python/slice.hpp>
#include <boost/smart_ptr.hpp>
namespace py = boost::python;

// Eigen includes
#include <Reaktoro/Math/Eigen/Core>

namespace Reaktoro {

template<typename T> auto numtype(T) -> NPY_TYPES { return NPY_DOUBLE; }
template<> auto numtype(int) -> NPY_TYPES { return NPY_INT; }
template<> auto numtype(long) -> NPY_TYPES { return NPY_LONG; }
template<> auto numtype(float) -> NPY_TYPES { return NPY_FLOAT; }
template<> auto numtype(double) -> NPY_TYPES { return NPY_DOUBLE; }

void IndexError(std::string msg)
{
    PyErr_SetString(PyExc_IndexError, msg.c_str());
    py::throw_error_already_set();
}

void SliceStepError()
{
    PyErr_SetString(PyExc_IndexError, "slice step size not supported");
    py::throw_error_already_set();
}

void TypeError(std::string msg)
{
    PyErr_SetString(PyExc_TypeError, msg.c_str());
    py::throw_error_already_set();
}

template <typename Scalar>
struct WrapperEigenVector
{
    typedef Eigen::Matrix<Scalar, -1, 1> VectorType;

    static auto init_default() -> boost::shared_ptr<VectorType>
    {
        return boost::shared_ptr<VectorType>(new VectorType());
    }

    static auto init_copy(const VectorType& copy) -> boost::shared_ptr<VectorType>
    {
        return boost::shared_ptr<VectorType>(new VectorType(copy));
    }

    static auto init_with_rows(int rows) -> boost::shared_ptr<VectorType>
    {
        return boost::shared_ptr<VectorType>(new VectorType(rows));
    }

    static auto init_with_rows_and_value(int rows, const Scalar& val) -> boost::shared_ptr<VectorType>
    {
        boost::shared_ptr<VectorType> vec = init_with_rows(rows);
        std::fill_n(vec->data(), rows, val);
        return vec;
    }

    static auto init_with_array(const py::numeric::array& array) -> boost::shared_ptr<VectorType>
    {
        long rows = len(array);
        boost::shared_ptr<VectorType> vec = init_with_rows(rows);
        for(long i = 0; i < rows; ++i)
        {
            py::extract<Scalar> ext(array[i]);
            if(!ext.check())
            {
                vec.reset();
                py::throw_error_already_set();
            }
            (*vec)[i] = ext();
        }
        return vec;
    }

    template <typename Sequence>
    static auto init_with_sequence(const Sequence& seq) -> boost::shared_ptr<VectorType>
    {
        return init_with_array(py::numeric::array(seq));
    }

    static auto size(VectorType& self) -> int
    {
        return self.size();
    }

    static auto get_slice_data(VectorType& self, const py::slice& slice, long& start, long& stop) -> void
    {
        start = 0;
        stop  = self.rows();
        if(!slice.start().is_none())
            start = py::extract<long>(slice.start());
        if(!slice.stop().is_none())
            stop = py::extract<long>(slice.stop());
        if(start < 0 || stop < start || stop > self.rows())
            IndexError("slice index out of range");
    }

    static auto set_item(VectorType& self, int i, const Scalar& scalar) -> void
    {
        self[i] = scalar;
    }

    static auto get_item(VectorType& self, int i) -> Scalar&
    {
        return self[i];
    }

    static auto set_slice_with_scalar(VectorType& self, const py::slice& slice, const Scalar& scalar) -> void
    {
        long start, stop;
        get_slice_data(self, slice, start, stop);
        for(long i = start; i < stop; ++i)
            self[i] = scalar;
    }

    static auto set_slice_with_eigen_vector(VectorType& self, const py::slice& slice, const VectorType& vec) -> void
    {
        long start, stop;
        get_slice_data(self, slice, start, stop);
        const long rows = stop - start;
        if(vec.rows() != rows)
            IndexError("mismatch number of rows");
        self.segment(start, rows) = vec;
    }

    static auto set_slice_with_array(VectorType& self, const py::slice& slice, const py::numeric::array& array) -> void
    {
        long start, stop;
        get_slice_data(self, slice, start, stop);
        for(long i = start; i < stop; ++i)
        {
            py::extract<Scalar> ext(array[i]);
            if(!ext.check())
                py::throw_error_already_set();
            self[i] = ext();
        }
    }

    template <typename Sequence>
    static auto set_slice_with_sequence(VectorType& self, const py::slice& slice, const Sequence& seq) -> void
    {
        set_slice_with_array(self, slice, py::numeric::array(seq));
    }

    static auto get_slice(VectorType& self, const py::slice& slice) -> py::object
    {
        long start, stop;
        get_slice_data(self, slice, start, stop);
        const long rows = stop - start;
        return py::object(VectorType(self.segment(start, rows)));
    }

    static auto array(VectorType& self) -> py::numeric::array
    {
        const int rows = self.rows();
        npy_intp dims[1] = {rows};
        py::object obj(py::handle<>(py::incref(PyArray_SimpleNew(1, dims, numtype(Scalar())))));
        Scalar* data = (Scalar*)PyArray_DATA((PyArrayObject*)obj.ptr());
        int k = 0; for(long i = 0; i < rows; ++i)
            data[k++] = self[i];
        return py::extract<py::numeric::array>(obj);
    }

    static auto data(VectorType& self) -> py::numeric::array
    {
        const int size = self.size();
        npy_intp dims[1] = {size};
        py::object obj(py::handle<>(py::incref(
            PyArray_SimpleNewFromData(1, dims, numtype(Scalar()), self.data()))));
        return py::extract<py::numeric::array>(obj);
    }

    static auto begin(VectorType& self) -> Scalar*
    {
        return self.data();
    }

    static auto end(VectorType& self) -> Scalar*
      {
        return self.data() + self.rows() * self.cols();
    }
};

template<typename Scalar>
struct WrapperEigenMatrix
{
    using MatrixType = Eigen::Matrix<Scalar, -1, -1>;

    static auto init_default() -> boost::shared_ptr<MatrixType>
    {
        return boost::shared_ptr<MatrixType>(new MatrixType());
    }

    static auto init_copy(const MatrixType& copy) -> boost::shared_ptr<MatrixType>
    {
        return boost::shared_ptr<MatrixType>(new MatrixType(copy));
    }

    static auto init_with_rows_cols(int rows, int cols) -> boost::shared_ptr<MatrixType>
    {
        return boost::shared_ptr<MatrixType>(new MatrixType(rows, cols));
    }

    static auto init_with_rows_cols_and_value(int rows, int cols, const Scalar& val) -> boost::shared_ptr<MatrixType>
    {
        boost::shared_ptr<MatrixType> mat = init_with_rows_cols(rows, cols);
        mat->fill(val);
        return mat;
    }

    static auto init_with_array(const py::numeric::array& array) -> boost::shared_ptr<MatrixType>
    {
        py::tuple shape = py::extract<py::tuple>(array.attr("shape"));
        long rows, cols;
        if(len(shape) == 1)
        {
            rows = 1;
            cols = py::extract<long>(shape[0]);
        }
        else
        {
            rows = py::extract<long>(shape[0]);
            cols = py::extract<long>(shape[1]);
        }
        boost::shared_ptr<MatrixType> mat = init_with_rows_cols(rows, cols);
        for(long i = 0; i < rows; ++i) for(long j = 0; j < cols; ++j)
            (*mat)(i, j) = py::extract<Scalar>(array(i, j));
        return mat;
    }

    template<typename Sequence>
    static auto init_with_sequence(const Sequence& seq) -> boost::shared_ptr<MatrixType>
    {
        return init_with_array(py::numeric::array(seq));
    }

    static auto size(MatrixType& self) -> int
    {
        return self.size();
    }

    static auto get_slice_data(MatrixType& self, const py::object& arg, long& start, long& stop, long max) -> void
    {
        py::extract<py::slice> xslice(arg);
        if(!xslice.check())
        {
            start = py::extract<long>(arg);
            stop = start + 1;
        }
        else
        {
            py::slice slice = xslice();
            if(!slice.step().is_none())
                SliceStepError();
            start = 0;
            stop = max;
            if(!slice.start().is_none())
                start = py::extract<long>(slice.start());
            if(!slice.stop().is_none())
                stop = py::extract<long>(slice.stop());
        }
        if(start < 0 || stop < start || stop > max)
            IndexError("slice index out of range");
    }

    static auto set_item(MatrixType& self, const py::tuple& tuple, const Scalar& scalar) -> void
    {
        long row_start, row_stop;
        long col_start, col_stop;
        get_slice_data(self, tuple[0], row_start, row_stop, self.rows());
        get_slice_data(self, tuple[1], col_start, col_stop, self.cols());
        if(row_stop - row_start == 1 && col_stop - row_start == 1)
            self(row_start, col_start) = scalar;
        else
            for(long i = row_start; i < row_stop; ++i)
                for(long j = col_start; j < col_stop; ++j)
                    self(i, j) = scalar;
    }

    static auto get_item(MatrixType& self, const py::tuple& tuple) -> py::object
    {
        long row_start, row_stop;
        long col_start, col_stop;
        get_slice_data(self, tuple[0], row_start, row_stop, self.rows());
        get_slice_data(self, tuple[1], col_start, col_stop, self.cols());
        const long row_size = row_stop - row_start;
        const long col_size = col_stop - col_start;
        if(row_size == 1 && col_size == 1)
            return py::object(self(row_start, col_start));
        else
            return py::object(MatrixType(self.block(row_start, col_start, row_size, col_size)));
    }

    static auto set_slice_with_eigen_matrix(MatrixType& self, const py::tuple& tuple, const MatrixType& mat) -> void
    {
        long row_start, row_stop;
        long col_start, col_stop;
        get_slice_data(self, tuple[0], row_start, row_stop, self.rows());
        get_slice_data(self, tuple[1], col_start, col_stop, self.cols());
        self.block(row_start, col_start, mat.rows(), mat.cols()) = mat;
    }

    static auto set_slice_with_array(MatrixType& self, const py::tuple& tuple, const py::numeric::array& array) -> void
    {
        long row_start, row_stop;
        long col_start, col_stop;
        get_slice_data(self, tuple[0], row_start, row_stop, self.rows());
        get_slice_data(self, tuple[1], col_start, col_stop, self.cols());
        long rows = row_stop - row_start;
        long cols = col_stop - col_start;
        for(long i = row_start; i < rows; ++i) for(long j = col_start; j < cols; ++j)
            self(i, j) = py::extract<Scalar>(array(i, j));
    }

    template<typename Sequence>
    static auto set_slice_with_sequence(MatrixType& self, const py::tuple& tuple, const Sequence& seq) -> void
    {
        set_slice_with_array(self, tuple, py::numeric::array(seq));
    }

    static auto array(MatrixType& self) -> py::numeric::array
    {
        const int rows = self.rows();
        const int cols = self.cols();
        npy_intp dims[2] = {rows, cols};
        py::object obj(py::handle<>(py::incref(PyArray_SimpleNew(2, dims, numtype(Scalar())))));
        Scalar* data = (Scalar*)PyArray_DATA((PyArrayObject*)obj.ptr());
        int k = 0; for(long i = 0; i < rows; ++i) for(long j = 0; j < cols; ++j)
            data[k++] = self(i, j);
        return py::extract<py::numeric::array>(obj);
    }

    static auto data(MatrixType& self) -> py::numeric::array
    {
        const int size = self.size();
        npy_intp dims[1] = {size};
        py::object obj(py::handle<>(py::incref(
            PyArray_SimpleNewFromData(1, dims, numtype(Scalar()), self.data()))));
        return py::extract<py::numeric::array>(obj);
    }

    static auto begin(MatrixType& self) -> Scalar*
    {
        return self.data();
    }

    static auto end(MatrixType& self) -> Scalar*
    {
        return self.data() + self.rows() * self.cols();
    }
};


template<typename Scalar>
struct EigenVector_from_python_array
{
    using VectorType = Eigen::Matrix<Scalar, -1, 1>;

    EigenVector_from_python_array()
    {
        py::converter::registry::push_back(&convertible, &construct, py::type_id<VectorType>());
    }

    static void* convertible(PyObject* obj_ptr)
    {
        if(!PyArray_Check(obj_ptr))
            return 0;
        return obj_ptr;
    }

    static void construct(PyObject* obj_ptr, py::converter::rvalue_from_python_stage1_data* data)
    {
        py::extract<py::numeric::array> x(obj_ptr);
        if(!x.check())
            py::throw_error_already_set();

        py::numeric::array l = x();

        void *storage =
            ((py::converter::rvalue_from_python_storage<VectorType>*) data)->
                storage.bytes;

        new (storage) VectorType();

        VectorType& v = *reinterpret_cast<VectorType*>(storage);

        v.resize(len(l));
        for(int idx = 0; idx < v.rows(); ++idx)
        {
            py::extract<Scalar> ext(l[idx]);
            if(!ext.check())
            {
                v.~VectorType();
                py::throw_error_already_set();
            }

            v[idx] = ext();
        }

        data->convertible = storage;
    }
};

template<typename Scalar>
void init_converter_EigenVector_from_python_array()
{
    EigenVector_from_python_array<Scalar>();
}

template <typename Scalar>
void export_EigenVectorType(const char* class_name)
{
    using VectorType = Eigen::Matrix<Scalar, -1, 1>;
    using WrapperType = WrapperEigenVector<Scalar>;
    using Index = typename VectorType::Index;

    auto res =
    py::class_<VectorType, boost::shared_ptr<VectorType>>(class_name)
        .def("__init__", py::make_constructor(&WrapperType::init_default))
        .def("__init__", py::make_constructor(&WrapperType::init_copy))
        .def("__init__", py::make_constructor(&WrapperType::init_with_rows))
        .def("__init__", py::make_constructor(&WrapperType::init_with_rows_and_value))
        .def("__init__", py::make_constructor(&WrapperType::init_with_array))
        .def("__init__", py::make_constructor(&WrapperType::template init_with_sequence<py::list>))
        .def("__init__", py::make_constructor(&WrapperType::template init_with_sequence<py::tuple>))
        .def("__len__", WrapperType::size)
        .def("__setitem__", WrapperType::set_item, py::with_custodian_and_ward<1,2>()) // to let container keep value
        .def("__getitem__", WrapperType::get_item, py::return_value_policy<py::copy_non_const_reference>())
        .def("__setitem__", WrapperType::set_slice_with_scalar)
        .def("__setitem__", WrapperType::set_slice_with_eigen_vector)
        .def("__setitem__", WrapperType::set_slice_with_array)
        .def("__setitem__", WrapperType::template set_slice_with_sequence<py::list>)
        .def("__setitem__", WrapperType::template set_slice_with_sequence<py::tuple>)
        .def("__getitem__", WrapperType::get_slice)
        .def("__iter__", py::range(&WrapperType::begin, &WrapperType::end))
        .def("__array__", WrapperType::data)
        .def("array", WrapperType::array)
        .def("data", WrapperType::data)
        .def("size", WrapperType::size)
        .def("rows", &VectorType::rows)
        .def("cols", &VectorType::cols)
        .def("resize", static_cast<void(VectorType::*)(Index)>(&VectorType::resize))
        .def("fill", &VectorType::fill)
        .def(py::self_ns::str(py::self_ns::self))
        ;

    init_converter_EigenVector_from_python_array<Scalar>();
}

template<typename Scalar>
auto export_EigenMatrixType(const char* class_name) -> void
{
    using MatrixType = Eigen::Matrix<Scalar, -1, -1>;
    using WrapperType = WrapperEigenMatrix<Scalar>;
    using Index = typename MatrixType::Index;

    auto res =
    py::class_<MatrixType, boost::shared_ptr<MatrixType>>(class_name)
        .def("__init__", py::make_constructor(&WrapperType::init_default))
        .def("__init__", py::make_constructor(&WrapperType::init_copy))
        .def("__init__", py::make_constructor(&WrapperType::init_with_rows_cols))
        .def("__init__", py::make_constructor(&WrapperType::init_with_rows_cols_and_value))
        .def("__init__", py::make_constructor(&WrapperType::init_with_array))
        .def("__init__", py::make_constructor(&WrapperType::template init_with_sequence<py::list>))
        .def("__init__", py::make_constructor(&WrapperType::template init_with_sequence<py::tuple>))
        .def("__len__", &WrapperType::size)
        .def("__setitem__", &WrapperType::set_item, py::with_custodian_and_ward<1, 2>()) // to let container keep value
        .def("__getitem__", &WrapperType::get_item)
        .def("__setitem__", &WrapperType::set_slice_with_eigen_matrix)
        .def("__setitem__", &WrapperType::set_slice_with_array)
        .def("__setitem__", &WrapperType::template set_slice_with_sequence<py::list>)
        .def("__setitem__", &WrapperType::template set_slice_with_sequence<py::tuple>)
        .def("__iter__", py::range(&WrapperType::begin, &WrapperType::end))
        .def("__array__", &WrapperType::data)
        .def("array", &WrapperType::array)
        .def("data", &WrapperType::data)
        .def("size", &WrapperType::size)
        .def("rows", &MatrixType::rows)
        .def("cols", &MatrixType::cols)
        .def("resize", static_cast<void(MatrixType::*)(Index,Index)>(&MatrixType::resize))
        .def("fill", &MatrixType::fill)
        .def(py::self_ns::str(py::self_ns::self))
        ;
}

void export_EigenVector()
{
    export_EigenVectorType<double>("VectorXd");
    export_EigenVectorType<float>("VectorXf");
    export_EigenVectorType<int>("VectorXi");
}

void export_EigenMatrix()
{
    export_EigenMatrixType<double>("MatrixXd");
    export_EigenMatrixType<float>("MatrixXf");
    export_EigenMatrixType<int>("MatrixXi");
}

auto export_Eigen() -> void
{
    import_array();
    export_EigenVector();
    export_EigenMatrix();
}

} // namespace Reaktoro
