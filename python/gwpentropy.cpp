#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "../gwpentropy.h"

namespace py = pybind11;

std::tuple<double, double> calc(py::array_t<double> buf, int n, int w, double q) {
    double H, C;
    double* buf_ptr = static_cast<double*>(buf.request().ptr);
    gwpentropy(&H, &C, buf_ptr, n, w, q);
    return std::tuple<double, double>(H, C);
}

PYBIND11_MODULE(TORCH_EXTENSION_NAME, m)
{
    m.def("calc", &calc, "calc");
}