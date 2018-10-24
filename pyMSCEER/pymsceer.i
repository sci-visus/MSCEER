/*pymsceer.i*/
%module pymsceer

%include "std_vector.i"
namespace std {
    %template(VectorIndex) std::vector<long long>;
    %template(VectorF) std::vector<float>;
};

%{
#include "pymsceer.h"
#include "py_basic_types.h"
#include "py_grad_builder.h"
#include "py_mesh.h"
#include "py_msc.h"
%}

%include "pymsceer.h"
%include "py_basic_types.h"
%include "py_grad_builder.h"
%include "py_mesh.h"
%include "py_msc.h"
