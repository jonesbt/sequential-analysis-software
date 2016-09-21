%module sequential_analysis_wrap
%{
#include "sequential_analysis.hpp"
%}


%include "std_string.i"
%include "std_vector.i"
// Instantiate templates used by example
namespace std {
  %template(IntVector) vector<int>;
  %template(IntMatrix) vector< vector<int> >;
  %template(DoubleVector) vector<double>;
  %template(StringVector) vector<string>;
}

%include "sequential_analysis.hpp"
