Sequential Analysis Software
============================

This repository contains 2 implementations of a Bayesian sequential procedure to
use when operating particle-tracking models to simulate multiphase flows. The
first implementation is in R. The second implemenation is writted in C++ and
includes a SWIG configuration file to generate bindings for other languages. The
two implementations may give minor differences in the results due to differences
between the returned values from rmultinom() in R and gsl_ran_multinomial() from
the GNU Scientific Library.

## R implementation
The R implementation may be installed using the below command within R. It
depends only on the gtools library, which may be obtained from CRAN.
```
install.packages('gtools')
devtools::install_github('btjones16/sequential-analysis-software')
```

## C++ implementation
The C++ implementation consists of 1 source code file, 1 header file, a
makefile, and a SWIG file. In addition, a second source code file is provided
that implements unit tests for the library.

### Installation
##### Linux / OS X
Paste the below command at a Terminal prompt. The script explains any changes
that it will make, then pauses prior to installing any new software.
```
curl -O -fsSL https://raw.githubusercontent.com/btjones16/sequential-analysis-software/master/install.sh && /bin/bash install.sh
```

#### Windows
For now, we recommend that Windows use the R package unless they are familiar
with compiling C++ code. If you wish to install the C++ library anyway, you'll
need to install the dependencies listed below, then modify src/makefile to work
with your system.

### Dependencies
- GNU Scientific Library (http://www.gnu.org/software/gsl/)
- BOOST.Math (http://www.boost.org/doc/libs/1_60_0/libs/math/doc/html/index.html)
- BOOST.Test (recommended, http://www.boost.org/doc/libs/1_60_0/libs/test/doc/html/index.html)

### Other languages
The C++ implementation includes a SWIG configuration to generate bindings for
use with other languages. For example, the below command will generate a Python
module that allows access to the public interface of the sequential analysis
routine. Generating bindings for other languages is similar. Please see the
SWIG documentation for details on modifying the target for other languages.

```
# Build the package.
cd src && make python
# Import the package and create a connectivity matrix.
python
import sequential_analysis_wrap
cm = sequential_analysis_wrap.ConnectivityMatrix(
    ['Origin 0', 'Origin 1'],
    ['Destination 0', 'Destination 1', 'Destination 2']
)
```

http://www.swig.org/
