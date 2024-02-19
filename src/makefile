.DEFAULT_GOAL = null

CXX = g++
LD = g++
CXXFLAGS = -Wall -std=c++14
LDFLAGS = 

CXXFLAGS += $(shell pkg-config --cflags gsl) -fPIC
LDFLAGS += $(shell pkg-config --libs gsl)

BOOST_TEST_CXXFLAGS = -DBOOST_TEST_DYN_LINK
BOOST_TEST_LDFLAGS = -L/usr/local/lib -lboost_unit_test_framework

PYTHON_CXXFLAGS = $(shell pkg-config --cflags python3)
PYTHON_LDFLAGS = $(shell pkg-config --libs python3)

clean:
	rm -f *.o
	rm -f libsequentialanalysis.so
	rm -f test_sequential_analysis
	rm -f sequential_analysis.py
	rm -f sequential_analysis_wrap.*

lib:
	echo $(CXXFLAGS)
	echo $(LDFLAGS)
	$(CXX) $(CXXFLAGS) sequential_analysis.cpp $(LDFLAGS) -shared -o libsequentialanalysis.so

null:
	echo ""

python:
	swig -v -c++ -python sequential_analysis.i
	$(CXX) $(CXXFLAGS) -c -o sequential_analysis.cpp.o sequential_analysis.cpp
	$(CXX) $(CXXFLAGS) $(PYTHON_CXXFLAGS) -c -o sequential_analysis_wrap.cxx.o sequential_analysis_wrap.cxx
	$(LD) -shared sequential_analysis.cpp.o sequential_analysis_wrap.cxx.o $(LDFLAGS) $(PYTHON_LDFLAGS) -o _sequential_analysis_wrap.so

r:
	swig -c++ -r -o sequential_analysis_wrap.cpp sequential_analysis.i
	PKG_LIBS="sequential_analysis.cpp $(LDFLAGS)" R CMD SHLIB sequential_analysis_wrap.cpp

test:
	$(CXX) $(CXXFLAGS) -c -o sequential_analysis.cpp.o sequential_analysis.cpp
	$(CXX) $(CXXFLAGS) $(BOOST_TEST_CXXFLAGS) -c -o test_sequential_analysis.cpp.o test_sequential_analysis.cpp
	$(LD) test_sequential_analysis.cpp.o sequential_analysis.cpp.o $(BOOST_TEST_LDFLAGS) $(LDFLAGS) -o test_sequential_analysis

check: test
	./test_sequential_analysis
