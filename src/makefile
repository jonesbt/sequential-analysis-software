.DEFAULT_GOAL = null

CXX = clang++
LD = clang++
CXXFLAGS = -Wall -std=c++11 
LDFLAGS = 

CXXFLAGS += $(shell pkg-config --cflags gsl)
LDFLAGS += $(shell pkg-config --libs gsl)

BOOST_TEST_CXXFLAGS = -DBOOST_TEST_DYN_LINK
BOOST_TEST_LDFLAGS = -L/usr/local/lib -lboost_unit_test_framework

PYTHON_CXXFLAGS = $(shell pkg-config --cflags python)
PYTHON_LDFLAGS = $(shell pkg-config --libs python)

clean:
	rm -f *.o
	rm -f libsequentialanalysis.so
	rm -f test_sequential_analysis
	rm -f sequential_analysis.py
	rm -f sequential_analysis_wrap.*

lib:
	$(CXX) $(CXXFLAGS) sequential_analysis.cpp $(LDFLAGS) -shared -o libsequentialanalysis.so

null:
	echo ""

python:
	swig -v -c++ -python sequential_analysis.i
	$(CXX) $(CXXFLAGS) -c -o sequential_analysis.cpp.o sequential_analysis.cpp
	$(CXX) $(CXXFLAGS) $(PYTHON_CXXFLAGS) -c -o sequential_analysis_wrap.cxx.o sequential_analysis_wrap.cxx
	$(LD) $(LDFLAGS) $(PYTHON_LDFLAGS) -shared sequential_analysis.cpp.o sequential_analysis_wrap.cxx.o -o _sequential_analysis_wrap.so

r:
	swig -c++ -r -o sequential_analysis_wrap.cpp sequential_analysis.i
	PKG_LIBS="sequential_analysis.cpp $(LDFLAGS)" R CMD SHLIB sequential_analysis_wrap.cpp

test:
	$(CXX) $(CXXFLAGS) -c -o sequential_analysis.cpp.o sequential_analysis.cpp
	$(CXX) $(CXXFLAGS) $(BOOST_TEST_CXXFLAGS) -c -o test_sequential_analysis.cpp.o test_sequential_analysis.cpp
	$(LD) $(LDFLAGS) $(BOOST_TEST_LDFLAGS) test_sequential_analysis.cpp.o sequential_analysis.cpp.o -o test_sequential_analysis

