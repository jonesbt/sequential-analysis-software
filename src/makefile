.DEFAULT_GOAL = lib

CXX = clang++
CXXFLAGS = -Wall
BOOST_TEST_CXXFLAGS = -DBOOST_TEST_DYN_LINK
BOOST_TEST_LDFLAGS = -L/usr/local/lib -lboost_unit_test_framework

clean:
	rm -f libsequentialanalysis.so
	rm -f test_sequential_analysis

lib:
	$(CXX) $(CXXFLAGS) sequential_analysis.cpp $(LDFLAGS) -shared -o libsequentialanalysis.so

test:
	$(CXX) $(CXXFLAGS) $(BOOST_TEST_CXXFLAGS) -c -o sequential_analysis.cpp.o sequential_analysis.cpp
	$(CXX) $(CXXFLAGS) $(BOOST_TEST_CXXFLAGS) -c -o test_sequential_analysis.cpp.o test_sequential_analysis.cpp
	$(CXX) $(LDFLAGS) $(BOOST_TEST_LDFLAGS) test_sequential_analysis.cpp.o sequential_analysis.cpp.o -o test_sequential_analysis

