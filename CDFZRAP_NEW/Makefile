CXX = c++

all:
	$(CXX) -c combine.cpp $(shell pkg-config lhapdf --cflags) $(shell applgrid-config --cxxflags)
	$(CXX) combine.o $(shell pkg-config lhapdf --cflags --libs) $(shell applgrid-config --cxxflags --ldflags) -o combine
