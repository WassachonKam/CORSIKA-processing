ccsrc = $(wildcard *.cpp)
obj = $(ccsrc:.cpp=.o)

LDFLAGS = -std=c++11 -lm
CXXFLAGS =  -O0 -fbounds-check -ggdb -Wall -lz

corsikaReader: $(obj)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

.PHONY: clean
clean:
	rm $(obj) corsikaReader

