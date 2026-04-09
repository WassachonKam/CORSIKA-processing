# ccsrc = $(wildcard *.cpp) cnpy/cnpy.cpp
# obj   = $(ccsrc:.cpp=.o)

# CXXFLAGS = -std=c++11 -O0 -fbounds-check -ggdb -Wall -Icnpy
# LDFLAGS  = -lm -lz   # keep zlib


# corsikaReader: $(obj)
# 	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# .PHONY: clean
# clean:
# 	rm -f $(obj) corsikaReader





# รายชื่อไฟล์ Source Code
ccsrc = corsikaReader.cpp cnpy/cnpy.cpp

# ค่า Default
PARTICLE ?= proton
E_VAL    ?= lgE_16.0
ZENITH   ?= sin2_0.0

# Flags
# -I. และ -I./cnpy เพื่อให้หา Header files เจอ
CXXFLAGS += -O3 -Wall -I. -I./cnpy
CXXFLAGS += -DPRIMARY=\"$(PARTICLE)\" -DENERGY=\"$(E_VAL)\" -DTHETA=\"$(ZENITH)\"
LIBS      = -lz

# Target หลัก
corsikaReader: $(ccsrc)
	$(CXX) $(CXXFLAGS) $(ccsrc) -o $@ $(LIBS)

clean:
	rm -f corsikaReader