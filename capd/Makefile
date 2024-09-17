CXX ?= g++
CFLAGS = -Wall -Wextra -O2

BUILD_DIR = build

# Should point to the capd_config binary that comes with the CAPD
# install. If this is already on your path the default value can be
# used. Otherwise you need to set this when calling make.
CAPD_CONFIG ?= capd-config

TARGETS = build/ginzburg build/ginzburg-curve

all: directories $(TARGETS)

build/ginzburg: src/ginzburg.cpp directories
	$(CXX) $(CFLAGS) $< -o $@ `$(CAPD_CONFIG) --cflags --libs`

build/ginzburg-curve: src/ginzburg-curve.cpp directories
	$(CXX) $(CFLAGS) $< -o $@ `$(CAPD_CONFIG) --cflags --libs`

directories:
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR) $(TARGETS)
