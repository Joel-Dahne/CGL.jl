CXX ?= g++
CFLAGS = -Wall -Wextra -O2

BUILD_DIR = build

# Should point to the capd_config binary that comes with the CAPD
# install. If this is already on your path the default value can be
# used. Otherwise you need to set this when calling make.
CAPD_CONFIG ?= capd-config

# Get the CAPD library path
CAPD_LIB_PATH=$(shell $(CAPD_CONFIG) --libs-only-L | sed 's/-L//g' | xargs)

TARGETS = build/Q_zero

all: directories $(TARGETS)

build/Q_zero: src/Q_zero.cpp | directories
	$(CXX) $(CFLAGS) $< -o $@ `$(CAPD_CONFIG) --cflags --libs` -Wl,-rpath,$(CAPD_LIB_PATH)

directories:
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR) $(TARGETS)
