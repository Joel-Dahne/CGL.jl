INCS = -I$(HOME)/Programs/capd/include/ -I$(CURDIR)/src/

CC ?= g++
CFLAGS = -Wall -Wextra -O2

SRC_DIR = src
BUILD_DIR = build
# Should point to the capd_config binary that comes with the CAPD
# install. If this is already on your path the default value can be
# used. Otherwise you need to set this when calling make.
CAPD_CONFIG ?= capd-config

TARGETS = build/ginzburg

all: directories $(TARGETS)

build/ginzburg: src/ginzburg.cpp
	$(CC) $(CFLAGS) $^ -o $@ `$(CAPD_CONFIG) --cflags --libs`

directories:
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR) $(TARGETS)
