INCS = -I$(HOME)/Programs/capd/include/ -I$(CURDIR)/src/

CC = g++
CFLAGS = -Wall -Wextra -O2

SRC_DIR = src
BUILD_DIR = build

TARGETS = build/ginzburg build/ginzburg-jacobian

# Need to do
# export LD_LIBRARY_PATH=/home/joeldahne/Programs/capd/lib/:$LD_LIBRARY_PATH

all: directories $(TARGETS)

build/ginzburg: src/ginzburg.cpp
	$(CC) $(CFLAGS) $^ -o $@  `/home/joeldahne/Programs/capd/bin/capd-config --cflags --libs`

build/ginzburg-jacobian: src/ginzburg-jacobian.cpp
	$(CC) $(CFLAGS) $^ -o $@  `/home/joeldahne/Programs/capd/bin/capd-config --cflags --libs`

directories:
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR) $(TARGETS)
