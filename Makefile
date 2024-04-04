#
# Lightweight Matrix Library (LML) - Library Compiler
#
# File:     Makefile
# Author:   James Bray
# Repo:     https://github.com/James-Bray19/Lightweight-Matrix-Library
#
# This Makefile compiles the LML library for distribution 
# should there be any changes to the source code.
#

# compiler settings
CC = gcc
CFLAGS = -Wall -Wextra -pedantic -std=c11 -fPIC

# directories
SRC_DIR = src
INCLUDE_DIR = include
LIB_DIR = lib

# source files
SRC_FILES := $(wildcard $(SRC_DIR)/*.c)

# object files
OBJ_FILES := $(patsubst $(SRC_DIR)/%.c, $(LIB_DIR)/%.o, $(SRC_FILES))

# target library
DLL_NAME = liblml.dll
SO_NAME = liblml.so

TARGET_DLL = $(LIB_DIR)/$(DLL_NAME)
TARGET_SO = $(LIB_DIR)/$(SO_NAME)

# targets
all: $(TARGET_DLL) $(TARGET_SO)

# build the shared library (SO)
$(TARGET_SO): $(OBJ_FILES)
	$(CC) -shared -o $@ $^

# build the dynamic link library (DLL)
$(TARGET_DLL): $(OBJ_FILES)
	$(CC) -shared -o $@ $^

# compile source files
$(LIB_DIR)/%.o: $(SRC_DIR)/%.c | $(LIB_DIR)
	$(CC) $(CFLAGS) -c -o $@ $< -I$(INCLUDE_DIR)

# create the lib directory if it doesn't exist
$(LIB_DIR):
	mkdir -p $(LIB_DIR)

# clean generated files
clean:
	rm -rf $(LIB_DIR)/*.o $(TARGET_DLL) $(TARGET_SO)
