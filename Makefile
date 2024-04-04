

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
LIB_NAME = liblml.so
TARGET_LIB = $(LIB_DIR)/$(LIB_NAME)

# targets
all: $(TARGET_LIB)

# build the shared library
$(TARGET_LIB): $(OBJ_FILES)
	$(CC) -shared -o $@ $^

# compile source files
$(LIB_DIR)/%.o: $(SRC_DIR)/%.c | $(LIB_DIR)
	$(CC) $(CFLAGS) -c -o $@ $< -I$(INCLUDE_DIR)

# create the lib directory if it doesn't exist
$(LIB_DIR):
	mkdir -p $(LIB_DIR)

# clean generated files
clean:
	rm -rf $(LIB_DIR)/*.o $(TARGET_LIB)