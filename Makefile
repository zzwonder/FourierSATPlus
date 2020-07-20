CC     = g++
CFLAGS = -static -I cudd/cudd -I cudd/util -I cudd/ -I $(HPP_DIR) $(CUDD_INCLUSIONS) $(DLIB_INCLUSIONS) -std=c++11

#LFLAGS = -L cudd/cudd/.libs/ -lcudd -lm
LFLAGS = cudd/cudd/.libs/libcudd.a -lm -lpthread -lX11 -std=c++11
CUDD_DIR = cudd
CUDD_INCLUSIONS = -I $(CUDD_DIR) -I $(CUDD_DIR)/cudd -I $(CUDD_DIR)/epd -I $(CUDD_DIR)/mtr -I $(CUDD_DIR)/st
DLIB_INCLUSIONS = -I dlib

HPP_DIR = include
CPP_DIR = src
OBJ_DIR = obj
DLIB_CPP = dlib/dlib/all/source.cpp
_OBJ = read_formula.o bdd_gradient.o optimize.o
_OBJ_MAIN = $(_OBJ) main.o

OBJ = $(patsubst %, $(OBJ_DIR)/%, $(_OBJ_MAIN))

fsp: $(OBJ)
	$(CC) -o $@ $(DLIB_CPP) $^ $(LFLAGS) && echo

$(OBJ_DIR)/%.o: $(CPP_DIR)/%.cpp $(HPP_DIR)/%.h
	mkdir -p $(OBJ_DIR) && $(CC) -c -o $@ $< $(CFLAGS) && echo

