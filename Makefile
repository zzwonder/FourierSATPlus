CC     = g++
CFLAGS = -static -I cudd/cudd -I cudd/util -I cudd/ -I $(HPP_DIR) $(CUDD_INCLUSIONS) -std=c++11

#LFLAGS = -L cudd/cudd/.libs/ -lcudd -lm
LFLAGS = cudd/cudd/.libs/libcudd.a -lm
CUDD_DIR = cudd
CUDD_INCLUSIONS = -I $(CUDD_DIR) -I $(CUDD_DIR)/cudd -I $(CUDD_DIR)/epd -I $(CUDD_DIR)/mtr -I $(CUDD_DIR)/st

HPP_DIR = include
CPP_DIR = src
OBJ_DIR = obj
_OBJ = read_formula.o bdd_gradient.o
_OBJ_MAIN = $(_OBJ) main.o
OBJ = $(patsubst %, $(OBJ_DIR)/%, $(_OBJ_MAIN))

fsp: $(OBJ)
	$(CC) -o $@ $^ $(LFLAGS) && echo

$(OBJ_DIR)/%.o: $(CPP_DIR)/%.cpp $(HPP_DIR)/%.h 
	mkdir -p $(OBJ_DIR) && $(CC) -c -o $@ $< $(CFLAGS) && echo

