# makefile for artific_dif
SRC_FILES = $(wildcard *.f90)
OBJ_FILES = $(addprefix bin/,$(notdir $(SRC_FILES:.f90=.o)))

./bin/artific_dif: $(OBJ_FILES)
	gfortran -o $@ $^

./bin/%.o: %.f90
	gfortran -c -o $@ $<

# Utility targets
clean:
	rm -f ./bin/*.o
