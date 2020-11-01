
SRC=src
DST=obj
LIB=lib
BLD=build

OBJ=$(patsubst %,$(DST)/%,edcdisc.o edcgrn.o edcmain.o edcoutput.o getdata.o okada.o dc3d.o dc3d0.o edcmp.o)
OBJ_LIB=$(patsubst %,$(DST)/%,edcdisc.o edcgrn.o edcoutput.o getdata.o okada.o dc3d.o dc3d0.o edcmp.o loadgrn.o)
OBJ_TEST=$(patsubst %,$(DST)/%,edcdisc.o edcgrn.o edcoutput.o getdata.o okada.o dc3d.o dc3d0.o edcmp.o loadgrn.o testedcmp.o)

FFLAGS=-O3 -fopenmp

.PHONY: all clean createdir library binary

all: createdir library binary

$(LIB):
	@echo creating directory $(LIB)
	@mkdir -p $(LIB)

$(DST):
	@echo creating directory $(DST)
	@mkdir -p $(DST)

$(BLD):
	@echo creating directory $(BLD)
	@mkdir -p $(BLD)

$(DST)/%.o: $(SRC)/%.f $(SRC)/edcglobal.h
	$(F77) $(FFLAGS) -c $(filter-out $(SRC)/edcglobal.h,$^) -o $@

$(DST)/%.o: $(SRC)/%.f90
	$(F77) $(FFLAGS) -c $^ -o $@

$(LIB)/libedcmp.a: $(OBJ_LIB)
	ar cr $@ $^

$(BLD)/edcmp: $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o $@

testedcmp: $(OBJ_TEST)
	$(FC) $(FFLAGS) $^ -o $(BLD)/testedcmp

library: $(LIB)/libedcmp.a

binary: $(BLD)/edcmp

createdir: $(DST) $(BLD) $(LIB)

clean:
	rm -f $(DST)/*.o $(BLD)/edcmp $(LIB)/libedcmp.a
	rm -rf $(DST) $(BLD) $(LIB)

