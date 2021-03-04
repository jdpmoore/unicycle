
SRC=src
DST=obj
BLD=build

OBJ=$(patsubst %,$(DST)/%,axb.o edgmain.o bessj.o edgsh.o bessj0.o edgkern.o edgmoment.o bessj1.o edgmatinv.o edgsublay.o gemp.o edgmatrix.o getdata.o edghask.o edgbstab.o memcpy.o edgwvint.o edglayer.o edgpsv.o)

FFLAGS=-O3

.PHONY: all clean createdir bin

all: createdir bin

$(DST)/%.o: $(SRC)/%.f
	$(F77) $(FFLAGS) -c $^ -o $@

$(DST):
	@echo creating directory $(DST)
	@mkdir -p $(DST)

$(BLD):
	@echo creating directory $(BLD)
	@mkdir -p $(BLD)

$(BLD)/edgrn: $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o $@

bin: $(BLD)/edgrn

createdir: $(DST) $(BLD)

clean:
	rm -f $(DST)/*.o $(BLD)/edgrn
	rm -rf $(DST) $(BLD)

