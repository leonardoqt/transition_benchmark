ROOT_DIR=$(shell pwd)
ODIR  = $(ROOT_DIR)/obj
SDIR  = $(ROOT_DIR)/src

F90   = gfortran
FFLAG = -llapack
 
OBJ   = $(ODIR)/benchmark_system.o $(ODIR)/model_h.o

hop.x : $(SDIR)/main.f90 $(OBJ)
	$(F90) -I$(ODIR) -o $@ $^ $(FFLAG)

$(ODIR)/%.o : $(SDIR)/%.f90 | $(ODIR)/.
	$(F90) -c -J$(ODIR) -o $@ $< $(FFLAG)

%/. : 
	mkdir -p $(patsubst %/.,%,$@)
	
.PRECIOUS: %/.
.PHONY: clean

clean :
	rm -rf hop.x $(ODIR)
