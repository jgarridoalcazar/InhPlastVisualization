# $Id: Makefile,v 2.5 1998/03/17 23:14:51 jelasity Exp $

INTFO=$(patsubst usrintf/%.cc, usrintf/%.o, $(wildcard usrintf/*.cc))
UEGOO=$(patsubst uegolib/%.cc, uegolib/%.o, $(wildcard uegolib/*.cc))
UEGOH=$(wildcard include/*.h)

# -------------------------------------------------------------------
# ---- some macros to customize the optimizer and objective functions
# ---- additional libraries  won't be built just linked 

# PRJO defines the optimizers and the opbjective functions
# for search spaces; change if necessary
PRJO=projects/nrealopt.o projects/nbinval.o projects/nrealval.o

# PRJLIBS defines additional libraries that must be linked to the project
# defined by PRJO; change when needed
PRJLIBS=

# -------------------------------------------------------------------

CCC=mpic++
CFLAGS=-O2 -Iinclude 
CLIBS= -lm 

uego: $(INTFO) $(UEGOO) $(PRJO) 
	$(CCC) $(INTFO) $(UEGOO) $(PRJO) $(PRJLIBS) -o $@ $(CLIBS)

%.o: %.cc $(UEGOH)
	$(CCC) $(CFLAGS) -c $< -o $@

clean:
	rm uego $(UEGOO) $(INTFO)
	cd projects; rm *o
	
showdoc:
	nroff -man ../doc/ | more
	
