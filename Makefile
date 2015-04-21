F90=gfortran
LINK=$(F90)
LIBS=
LIBPATH=
FLAGS=-O3
#FLAGS=-g -Wall -Wtabs -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -fdiagnostics-color=always
EXE=main

all: $(EXE)

OBJS=main.o output.o integrate.o system.o kinds.o
main.o: output.o integrate.o system.o kinds.o
output.o: system.o kinds.o
integrate.o: system.o kinds.o
system.o: kinds.o

$(EXE): $(OBJS)
	@echo 'Linking [$(EXE)] from [$(OBJS)] using [$(LINK)]'
	@$(LINK) $(FLAGS) -o $(EXE) $(OBJS) $(LIBPATH) $(LIBS)

%.o: %.f90
	@echo 'Compiling [$@] from [$<] using [$(F90)]'
	@$(F90) $(FLAGS) -c $<

clean:
	@-rm $(EXE) $(OBJS) *.mod trajectory.*
