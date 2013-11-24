EIGEN_DIR = $(HOME)/lib/eigen3.2.0/
ODEINT_DIR = $(HOME)/lib/boost_1_53_0/
STLIB_DIR = $(HOME)/lib/stlib/src
CC = g++
CLIB = -L/usr/lib64 -lnetcdf_c++
CINC = -I/usr/include -I. -I $(EIGEN_DIR) -I $(ODEINT_DIR) -I $(STLIB_DIR)
CFLAG = -O3 -msse2 -std=c++0x -fopenmp
MAIN = Main
EXE = run
ADDONS = Include.hh NumericalMethod.hh Grid.hh Dynamics.hh ProgVariable.hh\

$(EXE): $(MAIN).o 
	$(CC) $(CFLAG) $(CLIB) -o $(EXE) $(<)
$(MAIN).o: $(MAIN).cc $(ADDONS)
	$(CC) $(CFLAG) $(CINC) -c $(<)
clean:
	@ rm -f run
	@ rm -f *.o
	@ rm -f start.in
.PHONY: clean
