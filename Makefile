all: dgs

dgs: dgs.cpp error.cpp main.cpp problem.cpp prolongation.cpp residual.cpp restriction.cpp utils.cpp vcycle.cpp correction.cpp uzawa.cpp iuzawa.cpp cg.cpp project.h
	g++ -o dgs -g dgs.cpp error.cpp main.cpp problem.cpp prolongation.cpp residual.cpp restriction.cpp utils.cpp vcycle.cpp correction.cpp uzawa.cpp iuzawa.cpp cg.cpp -lm -fopenmp -O2
