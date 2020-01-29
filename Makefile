all: dgs

dgs: dgs.cpp error.cpp main.cpp problem.cpp prolongation.cpp residual.cpp restriction.cpp utils.cpp vcycle.cpp correction.cpp uzawa.cpp iuzawa.cpp cg.cpp pcg.cpp pcg_uzawa.cpp vcycle_precondition.cpp project.h
	g++ -o dgs -g dgs.cpp error.cpp main.cpp problem.cpp prolongation.cpp residual.cpp restriction.cpp utils.cpp vcycle.cpp correction.cpp uzawa.cpp iuzawa.cpp cg.cpp pcg.cpp pcg_uzawa.cpp vcycle_precondition.cpp -lm -fopenmp -O2
