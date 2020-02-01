all: st mt

st: dgs.cpp error.cpp main.cpp problem.cpp prolongation.cpp residual.cpp restriction.cpp utils.cpp vcycle.cpp correction.cpp uzawa.cpp iuzawa.cpp cg.cpp pcg.cpp pcg_uzawa.cpp vcycle_precondition.cpp project.h
	g++ -o st -g dgs.cpp error.cpp main.cpp problem.cpp prolongation.cpp residual.cpp restriction.cpp utils.cpp vcycle.cpp correction.cpp uzawa.cpp iuzawa.cpp cg.cpp pcg.cpp pcg_uzawa.cpp vcycle_precondition.cpp -lm -O2
mt: dgs.cpp error.cpp main.cpp problem.cpp prolongation.cpp residual.cpp restriction.cpp utils.cpp vcycle.cpp correction.cpp uzawa.cpp iuzawa.cpp cg.cpp pcg.cpp pcg_uzawa.cpp vcycle_precondition.cpp project.h
	g++ -o mt -g dgs.cpp error.cpp main.cpp problem.cpp prolongation.cpp residual.cpp restriction.cpp utils.cpp vcycle.cpp correction.cpp uzawa.cpp iuzawa.cpp cg.cpp pcg.cpp pcg_uzawa.cpp vcycle_precondition.cpp -lm -fopenmp -O2

clean:
	rm st mt