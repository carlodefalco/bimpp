#MATLAB           =/Applications/MATLAB704/bin/matlab

include ${PETSC_DIR}/bmake/common/base

advdiffts: advdiffts.o  chkopts
	-$(CLINKER) -o advdiffts advdiffts.o $(PETSC_LIB)
	$(RM) advdiffts.o

examples: 
	${MATLAB} -nodisplay -nojvm -nodesktop -nosplash -r "buildexamples;quit"

runexample: advdiffts examples
	mpirun -np 2 advdiffts -prefix nproc2/ -preload false -epsilon 1e-6 -maxts 100 -deltat 0.1 -ts_type beuler -log_summary
	cd nproc2
	$(MATLAB) -r "plotseries;quit"