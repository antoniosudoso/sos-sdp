matlab_path = /usr/local/MATLAB/R2019b
sources = main.cpp matlab_util.cpp sdp_solver_util.cpp sdp_branch_and_bound.cpp kmeans_util.cpp Kmeans.cpp JobQueue.cpp util.cpp ThreadPool.cpp
objects = $(sources:.cpp=.o)

all: bb

bb: $(objects)
	g++ -L $(matlab_path)/extern/bin/glnxa64/ -o $@ $^ -larmadillo -pthread -llapack -lopenblas -lMatlabDataArray -lMatlabEngine -Wl,-rpath,$(matlab_path)/extern/bin/glnxa64/

%.o: %.cpp
		g++ -Wall -O2 -o $@ -c $< -I $(matlab_path)/extern/include/


.PHONY: clean
		
clean:
		rm -f bb *.o
