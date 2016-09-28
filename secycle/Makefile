simcell  : driver.o sims/data.hpp sims/geom.hpp sims/Vec.hpp sims/Mat.hpp
	c++ -ggdb sims/main.cpp -o simgen driver.o sims/potter.o sims/params.o sims/models.o sims/viewer.o sims/shaker.o sims/keeper.o sims/linker.o sims/bumper.o sims/sorter.o sims/fixers.o sims/cell.o sims/data.o sims/geom.o sims/util.o -lGL -lGLU -lglut -lpthread -lm 

driver.o : driver.cpp
	c++ -c -ggdb driver.cpp
