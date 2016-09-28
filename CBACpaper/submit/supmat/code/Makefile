main  : main.cpp potter.o driver.o looker.o params.o models.o viewer.o shaker.o keeper.o linker.o bumper.o sorter.o fixers.o util.o cell.o data.o geom.o util.hpp cell.hpp geom.hpp Vec.hpp Mat.hpp
	c++ -ggdb main.cpp -o main potter.o driver.o looker.o params.o models.o viewer.o shaker.o keeper.o linker.o bumper.o sorter.o fixers.o util.o cell.o data.o geom.o -lGL -lGLU -lglut -lpthread -lm


potter.o : potter.cpp util.hpp data.hpp geom.hpp Vec.hpp Mat.hpp
	c++ -c -ggdb potter.cpp

driver.o : driver.cpp util.hpp data.hpp geom.hpp Vec.hpp Mat.hpp
	c++ -c -ggdb driver.cpp

looker.o : looker.cpp util.hpp data.hpp geom.hpp Vec.hpp Mat.hpp
	c++ -c -ggdb looker.cpp

params.o : params.cpp util.hpp data.hpp geom.hpp Vec.hpp Mat.hpp
	c++ -c -ggdb params.cpp

models.o : models.cpp util.hpp data.hpp cell.hpp geom.hpp Vec.hpp Mat.hpp
	c++ -c -ggdb models.cpp

viewer.o : viewer.cpp util.hpp data.hpp cell.hpp geom.hpp Vec.hpp Mat.hpp
	c++ -c -ggdb viewer.cpp

shaker.o : shaker.cpp util.hpp data.hpp cell.hpp geom.hpp Vec.hpp Mat.hpp
	c++ -c -ggdb shaker.cpp

keeper.o : keeper.cpp util.hpp data.hpp cell.hpp geom.hpp Vec.hpp Mat.hpp
	c++ -c -ggdb keeper.cpp

linker.o : linker.cpp util.hpp data.hpp cell.hpp geom.hpp Vec.hpp Mat.hpp
	c++ -c -ggdb linker.cpp

bumper.o : bumper.cpp util.hpp data.hpp cell.hpp geom.hpp Vec.hpp Mat.hpp
	c++ -c -ggdb bumper.cpp

sorter.o : sorter.cpp util.hpp data.hpp cell.hpp geom.hpp Vec.hpp Mat.hpp
	c++ -c -ggdb sorter.cpp

fixers.o : fixers.cpp util.hpp data.hpp cell.hpp geom.hpp Vec.hpp Mat.hpp
	c++ -c -ggdb fixers.cpp

cell.o : cell.cpp util.hpp data.hpp cell.hpp geom.hpp Vec.hpp Mat.hpp
	c++ -c -ggdb cell.cpp

data.o : data.cpp util.hpp geom.hpp Vec.hpp Mat.hpp
	c++ -c -ggdb data.cpp

geom.o : geom.cpp util.hpp geom.hpp Vec.hpp Mat.hpp
	c++ -c -ggdb geom.cpp

util.o : util.cpp util.hpp
	c++ -c -ggdb util.cpp
