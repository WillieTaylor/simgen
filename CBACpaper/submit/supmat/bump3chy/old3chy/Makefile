sim3chy  : main/main.c driver.o bumper.o linker.o keeper.o shaker.o viewer.o utility.o main/common.h
	cc -ggdb main/main.c  -o sim3chy  driver.o bumper.o linker.o keeper.o shaker.o viewer.o utility.o util/wt/util.o util/wt/geom.o util/wt/sort.o -lGL -lGLU -lX11 -lglut -lXext -lXmu -lXi -lSDL -lpthread -lm

driver.o : driver.c main/common.h
	cc -c -ggdb driver.c

bumper.o : main/bumper.c main/common.h
	cc -c -ggdb main/bumper.c

linker.o : main/linker.c main/common.h
	cc -c -ggdb main/linker.c

keeper.o  : main/keeper.c main/common.h
	cc -c -ggdb main/keeper.c

shaker.o  : main/shaker.c main/common.h
	cc -c -ggdb main/shaker.c

viewer.o : main/viewer.c main/common.h
	cc -c -ggdb main/viewer.c

utility.o : main/utility.c main/common.h
	cc -c -ggdb main/utility.c
