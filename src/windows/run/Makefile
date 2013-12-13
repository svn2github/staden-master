all: sprun.exe

CC=$(COMPILER)gcc
WINDRES=$(COMPILER)windres

run.res.o: sprun.rc
	$(WINDRES) -o sprun.res.o sprun.rc 

sprun.exe: main.cpp run.res.o
	$(CC) -g -mwindows main.cpp -o sprun.exe sprun.res.o
