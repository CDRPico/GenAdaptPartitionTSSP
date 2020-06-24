##--------------------------------
# Created by CDRPico
# 09/06/2020 23:06
## -------------------------------

rm = /bin/rm -f
CC = g++
#To compile as static library
#PROGNAME = read_instances_TSSP.a
#To compile as program
PROGNAME = GAPM
INCLUDES = -I. -I/inc -I/home/cramirezp/GAPM/TSSPInstancesReader/inc -I/share/apps/cplex/include/ -I/share/apps/cplex129/concert/include/ -I/share/apps/gurobi902/linux64/include 
LIBS = -L/home/cramirezp/GAPM/TSSPInstancesReader -L/share/apps/gurobi902/linux64/lib/ -lgurobi90 -lgurobi_c++ -L/share/apps/cplex129/concert/lib/x86-64_linux/static_pic/  -L/share/apps/cplex129/cplex/lib/x86-64_linux/static_pic/ -lm -lconcert -lilocplex -lcplex -lpthread -ldl -lread_instances_TSSP

DEFS= 

DEFINES= -DIL_STD $(INCLUDES) $(DEFS) -DSYS_UNIX=1 -std=c++0x -Wall
CFLAGS= -g $(DEFINES)
CPPFLAGS= -g $(DEFINES)
CXXFLAGS='-D_GLIBCXX_USE_CXX11_ABI=0'

#To compile as static library
SRCS = src/BendersAPM_SFLP.cpp src/GAPM.cpp src/SFLP_GAPM.cpp src/UsefulFunctions.cpp main.cpp
#To compile as program
#SRCS = src/UsefulFunctions.cpp src/GenerateInstanceSFLP.cpp src/InstanceSFLP.cpp main.cpp

#To compile as static library
OBJS = src/BendersAPM_SFLP.o src/GAPM.o src/SFLP_GAPM.o src/UsefulFunctions.o main.o
#To compile as program
#OBJS = src/UsefulFunctions.o src/GenerateInstanceSFLP.o src/InstanceSFLP.o main.o

.c.o:
	$(rm) $@
	$(CC) $(CFLAGS) -c $*.c

all: $(PROGNAME)

$(PROGNAME) : $(OBJS)
	#To compile as static library
	#ar rcs $(PROGNAME) $(OBJS) $(LIBS)
	#To compile as program
	$(CC) $(CFLAGS) -o $(PROGNAME) $(OBJS) $(LIBS)
clean:
	$(rm) $(OBJS) $(PROGNAME) core *~
