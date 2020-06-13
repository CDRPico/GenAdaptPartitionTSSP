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
INCLUDES = -I. -I/inc -I/home/cdrpico/cppProjects/tsspInstancesReader/TSSPInstancesReader/inc -I/opt/ibm/ILOG/CPLEX_Studio1210/concert/include -I/opt/ibm/ILOG/CPLEX_Studio1210/cplex/include
LIBS = -L/home/cdrpico/cppProjects/tsspInstancesReader/TSSPInstancesReader -L/opt/ibm/ILOG/CPLEX_Studio1210/concert/lib/x86-64_linux/static_pic/ -L/opt/ibm/ILOG/CPLEX_Studio1210/cplex/lib/x86-64_linux/static_pic/ -lm -lconcert -lilocplex -lcplex -lpthread -ldl -lread_instances_TSSP

DEFS= 

DEFINES= -DIL_STD $(INCLUDES) $(DEFS) -DSYS_UNIX=1 -std=c++0x -Wall
CFLAGS= -g $(DEFINES)
CPPFLAGS= -g $(DEFINES)
CXXFLAGS='-D_GLIBCXX_USE_CXX11_ABI=0'

#To compile as static library
SRCS = src/GAPM.cpp src/SFLP_GAPM.cpp src/UsefulFunctions.cpp main.cpp
#To compile as program
#SRCS = src/UsefulFunctions.cpp src/GenerateInstanceSFLP.cpp src/InstanceSFLP.cpp main.cpp

#To compile as static library
OBJS = src/GAPM.o src/SFLP_GAPM.o src/UsefulFunctions.o main.o
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
