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
INCLUDES = -I. -I/inc -I/home/cdrpico/Documents/cppProjects/TSSPInstancesReader/inc -I/opt/ibm/ILOG/CPLEX_Studio201/concert/include -I/opt/ibm/ILOG/CPLEX_Studio201/cplex/include -I/opt/gurobi910/linux64/include
LIBS = -L/opt/gurobi910/linux64/lib/ -lgurobi_c++ -lgurobi91 -L/home/cdrpico/Documents/cppProjects/TSSPInstancesReader -L/opt/ibm/ILOG/CPLEX_Studio201/concert/lib/x86-64_linux/static_pic/ -L/opt/ibm/ILOG/CPLEX_Studio201/cplex/lib/x86-64_linux/static_pic/ -lm -lconcert -lilocplex -lcplex -lpthread -ldl -lread_instances_TSSP 

DEFS= 

DEFINES= -DIL_STD $(INCLUDES) $(DEFS) -DSYS_UNIX=1 -std=c++0x -Wall
CFLAGS= -g $(DEFINES)
CPPFLAGS= -g $(DEFINES)
CXXFLAGS='-D_GLIBCXX_USE_CXX11_ABI=0'

#To compile as static library
SRCS = src/GAPM.cpp src/SFLP_GAPM.cpp src/UsfFunctions.cpp src/BendersAPM_CP.cpp src/BendersAPM_SFLP.cpp src/CapPlan_GAPM.cpp src/InstancesSFCMFP.cpp src/OuterBendersSFLP.cpp src/SFCMFP_GAPM.cpp main.cpp

#To compile as static library
OBJS = src/GAPM.o src/SFLP_GAPM.o src/UsfFunctions.o src/BendersAPM_CP.o src/BendersAPM_SFLP.o src/CapPlan_GAPM.o src/InstancesSFCMFP.o src/OuterBendersSFLP.o src/SFCMFP_GAPM.o main.o

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
