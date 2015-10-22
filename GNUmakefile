#-------------------------------------------------
# Author: ZhuHaiWen
# Email : zhuhw@ihep.ac.cn
# Last modified: 2014-08-25 10:16
# Filename: GNUmakefile
# Phone : 18627814293: 
#-------------------------------------------------

TARGET := main.exe

HEAD := ./include
HEADFILE := $(wildcard ${HEAD}/*.h)

SRC := ./src
SRCFILE := $(wildcard ${SRC}/*.cpp)

TEMP := ./temp
OBJFILE := $(patsubst ${SRC}/%.cpp,${TEMP}/%.o,${SRCFILE})

CC := g++
CPPFLAGS := -I$(HEAD)
LDLIBS :=
FFT3 := -lm -lstdc++ -lfftw3
# ---------------For Root----------------
# CPPFLAGS += $(shell root-config --cflags)
# LDLIBS += $(shell root-config --glibs --libs)
# --------------------------------------------
# ---------------For Geant4--------------------------------------------
# CPPFLAGS += $(shell geant4-config --cflags-without-gui) -I${G4INCLUDE}
# LDLIBS += $(shell geant4-config --libs-without-gui)
#----------------------------------------------------
.PHONY : all clean

all : $(TARGET)

clean:
	rm  -rf $(OBJFILE)
	rm -rf $(TARGET)
	rm -rf main.o
	rm -rf *~
	
$(TEMP)/%.o : $(SRC)/%.cpp $(HEADFILE)
	$(CC) -o $@ -c $< $(CPPFLAGS)

main.o : main.cpp $(HEADFILE)
	$(CC) -o $@ -c $< $(CPPFLAGS)

$(TARGET) : main.o $(OBJFILE)
	$(CC) -o $@ $^ $(LDLIBS) $(FFT3)
	rm -rf *.o
