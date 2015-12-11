if [[-z "$CXX"]]; then CXX='c++'
CXX=clang++-3.5
CFLAGS=-c -std=c++11 
LIBS=-lboost_program_options
VPATH=./src
OBJECTS3 = utils.o network.o net3.o
OBJECTS4 = utils.o network.o net4.o
OBJECTS_DEND = utils.o utils_dendlatt.o randnet_dendritic.o
OBJECTS_LATT = utils.o utils_dendlatt.o randnet_lattice.o 

all: net3 net4 dend

net3: $(OBJECTS3)
	$(CXX) $(OBJECTS3) -o net3 $(LIBS) 

net4: $(OBJECTS4)
	$(CXX) $(OBJECTS4) -o net4 $(LIBS) 

dend: $(OBJECTS_DEND)
	$(CXX) $(OBJECTS_DEND) -o randnet_dendritic $(LIBS)

latt: $(OBJECTS_LATT)
	$(CXX) $(OBJECTS_LATT) -o randnet_lattice $(LIBS)

%.o: %.c++
	$(CXX) $(CFLAGS) $<

clean:
	rm -f *.o 
