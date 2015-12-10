if [[-z "$CXX"]]; then CXX='c++'
CXX=clang++-3.5
CFLAGS=-c -std=c++11 
LIBS=-lboost_program_options
VPATH=./src
OBJECTS3 = net3.o
OBJECTS4 = net4.o

all: net3 net4

net3: $(OBJECTS3)
	$(CXX) $(OBJECTS3) -o net3 $(LIBS) 

net4: $(OBJECTS4)
	$(CXX) $(OBJECTS4) -o net4 $(LIBS) 

%.o: %.c++
	$(CXX) $(CFLAGS) $<

clean:
	rm -f *.o 
