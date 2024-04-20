CPPFLAGS = -I../inc -DNDEBUG -I/home/galati_j/PACS_TEST/pacs-examples/Examples/include
CXX = g++ -std=c++20
#LDFLAGS ?= -L/home/galati_j/PACS_TEST/pacs-exampleS/Examples/lib
#LIBS  ?= -lmuparser


EXEC = main
SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o)
HEADER=$(wildcard *.hpp)


.phony = all clean distclean

.DEFAULT_GOAL = all

all: $(EXEC)

clean:
	$(RM) *.o

distclean: clean
	$(RM) main


$(EXEC): $(OBJS)
	$(CXX) $(OBJS) -o $(EXEC)

main.o: $(SRCS)
	$(CXX) $(SRCS) $(CPPFLAGS) -c main.o


MyMatrix.o: MyMatrix.cpp
	$(CXX) MyMatrix.cpp $(CPPFLAGS) -c MyMatrix.o




