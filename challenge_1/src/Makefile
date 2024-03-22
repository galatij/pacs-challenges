CPPFLAGS = -I../inc -DNDEBUG -I/home/galati_j/PACS_TEST/pacs-examples/Examples/include
CXX = g++ -std=c++17
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


gradient_method.o: parameters.cpp read_parameters.cpp gradient_method.cpp
	$(CXX) parameters.cpp read_parameters.cpp gradient_method.cpp $(CPPFLAGS) -c gradient_method.o

compute_gradients.o: compute_gradients.cpp derivatives.cpp
	$(CXX) compute_gradients.cpp derivatives.cpp $(CPPFLAGS) -c compute_gradients.o

read_parameters.o: read_parameters.cpp parameters.cpp
	$(CXX) parameters.cpp read_parameters.cpp $(CPPFLAGS) -c read_parameters.o

parameters.o: parameters.cpp
	$(CXX) parameters.cpp -I../inc -c parameters.o



