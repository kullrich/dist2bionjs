CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra

SRCS = njs.cpp bionjs.cpp dist2bionjs.cpp
OBJS = $(SRCS:.cpp=.o)
EXEC = dist2bionjs

.PHONY: all clean

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(EXEC)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(EXEC)
