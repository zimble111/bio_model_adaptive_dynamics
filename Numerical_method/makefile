CXX = g++
CXXFLAGS = -std=c++14 -Wall -pthread -Wextra -pedantic -Wconversion -g3 -march=native -O2 -isystem /Users/pavelgnilomedov/include/eigen-3.4-rc1 -isystem src/clipp.h -isystem /Users/pavelgnilomedov/gsl-master -isystem /usr/local/opt/gsl/include

SRC = src/
build_dir = build/
LIBS += -lpthread -L/usr/local/opt/gsl/lib -lgsl -lgslcblas

vpath %.cc $(SRC)
vpath %.h $(SRC)

all: $(build_dir)a

$(build_dir)a: $(addprefix $(build_dir),main.o conv.o)
	@mkdir -p $(build_dir)
	$(CXX) $^ $(CXXFLAGS) $(LIBS) -o $@

$(build_dir)%.o: %.cc %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(build_dir)%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) -rf $(build_dir)*.o $(build_dir)a

.PHONY: all clean
