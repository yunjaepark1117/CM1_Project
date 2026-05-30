CXX := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -I/opt/homebrew/include
SFML_LIBS := -L/opt/homebrew/lib -lsfml-graphics -lsfml-window -lsfml-system

.PHONY: all clean

all: sim spiral batch

sim: sim.cpp Body.h Vec2.h GravitySimulator.h InitialCondition.h DistributionGrid.h SpiralAnalysis.h
	$(CXX) $(CXXFLAGS) sim.cpp $(SFML_LIBS) -o sim

spiral: spiral.cpp Body.h Vec2.h GravitySimulator.h InitialCondition.h DistributionGrid.h SpiralAnalysis.h
	$(CXX) $(CXXFLAGS) spiral.cpp $(SFML_LIBS) -o spiral

batch: batch.cpp Body.h Vec2.h GravitySimulator.h InitialCondition.h DistributionGrid.h SpiralAnalysis.h
	$(CXX) $(CXXFLAGS) batch.cpp $(SFML_LIBS) -o batch

clean:
	rm -f sim spiral batch
