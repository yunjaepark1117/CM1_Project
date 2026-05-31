CXX := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -I/opt/homebrew/include
SFML_LIBS := -L/opt/homebrew/lib -lsfml-graphics -lsfml-window -lsfml-system

.PHONY: all clean

all: sim spiral batch merge_data

sim: sim.cpp Body.h Vec2.h GravitySimulator.h InitialCondition.h DistributionGrid.h SpiralAnalysis.h
	$(CXX) $(CXXFLAGS) sim.cpp $(SFML_LIBS) -o sim

spiral: spiral.cpp Body.h Vec2.h GravitySimulator.h InitialCondition.h DistributionGrid.h SpiralAnalysis.h
	$(CXX) $(CXXFLAGS) spiral.cpp $(SFML_LIBS) -o spiral

batch: batch.cpp Body.h Vec2.h GravitySimulator.h InitialCondition.h DistributionGrid.h SpiralAnalysis.h
	$(CXX) $(CXXFLAGS) batch.cpp $(SFML_LIBS) -o batch

merge_data: merge_data.cpp
	$(CXX) $(CXXFLAGS) merge_data.cpp -o merge_data

clean:
	rm -f sim spiral batch merge_data
