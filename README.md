# CM1_Project
Archive for simulation codes used in Classical Mechanics 1 2026-1 Term Project


## Requirements

- C++17 or higher version
- Make
- SFML 3
- macOS Homebrew

SFML is needed for simulation display and saving snapshots.

```bash
brew install sfml
```

## How to compile:

## Build

## To compile all programs:

```bash
make
```

## Four programs are compiled as a result:
- `sim`: N-body simulator + spiral analyzer with SFML display
- `spiral`: N-body simulator + spiral analyzer without SFML display, simulation logs are saved at ./spiralSimulLogs1, ... ,./spiralSimulLogs4 according to the configuration chosen. Also saves key data in ./spiralSimulLogs*/analysis_index.csv
- `batch`: runs spiral simulator with multiple initial conditions automatically
- `merge_data`: merges data in ./spiralSimulLogs*/analysis_index.csv into ./Data/spiralSimulLogs1_merged.csv, ..., ./Data/spiralSimulLogs4_merged.csv

## To Run a compile & program:
```bash
make
./sim
```

```bash
make
./spiral
```

...

## The breif comment regarding the phenomenon mentioned in section 3-2 of the project paper can be found at 
./comment.txt




## Refer to the term project paper for more details!
