# district_mbo

This is code simulates volume-preserving motion by mean curvature on U.S. state maps. In other words, it smooths district boundaries of maps while keeping the population in each district roughly equal. The flow can be started from either a prescribed set of district lines or random. A fuller description of the method involved can be found in a paper to appear on Arxiv. 


## Getting Started

We use a Python wrapper to call the C executable where the flow is computed. 

### C Dependencies

The C code requires [VLFEAT](http://www.vlfeat.org) to compile. Instructions for installing from the command line can be found [here](http://www.vlfeat.org/install-shell.html). The directory where the source is unpacked is VLFEATROOT, and should be added to your .bash_profile. The compand to compile the district_mbo.c script on a Mac is

```
cc -O3 -g -Wall -o district_mbo district_mbo.c -I$VLFEATROOT -lvl -L/$VLFEATROOT/bin/maci64
```

### Python set-up

The required Python packages are
* Shapely
* imageio 
* numpy
* Pillow
* matplotlib
* pyshp

They're listed in requirements.txt and can be installed with

```
pip install requirements.txt
```

## Demos 

If you're interested in seeing what the flow looks like from random, or using the algorithm to compute the "compactness energy" for a given districting, there are several demos you can try. 

### Making a GIF
To make a gif showing the flow from random, call 

```
python demo_gif.py [2-Letter State Abbreviation]
```

e.g.,

```
python demo_gif.py VA
```

This will save an output gif to VA_output.gif. You can also pass along additional parameters to generate the flow with different properties or from an existing map. The full argument list is 

```
python demo_gif.py [state] [filename.shp (in data/)] [# nearest neighbors] [centroid coefficient] [temperature] [annealing rate]
```
where 
* state : The state you want to run
* filename.shp : (Optional) The mapping you want to initialize from (without this, you initialize from random). This should be in the data/ folder.
* nearest neighbor count : (Optional) Number of neighboring units considered. Generally should be > 100
* centroid coefficient : (Optional) Term to penalize split districts; increasing it reduces the chance of districts splitting. Generally should be > 1. 
* temperature : (Optional) Randomness term. Default set to 0.1
* annealing rate : (Optional) Multiplier that decreases temperature each iteration. Default is 0.95.


### Energy Histogram
To make a histogram comparing a districting's compactness energy to the energies of the local minima arrived at by the algorithm, 

```
python demo_energy_hist.py [2-Letter State Abbreviation] [Districting Shapefile]
```

e.g.,

```
python demo_energy_hist.py PA "Remedial Plan Shapefile.shp"
```

where Remedial Plan Shapefile.shp is in the data/ folder. This will save an image of the histogram comparing the energies from 30 runs of the algorithm to the energy in NewMapping.shp to VA_energy_hist.png


gif to VA_output.gif. You can also pass along additional parameters to generate the flow with different properties or from an existing map. The full argument list is 

```
python demo_energy_hist.py [state] [filename.shp (in data/)] [# nearest neighbors] [centroid coefficient] [num runs]

```
where 
* state : The state you want to run
* filename.shp : (Optional) The mapping you want to initialize from. This should be in the data/ folder.
* nearest neighbor count : (Optional) Number of neighboring units considered. Generally should be > 100
* centroid coefficient : (Optional) Term to penalize split districts; increasing it reduces the chance of districts splitting. Generally should be > 1. 
* number of runs : (Optional) Number of times to run the algorithm to produce the comparison scores. Default is 30.


### Flow
To make a static image showing flow, call

```
python demo_flow.py [2-Letter State Abbreviation]
```

e.g.,

```
python demo_flow.py PA
```

This will save an image showing 8 steps of the flow, initialized from the current mapping. The full comand syntax is 


```
python demo_nearby_maps.py [state] [filename.shp (in data/)] [# nearest neighbors] [centroid coefficient] 

```
where 
* state : The state you want to run
* filename.shp : (Optional) The mapping you want to initialize from. This should be in the data/ folder.
* nearest neighbor count : (Optional) Number of neighboring units considered. Generally should be > 100
* centroid coefficient : (Optional) Term to penalize split districts; increasing it reduces the chance of districts splitting. Generally should be > 1. 

### Figures

Code to reproduce all the figures is in make_figures.py. You can run them with 

```
python make_figures.py
```

If you don't want to generate all the figures, you can avoid out the ones you don't want by opening make_figures.py and commenting out the calls to generate any figure you're not interested in. 

## Data

All the data is taken from openly available datasets. 

* [U.S. Census district data](https://www.census.gov/geo/maps-data/data/cbf/cbf_cds.html)
* [U.S. Census tract data](https://www.census.gov/geo/maps-data/data/cbf/cbf_tracts.html)
* [American Community Survey (Population Data](http://www.census.gov/programs-surveys/acs/data.html)
* [PA Remedial Mapping 2018](http://www.pacourts.us/news-and-statistics/cases-of-public-interest/league-of-women-voters-et-al-v-the-commonwealth-of-pennsylvania-et-al-159-mm-2017)

## Authors

* **Olivia Walch** -- [Website](http://www.oliviawalch.com)
* **Matt Jacobs** - [Website](http://www.math.ucla.edu/~majaco/)


## License

This software is open source and under an MIT license. 

