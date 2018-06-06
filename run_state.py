import shapefile
from shapely.geometry import shape
import csv
import sys
import matplotlib.pyplot as plt
import numpy as np
from random import randint
from subprocess import call
from array import array
from shapely.geometry import Polygon
from shapely.geometry.multipolygon import MultiPolygon
import glob
import math
import visualize_maps
import config

## Reads Congressional district files
def read_districts(state,alt_map):
    
    # If an alternative map is provided
    if(len(alt_map) > 0):
        district_data = shapefile.Reader('data/' + alt_map)
        districts = []
        for tract in district_data.shapeRecords():
            district_poly = shape(tract.shape.__geo_interface__)
            districts.append(district_poly)
        
        return districts
    
    # Otherwise, default to 2016 Census data for districting
    district_data = shapefile.Reader('data/cb_2016_us_cd115_500k/cb_2016_us_cd115_500k.shp')

    districts = []
    for tract in district_data.shapeRecords():
        if(tract.record[0] == config.state_dict[state]):
            district_poly = shape(tract.shape.__geo_interface__)
            districts.append(district_poly)

    return districts

## Finds district containing a tract
def shape_to_district(tract,districts):
    ind = -1
    shp_geom = shape(tract.shape.__geo_interface__)
    centroid = shp_geom.centroid

    for i in range(0,len(districts)):
        
        if(districts[i].contains(centroid)):
            ind = i
            break

    return ind

## Gets shapefile for state
def data_for_state(state):
    
    # Can be replaced with 2017 data
    return shapefile.Reader('data/cb_2016_' + config.state_dict[state] + '_tract_500k/cb_2016_' + config.state_dict[state] +'_tract_500k.shp')



## Loads all relevant data for a state
def read(state, alt_map):
    
    print 'Reading ' + state + ' data...'

    shape_data = data_for_state(state)
    
    id_to_centroids = {} # GEOID -> centroid dictionary
    id_to_district = {} # GEOID -> district index in state dictionary

    num_dists = 0
    district_keys = {} # district index in US list -> district index in state

    # Load district data for state
    districts = read_districts(state,alt_map)

    # For every tract in shape_data, find containing district:
    for tract in shape_data.shapeRecords():
        aff_geo_id = tract.record[3] # GEOID
        shp_geom = shape(tract.shape.__geo_interface__)
        
        id_to_centroids[aff_geo_id] = shp_geom.centroid # Store centroid

        ind = shape_to_district(tract,districts) # Index of containing district in US
        
        if(ind == -1): # If the centroid failed to get a district
            id_to_district[aff_geo_id] = 0
        else:
            if(ind in district_keys):  # If district already re-mapped to in-state index
                id_to_district[aff_geo_id] = district_keys[ind]
            else: # Otherwise, create a new rank and add to the dictionary
                district_keys[ind] = num_dists
                id_to_district[aff_geo_id] = district_keys[ind]
                num_dists = num_dists + 1


    # Extract data for census units + containing districts
    unit_data = []
    unit_districts = []

    num_valid_units = 0
    tot_pop = 0
    district_pops = [0] * num_dists
    
    # American Community Survey (ACS) data from Census
    # Contains population data for entire US.
    with open('data/ACS_16_5YR_DP05_with_ann.csv') as csvDataFile:
        csv_reader = csv.reader(csvDataFile)
        
        for row in csv_reader:
            key = row[0] # GEOID
            if key in id_to_centroids: # If we have spatial data for this key...
                
                centroid = id_to_centroids[key] # Get centroid
                
                pop = int(row[3]) # Get population
                
                # Append GEOID, centroid, and population to array
                unit_data.append(int(key[config.id_0:]))
                unit_data.append(centroid.x)
                unit_data.append(centroid.y)
                unit_data.append(pop)
                
                # Append GEOID and district index to file
                unit_districts.append(int(key[config.id_0:]))
                unit_districts.append(id_to_district[key])
                district_pops[id_to_district[key]] += pop
                
                tot_pop += pop
                num_valid_units = num_valid_units + 1

    call('rm ' + config.temp_folder + state + '*', shell=True)

    # Write unit data to file
    output_file = open(config.temp_folder + state + config.unit_data_suffix, 'wb')
    float_array = array('d', unit_data)
    float_array.tofile(output_file)
    output_file.close()

    # Write containing district data to file
    output_file = open(config.temp_folder + state + config.unit_district_suffix, 'wb')
    float_array = array('d', unit_districts)
    float_array.tofile(output_file)
    output_file.close()

    return [num_valid_units, num_dists, tot_pop]



# Main function without re-loading data
def run_with_data(state,k,max_iter,initial_state,ms_param,stopCrit,lb_frac,temp,annealing,verbose,driving_distance,mode,alt_map,p):

    print state + ' has ' + str(p[0]) + ' tracts and ' + str(p[1]) + ' districts.'
    print '\nRunning ' + state + '...'
    lowerBound = round(lb_frac*p[2]/p[1])

    if(verbose):
        print '\nMinimum district population set to ' + str(lowerBound) + '.'
    
    
    # Remove old files iteration files
    call('rm ' + config.temp_folder + state + '*step*', shell=True)


    # Call main algorithm
    command = ['./district_mbo', config.temp_folder + state + config.unit_data_suffix, str(p[0]),str(p[1]),str(k),str(max_iter),str(initial_state),str(ms_param),str(stopCrit),str(lowerBound), config.temp_folder + state + config.unit_district_suffix,str(temp),str(annealing),str(verbose),str(driving_distance)]
    call(command)

    # Count number of data files saved
    num_iter = len(glob.glob(config.temp_folder + state + '*step*'))
    print str(num_iter) + ' iterations completed.'
    
    # Remove old images
    call('rm output/flow*', shell=True)

    if mode != config.MODE_NONE:
        visualize_maps.make_pics(state,data_for_state(state),p[1],num_iter,mode)
    
    return num_iter


# Main function
def run(state,k,max_iter,initial_state,ms_param,stopCrit,lb_frac,temp,annealing,verbose,driving_distance,mode,alt_map):

    p = read(state,alt_map)
    num_iter = run_with_data(state,k,max_iter,initial_state,ms_param,stopCrit,lb_frac,temp,annealing,verbose,driving_distance,mode,alt_map,p)

    return num_iter


# If running straight from the command line:
if __name__ == '__main__':

    if(len(sys.argv) == 1):
        state = 'VA' # Run VA by default
    else:
        state = sys.argv[1]


    ### Parameters for auction dynamics algorithm ###
    k = 150                                              # Number of nearest neighbors
    max_iter = 300                                       # Maximum number of iterations to run
    initial_state = config.INIT_CURRENT                  # INIT_RAND: start from random, INIT_CURRENT: start with 2016 districts, INIT_ALT: alternative map
    ms_param = 1                                         # Mumford-Shah parameter -- higher values means less likely splits
    stopCrit = 0.00                                      # Stopping criteria
    lb_frac = 0.985                                      # Lower bound on the population differences between
    temp =  .1                                           # Temperature
    annealing = 0.985                                    # Multiplicative annealing term
    verbose = 0                                          # 1: print during ms_mbo call, 0: suppress output
    driving_distance = 0                                 # Make 1 to use driving distance for states where that data is available
    mode = config.MODE_BEGIN_END                         # Generate visualization with mode : 0 (first and last), 1 (log sampling), 2 (all)
    alt_map = ''                                         # Alternative mapping file

    run(state,k,max_iter,initial_state,ms_param,stopCrit,lb_frac,temp,annealing,verbose,driving_distance,mode,alt_map)

