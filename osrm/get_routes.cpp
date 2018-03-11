#include "osrm/match_parameters.hpp"
#include "osrm/nearest_parameters.hpp"
#include "osrm/route_parameters.hpp"
#include "osrm/table_parameters.hpp"
#include "osrm/trip_parameters.hpp"
#include "osrm/coordinate.hpp"
#include "osrm/engine_config.hpp"
#include "osrm/json_container.hpp"

#include "osrm/osrm.hpp"
#include "osrm/status.hpp"

#include <exception>
#include <iostream>
#include <string>
#include <utility>

#include <chrono>
#include <cstdlib>

/* Usage: ./get_routes [map.osrm] [block_data] [block_count] */

typedef struct{
    double id; /* GEOID for block */
    double latitude;
    double longitude;
    double pop;
}block_data;

/* All-purpose function to write data to file */
static void overwrite_data_to_file(void *data, const char *fileName, size_t numBytes){
    FILE *f=fopen(fileName,"w");
    fwrite(data,1,numBytes,f);
    fclose(f);
}

/* All-purpose function to read data from file */
static void *get_file_data(const char *fileName){
    
    FILE *f = fopen(fileName,"r");
    size_t fsize;
    fseek(f,0L,SEEK_END);
    fsize = ftell(f);
    fseek(f,0L,SEEK_SET);
    void *data = malloc(fsize);
    fread(data,1,fsize,f);
    fclose(f);
    return data;
}


int main(int argc, const char *argv[])
{
    
    if (argc < 4)
    {
        std::cerr << "Usage: " << argv[0] << " [map.osrm] [block_data] [block_count]\n";
        return EXIT_FAILURE;
    }
    
    /* Read in the block data */
    block_data *stateData = (block_data*)get_file_data(argv[2]);
    int pcount = atoi(argv[3]);
    
    /* Initialize output holder */
    double *output = (double *)calloc(pcount*pcount,sizeof(double));
    
    /* Adapted from OSRM sample code */
    using namespace osrm;
    
    EngineConfig config;
    
    config.storage_config = {argv[1]};
    config.use_shared_memory = false;
    
    
    /* From OSRM sample code:
    // We support two routing speed up techniques:
    // - Contraction Hierarchies (CH): requires extract+contract pre-processing
    // - Multi-Level Dijkstra (MLD): requires extract+partition+customize pre-processing */
    
    // config.algorithm = EngineConfig::Algorithm::CH; // 
    config.algorithm = EngineConfig::Algorithm::MLD;
    
    const OSRM osrm{config};
    
    TableParameters tableParameters;
    
    for(int i = 0; i < pcount; i++){
        tableParameters.coordinates.push_back({util::FloatLongitude{stateData[i].longitude}, util::FloatLatitude{stateData[i].latitude}});
    }
    
    for(int source = 0; source < pcount; source++){
        std::cout << (double)source/(double)pcount << " % complete" << std::endl;
        tableParameters.sources.clear();

        tableParameters.sources.push_back(source); // Single source
        //tableParameters.destinations.push_back(dest); // Leaving blank does 1 -> All
        
        json::Object tableResult;
        
        const auto rc = osrm.Table(tableParameters, tableResult);
        const auto code = tableResult.values.at("code").get<json::String>().value;
        const auto &durations_array = tableResult.values.at("durations").get<json::Array>().values;
        
        const auto durations_matrix = durations_array[0].get<json::Array>().values;
        
        for(int j = 0; j < durations_matrix.size(); j++){
            output[source*pcount + j] = durations_matrix.at(j).get<json::Number>().value;
        }
        
    }
    
    /* Save output to file */
    overwrite_data_to_file(output, "route_dist", pcount*pcount*sizeof(double));
    free(output);
    
    return 0;
    
}
