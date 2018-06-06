#define k_type double
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <vl/kdtree.h>
#include "maj_pop_heap.h"

#define SCALING_FACTOR 0.25
#define EPSILON_0 1
#define EPSILON_MIN 1e-8
#define TEMP_THRESH 0.0001   /* Temperature threshold for stopping */
#define NUM_COMPARISONS 10000

typedef struct{
    double *weights; /* Weight between two census units */
    double *degree; /* diagonal matrix where the diagonal entries are the row sum of the weight matrix */
    int *counts; /* counts[i] is the index where neighbors of i begin */
    int *neighbors; /* neighbors[counts[i] + j] is the identity of i's jth neighbor; weights follows same convention */
}weightedGraph;

typedef struct{
    double id; /* GEOID for district */
    double latitude;
    double longitude;
    double pop;
}unitGeoData;

typedef struct{
    double id; /* GEOID for census units */
    double district; /* Initial district identifier: Used when initializing from pre-set map */
}unitDistrictData;


typedef struct{ /* For sorting doubles while keeping track of indices*/
    double dist;
    int index;
}indexedDouble;

typedef struct{
    weightedGraph g; /* Graph for MBO */
    double *convolutionCoefficients; /* convolutionCoefficients[lcount*(i)) + j] : the district j coefficient for census unit (unit, tract, etc.) i  */
    int *labelMasses; /* labelMasses[lcount*(i)) + j] : the population of unit i assigned to district j */
    double *districtCentroids; /* Centroids of each district */
    double *unitCentroids; /* Centroids of each unit */
    int *districtPopulation; /* Populations of each district */
    int *unitPopulation; /* Populations of each unit */
    double msParameter; /* Mumford-Shah parameter: importance of distance to district centroid vs. district boundary length */
    double stoppingCriteria; /* Halt iterations if changed/totalPop <= stoppingCriteria */
    int pcount; /* Number of units */
    int lcount; /* Number of districts */
    int maxIters; /* Hard upper limit of iterations to run */
    double temperature; /* Randomness term; introduced in convolution step */
    double annealing; /* Multiplicative rate of annealing of temperature */
    int initializeRandom; /* 1 if want random initialization */
}mboStruct;


/* Initialize holders for data, verbose flag */
unitGeoData *stateData;
unitDistrictData* unitDistricts;
char *fileName;
int verbose;
double* drivingDistance;
int useDrivingDistance;
double* energyAtStep;


/* All-purpose function to write data to file */
static void writeDataToFile(void *data, const char *fileName, size_t numBytes){
    FILE *f=fopen(fileName,"w");
    fwrite(data,1,numBytes,f);
    fclose(f);
}

/* All-purpose function to read data from file */
static void *readDataFromFile(const char *fileName){
    
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

/* Save state to file */
void saveState(mboStruct mbos, int iteration){
    
    int pcount = mbos.pcount;
    int lcount = mbos.lcount;
    int *output = calloc(pcount*(lcount + 1),sizeof(int)); /* +1 to include the district ID */
    
    /* Zero-out popCount and compute by summing over all units */
    int popCount[lcount];
    memset(popCount,0,lcount*sizeof(int));
    
    int splitCount = 0;

    for(int i = 0; i < pcount; i++){
        output[(lcount + 1)*i] = stateData[i].id;
        
        int massYet = 0;
        int didSplit = 0;
        for(int d = 0; d < lcount; d++){
            
            output[(lcount + 1)*i + d + 1] = mbos.labelMasses[lcount*i + d];
            popCount[d] += mbos.labelMasses[lcount*i + d];
            
            
            if(mbos.labelMasses[lcount*i + d] > 0 && massYet == 1){
                didSplit = 1;
            }
            
            if(mbos.labelMasses[lcount*i + d] > 0 && massYet == 0){
                massYet = 1;
            }
            
        }
        splitCount += didSplit; /* Track splitting */
    }
    
    for(int d = 0; d < lcount; d++){
        if(verbose){
            printf("Pop %d: %d\n",d,popCount[d]);
        }
    }
    
    if(verbose){
        printf("Fraction of units split: %f \n",((double)splitCount)/((double) pcount));
    }
    
    /* Write "snapshot" of flow to file */
    char buffer[32];

    snprintf(buffer, sizeof(char) * 32, "%s_step_%i", fileName, iteration);
    writeDataToFile(output, buffer, pcount*(lcount + 1)*sizeof(int));
    free(output);
}



/* Compare function for sort that keeps track of index of sorted object */
int cmpIndexedDouble(const void *x, const void *y)
{
    indexedDouble xx = *(indexedDouble*)x, yy = *(indexedDouble*)y;
    if (xx.dist < yy.dist) return -1;
    if (xx.dist > yy.dist) return  1;
    return 0;
}



/* Compare function for sort */
int cmp(const void *x, const void *y)
{
    double xx = *(double*)x, yy = *(double*)y;
    if (xx < yy) return -1;
    if (xx > yy) return  1;
    return 0;
}


/* Build weight matrix from census unit longitude/latitude data */
weightedGraph createSymmetricAdjacencyMatrix(double *weights, int *neighbors, int pcount, int k){
    
    // Use medians to rescale distance
    double *medians = calloc(pcount,sizeof(double));
    
    for(int i = 0; i < pcount; i++){
        medians[i] = weights[i*k + k/2];
    }
    
    /* Build temporary asymmetric graph (tcounts : temporary counts) to store number of counts for each one */
    int *tcounts = calloc(pcount + 1,sizeof(int));
    for(int i = 0; i < pcount; i++){
        for(int j = 0; j < k; j++){
            int index = neighbors[k*i + j];
            tcounts[index]++; // Increment neighbor_count for neighbor at index
            tcounts[i]++; // Increment neighbor_count for i
        }
    }
    
    /* Convert # counts for each district to index of the weight/neighbor array for that district */
    for(int i = 0; i < pcount; i++){
        tcounts[i + 1] += tcounts[i];
    }
    
    int num = tcounts[pcount]; /* Total number of neighbors - initially 2*k*pcount, changes when we symmetrize, tcounts[pcount] = tcounts[pcount - 1] here */
    double *tempWeights = calloc(num,sizeof(double));
    int *tempIndices = calloc(num,sizeof(int));
    
    /* Build temporary graph */
    for(int i = 0; i < pcount; i++){
        for(int j = 0; j < k; j++){
            
            int index = neighbors[k*i + j]; /* Identify of i's jth neighbor */
            double m1 = medians[i]; /* Median distance for unit i, index */
            double m2 = medians[index];
            
            double factor = sqrt(m1*m2); /* Weight for i, index depends on both medians so points with high medians aren't penalized */
            
            tempIndices[tcounts[i] - 1] = index;
            tempWeights[tcounts[i] - 1] = weights[i*k + j]/factor;
            
            tcounts[i]--;
            
            tempIndices[tcounts[index] - 1] = i;
            tempWeights[tcounts[index] - 1] = weights[i*k + j]/factor;
            
            tcounts[index]--;
            
        }
    }
    
    
    /* tempIndices includes repeats: if i is a neighbor of j and j is a neighbor of i, j is counted twice in i's chunk of tempIndices
    The rest of the code in this function is to remove duplicates, resulting in a non-double counted symmetric graph */
    
    int *counts = calloc(pcount + 1,sizeof(int));
    int *seen = calloc(pcount,sizeof(int));
    
    
   /* Set-up counts to remove duplicates */
    for(int i = 0; i < pcount; i++){
        for(int j = tcounts[i]; j < tcounts[i + 1]; j++){
            int index = tempIndices[j];
            if(seen[index] != i + 1){ // +1 is to avoid mistaking calloc'd 0 for an index.
                seen[index] = i + 1;
                counts[i]++;
            }
        }
    }
    
    /* Shift counts to indices */
    for(int i = 0; i < pcount; i++){
        counts[i + 1] += counts[i];
    }
    
    /* Allocate correctly sized holders post-duplicate removal */
    double *gWeights = calloc(counts[pcount],sizeof(double));
    int *indices = calloc(counts[pcount],sizeof(int));
    
    /* Reset seen to zero */
    memset(seen,0,pcount*sizeof(int));
    
    /* Remove duplicates, storing unique values in appropriately sized holders */
    for(int i = 0; i < pcount; i++){
        for(int j = tcounts[i]; j < tcounts[i + 1]; j++){
            int index = tempIndices[j];
            if(seen[index] != i + 1){
                seen[index] = i + 1;
                indices[counts[i] - 1] = index;
                gWeights[counts[i] - 1] = tempWeights[j];
                counts[i]--;
            }
        }
    }
    
    /* Create the degree matrix - a diagonal matrix of row sums for each census unit */
    double *degree = calloc(pcount,sizeof(double));
    
    for(int i = 0; i < pcount;i++){
        for(int j = counts[i]; j < counts[i+1]; j++){
            degree[i] += gWeights[j];
        }
    }
    
    /* Store in graph struct */
    weightedGraph g;
    g.neighbors = indices;
    g.weights = gWeights;
    g.degree = degree;
    g.counts = counts;

    free(tcounts);
    free(medians);
    free(tempIndices);
    free(tempWeights);
    free(seen);
    return g;
    
}

/* Kernel for weights */
double kernel(double x){
    return exp(-x*x);
}

/* Calls createSymmetricAdjacencyMatrix then kernelizes the weights */
weightedGraph createSymmetricMatrix(double (* kernel)(double), double *weights, int *neighbors, int pcount, int k){
    weightedGraph g = createSymmetricAdjacencyMatrix(weights, neighbors, pcount,  k);
    
    /* Zero out degree matrix */
    memset(g.degree,0,pcount*sizeof(double));
    
    for(int i = 0; i < pcount; i++){
        for(int j = g.counts[i]; j < g.counts[i + 1]; j++){
            g.weights[j] = kernel(g.weights[j]);
            g.degree[i] += g.weights[j];
        }
    }
    
    return g;
}


/* Normalizing is an important step that improves speed and performance
 W_new = D^(-1/2) W_old D^(-1/2), D is the diagonal degree matrix where D(i,i) : the sum of all of i's weights */
void normalizeMatrix(weightedGraph g, int pcount){
   
    double *sums = calloc(pcount,sizeof(double));
    double *weights = g.weights;
    int *indices = g.neighbors;
    int *counts = g.counts;
    
    for(int i = 0; i < pcount; i++){
        for(int j = counts[i]; j < counts[i + 1]; j++){
            sums[i] += weights[j];
        }
    }
    
    for(int i = 0; i < pcount; i++){
        for(int j = counts[i]; j < counts[i + 1]; j++){
            int index = indices[j];
            weights[j]/=sqrt(sums[i]*sums[index]);
        }
    }
    
    free(sums);
    
    /* Zero out and recalculate degree matrix */
    memset(g.degree,0,pcount*sizeof(double));
    
    for(int i = 0; i < pcount; i++){
        for(int j = g.counts[i]; j < g.counts[i + 1]; j++){
            g.degree[i] += g.weights[j];
        }
    }
}


/* Build transpose of W */
weightedGraph createDualWeightedGraph(mboStruct mbos){
    
    int tot = 0;
    int pcount = mbos.pcount;
    
    weightedGraph g=mbos.g;
    
    int *counts = calloc(pcount + 1,sizeof(int)); /* Holder for number of neighbors for each point */
    int *nncounts = mbos.g.counts;
    
    
    for(int i = 0; i < pcount; i++){
        for(int j = nncounts[i]; j < nncounts[i + 1]; j++){
            int neighbor = mbos.g.neighbors[j];
            counts[neighbor]++; /* Count the number of non-zero entries in each column of W; i.e. the number of neighbors for each point */
            tot++;
        }
    }
    
    /* Convert # counts to index */
    for(int i = 0; i < pcount; i++){
        counts[i + 1] += counts[i];
    }
    
    counts[pcount] = tot;
    
    int *neighbors = calloc(tot,sizeof(int));
    double *weights = calloc(tot,sizeof(double));
    
    for(int i = 0; i < pcount; i++){
        for(int j = nncounts[i]; j < nncounts[i + 1]; j++){
            int neighbor = mbos.g.neighbors[j]; /* Get i's jth neighbor */
            double str = mbos.g.weights[j]; /* Get weight between i and jth neighbor*/
            neighbors[counts[neighbor] - 1] = i; /* Set i to be a neighbor of its jth neighbor in neighbors */
            weights[counts[neighbor] - 1] = str; /*  Set weight between i and its jth neighbor  */
            counts[neighbor]--;
        }
    }
    
    
    /* Zero out and recalculate degree matrix */
    /* Dual degree matrix is the degree matrix of W^TW */
    
    double *degree = calloc(pcount,sizeof(double));
    for(int i = 0; i < pcount; i++){
        for(int j = counts[i]; j < counts[i + 1]; j++){
            int neighbor = g.neighbors[j];
            degree[i] += weights[j]*g.degree[neighbor];
        }
    }
    
    /* Create struct */
    weightedGraph dual;
    dual.neighbors = neighbors;
    dual.weights = weights;
    dual.degree = degree;
    dual.counts = counts;
    
    return dual;
    
}

/* Clean up graph */
void destroyWeightedGraph(weightedGraph g){
    free(g.neighbors);
    free(g.weights);
    free(g.counts);
    free(g.degree);
}


/* Calculate the coefficients for the convolution */
void calculateCoefficients(mboStruct mbos, weightedGraph dual){
    
    int pcount = mbos.pcount;
    int lcount = mbos.lcount;
    
    weightedGraph g = mbos.g;
    int *labelMasses = mbos.labelMasses;
    int *unitPopulation = mbos.unitPopulation;
   
    /* We're carrying out W^TWX: temp is WX and each row of X corresponds to a point; each column to a label */
    /* W is pcount x pcount; X is pcount by lcount */
    double *temp = calloc(pcount*lcount,sizeof(double));
    
    for(int i = 0; i < pcount; i++){
        for(int j = g.counts[i]; j < g.counts[i + 1]; j++){
            int index = g.neighbors[j];
            
            for(int d = 0; d < lcount; d++){
                
                /* Convert labeling to probability */
                double simplex = labelMasses[lcount*index + d]/(unitPopulation[index]*1.0);
                temp[lcount*i + d] += g.weights[j]*simplex;
                
            }
        }
    }
    
   

    /* Get previous values and overwrite to zero */
    double *convolutionCoefficients = mbos.convolutionCoefficients;
    memset(convolutionCoefficients, 0, pcount*lcount*sizeof(double));
    
    
    /* Normalizing factor for the heat content energy  */
    double kernelWidth = sqrt(mbos.g.counts[pcount]/(pcount*1.0));

    /* Replace with rest of matrix vector multiplication W^T*(WX) */
    for(int i = 0; i < pcount; i++){
        for(int j = dual.counts[i]; j < dual.counts[i + 1]; j++){
            int index = dual.neighbors[j];
            for(int d = 0; d < lcount; d++){
                convolutionCoefficients[lcount*i + d] += temp[lcount*index + d]*dual.weights[j]/kernelWidth;
            }
        }
    }
    
    free(temp);
    
    /* Mumford-Shah step: Add the centroid distance penalty term to each entry of W^TWDX*/
    double *districtCentroids = mbos.districtCentroids;
    double *unitCentroids = mbos.unitCentroids;
    double msParameter = mbos.msParameter;

    for(int i = 0; i < pcount; i++){
        for(int d = 0; d < lcount; d++){
            /* Compute squared distance (latitude, longitude) between districtCentroids and individual unit centroids and subtract from coefficients */
            convolutionCoefficients[lcount*i + d] -= msParameter*(pow(districtCentroids[2*d] - unitCentroids[2*i],2) + pow(districtCentroids[2*d + 1] - unitCentroids[2*i + 1],2));
        }
    }
}

/* Calculate heat content energy from the convolution values evaluated **at the current configuration** */
/* Sum_i of integral of (1-u_i) times (K * u_i) */
double calculateEnergy(mboStruct mbos, weightedGraph dual){
    
    int pcount = mbos.pcount;
    int lcount = mbos.lcount;
    int *labelMasses = mbos.labelMasses;
    int *unitPopulation = mbos.unitPopulation;
    
    double *degree = dual.degree;
    double kernelWidth = sqrt(mbos.g.counts[pcount]/(pcount*1.0));
    
    double energy = 0;
    
    
    double *convolutionCoefficients = mbos.convolutionCoefficients;
    
    /* Energy calculation step. Heat content energy = (sum over x, y, i) W(x,y)*u_i(x)*(1 - u_i(y))
     Auction uses the negative of the gradient.  Need to compensate for this when calculating the energy*/
    
    for(int i = 0; i < pcount; i++){ // Sum over x
        for(int d = 0; d < lcount; d++){ // Sum over i
            double simplex = labelMasses[lcount*i + d]/(unitPopulation[i]*1.0); // u_i(x)
            double degreeTerm = degree[i]/kernelWidth; // Accounts for the 1 in (1 - u_i(y))
            energy += simplex*(degreeTerm - convolutionCoefficients[lcount*i + d]); // Convolution coefficients contain (sum over y) W(x,y)*u_i(y)
        }
    }
    
    /* Normalize by number of census units in state */
    energy/=pcount;
    
    return energy;
}



/* Update information after an auction step */
void updateDistrictInformation(mboStruct mbos){
    
    int pcount = mbos.pcount;
    int lcount = mbos.lcount;
    int *labelMasses = mbos.labelMasses;
    double *unitCentroids = mbos.unitCentroids;
    int *unitPopulation = mbos.unitPopulation;
    int *districtPopulation = mbos.districtPopulation;
    
    /* Clear district centroids + populations */
    double *districtCentroids = mbos.districtCentroids;
    memset(districtCentroids, 0, 2*lcount*sizeof(double));
    memset(districtPopulation, 0, lcount*sizeof(int));
    
    /* Recompute district centroids + populations by adding up over all census units */
    double unitCount[lcount];
    memset(unitCount,0,lcount*sizeof(double));
    
    /* Weight centroids by amount of units in each district  */
    for(int i = 0; i < pcount; i++){
        double bx = unitCentroids[2*i];
        double by = unitCentroids[2*i + 1];
        for(int d = 0; d < lcount; d++){
            districtCentroids[2*d] += bx*labelMasses[lcount*i + d];
            districtCentroids[2*d + 1] += by*labelMasses[lcount*i + d];
            districtPopulation[d] += labelMasses[lcount*i + d];
            unitCount[d] += labelMasses[lcount*i + d]/(unitPopulation[i]*1.0);
        }
    }
    
    /* Weight centroids by population */
    for(int d = 0; d < lcount; d++){
        districtCentroids[2*d]/=districtPopulation[d];
        districtCentroids[2*d + 1]/=districtPopulation[d];
    }
}


/* Reverse auction phase */
void reverseAuctionPhase(mboStruct mbos, pop_heap *heapHolder, double *incentives, double epsilon, int lowerBound){
   
    double *convolutionCoefficients = mbos.convolutionCoefficients;
    int *unitPopulation = mbos.unitPopulation;
    int *districtPopulation = mbos.districtPopulation;
    int *labelMasses = mbos.labelMasses;
    int pcount = mbos.pcount;
    int lcount = mbos.lcount;

    /* Loop over all districts */
    for(int d = 0; d < lcount; d++){
        
        // If district has insufficient population...
        if(districtPopulation[d] < lowerBound){
            
            /* Get the gap between the lower bound and the district's current population and reset heap */
            int gap = lowerBound - districtPopulation[d];
            pop_heap_reset(&heapHolder[d]);
            
            /* Loop over all points*/
            for(int i = 0; i < pcount; i++){

                /* If we don't have the full unit assigned to this district already... */
                if(labelMasses[lcount*i + d] < unitPopulation[d]){
                    
                    /* Loop over all the districts again */
                    for(int k = 0; k < lcount; k++){
                        
                        /* If we find another one with non-zero mass... */
                        if(k != d && labelMasses[lcount*i + k] > 0){
                            
                            /* Compute the gap between how badly census unit i wants to be in district k and district d */
                            double cost = (convolutionCoefficients[lcount*i + k] + incentives[k]) - (convolutionCoefficients[lcount*i + d] + incentives[d]);
                            
                            /* Get the mass (population) at this entry. */
                            int mass = labelMasses[lcount*i+k];
                            int lbl = k;
                            
                            double key = -cost;
                            
                            /* District makes "decision" on whether or nor the census unit is a good enough candidate to justify paying its price */
                            /* Before reaching lower bound, candidates are always accepted by the district; after "sateity", this step looks for "good deals" */
                            
                            pop_heap_consider_candidate(&heapHolder[d], i, key, mass, lbl, gap);
                            
                        }
                    }
                }
                
            }
            
            double cost = -FLT_MAX;
            
            /* While each district still has "cheapest units" held, transfer the mass of those units to that district */
            while(heapHolder[d].count > 0){
                
                pop_node node = extract_min(&heapHolder[d]);
                
                int index = node.originalIndex;
                int mass = node.mass;
                int lbl = node.lbl;
                
                if(-node.key > cost){
                    cost = -node.key;
                }
                
                /* Move the mass from lbl to d */
                labelMasses[lcount*index + lbl] -= mass;
                labelMasses[lcount*index + d] += mass;
                
                districtPopulation[lbl] -= mass;
                districtPopulation[d] += mass;
                
            }

            /* Update incentives that needed to be offered to get the census units */
            incentives[d] += fmax(cost + epsilon,0);
        
        }
    }
}

/* Auction set-up */
void prepareAuctionPhase(mboStruct mbos, double *incentives){
    
    double *convolutionCoefficients = mbos.convolutionCoefficients;
    int *unitPopulation = mbos.unitPopulation;
    int *districtPopulation = mbos.districtPopulation;
    int *labelMasses = mbos.labelMasses;
    int pcount = mbos.pcount;
    int lcount = mbos.lcount;

    /* Zero out population, masses */
    memset(labelMasses, 0, pcount*lcount*sizeof(int));
    memset(districtPopulation, 0, lcount*sizeof(int));
    
    
    /* Assign every census unit to the district they want most based on current incentives */
    for(int i = 0; i < pcount; i++){
        double max = -FLT_MAX;
        int ml = 0; /* Winning label */
        for(int d = 0; d < lcount; d++){
            if(convolutionCoefficients[lcount*i + d] + incentives[d] > max){
                max = convolutionCoefficients[lcount*i + d] + incentives[d];
                ml = d;
            }
        }
        
        labelMasses[lcount*i + ml] = unitPopulation[i];
        districtPopulation[ml] += unitPopulation[i];
    }
}


/* Check if each district satisfies lower population bound */
int feasibilityTest(int *districtPopulation, int lowerBound, int lcount){

    int feasible = 1;
    for(int d = 0; d < lcount; d++){
        if(districtPopulation[d] < lowerBound){
            feasible = 0;
        }
    }
    
    return feasible;
}

/* Creates two Gaussian random variables of mean mu and std_dev sigma */
void boxMuellerTransform(double *z1, double *z2, double mu, double sigma){
    double r1 = rand()/(RAND_MAX*1.0);
    double r2 = rand()/(RAND_MAX*1.0);
    double lsq = sqrt(-2*log(r1));
    *z1 = mu + sigma*lsq*cos(2*M_PI*r2);
    *z2 = mu + sigma*lsq*sin(2*M_PI*r2);
}

/* Add noise to convolution coefficients */
void addTemperature(mboStruct mbos){
    double temperature = mbos.temperature;
   
    if(verbose){
        printf("Temperature: %f\n", temperature);
    }
    
    double noise1, noise2; // Noise has to be generated in pairs

    int pcount = mbos.pcount;
    int lcount = mbos.lcount;
    
    double *convolutionCoefficients = mbos.convolutionCoefficients;
    
    for(int i = 0; i < pcount; i++){
        for(int d = 0; d < lcount; d++){
    
            boxMuellerTransform(&noise1,&noise2,0,temperature);
            convolutionCoefficients[lcount*i + d] += noise1;
            
        }
    }
}


/* Runs auction */
void districtAuction(mboStruct mbos, int lowerBound){
    
    /* At every step we lower the tolerance for error by multiplying by SCALING_FACTOR */
    int *districtPopulation = mbos.districtPopulation;
    double epsilon = EPSILON_0;
    int lcount = mbos.lcount;

    
    /* Initialize heap */
    pop_heap heapHolder[lcount];
    
    for(int d = 0; d < lcount; d++){
        heapHolder[d] = pop_heap_create_empty_heap(lowerBound);
    }
    
    /* Initialize incentives */
    double incentives[lcount];
    memset(incentives, 0, lcount*sizeof(double));
    
    
    /* Get solution with error at most EPSILON_MIN */
    while(epsilon > EPSILON_MIN){

        int feasible = 0;
        
        prepareAuctionPhase(mbos, incentives);
        feasible = feasibilityTest(districtPopulation, lowerBound, lcount);
        
        /* Run reverse auction until feasible solution found */
        while(!feasible){
            reverseAuctionPhase(mbos, heapHolder, incentives, epsilon, lowerBound);
            feasible = feasibilityTest(districtPopulation, lowerBound, lcount);
        }
        
        /* Feasible solution found! Try again with lower error tolerance. */
        epsilon *= SCALING_FACTOR;
    }
    
    /* Clean-up */
    for(int d = 0; d < lcount; d++){
        pop_heap_destroy_heap(&heapHolder[d]);
    }
}


/* Tracks change after one full MBO iteration; used to determine when to stop */
int massChanged(int *labelMasses, int *oldLabelMasses, int pcount, int lcount){
    
    int change = 0;
    for(int i = 0; i < pcount*lcount; i++){
        change += abs(labelMasses[i] - oldLabelMasses[i]);
    }
    return change;
}


/* Main function: Runs Auction dynamics + Centroid distance penalty (Mumford-Shah) */
void runMBO(mboStruct mbos, int lowerBound){
    
    int pcount = mbos.pcount;
    int lcount = mbos.lcount;
    int maxIters = mbos.maxIters;
    
    /* Create dual */
    weightedGraph dual = createDualWeightedGraph(mbos);
    
    int *labelMasses = mbos.labelMasses;
    int *oldLabelMasses = calloc(pcount*lcount,sizeof(int));

    memcpy(oldLabelMasses,labelMasses,pcount*lcount*sizeof(int));
    int *unitPopulation = mbos.unitPopulation;
    
    double stoppingCriteria = mbos.stoppingCriteria;
    
    int totalPop = 0;
    for(int i = 0; i < pcount; i++){
        totalPop += unitPopulation[i];
    }
    
    /* 'best' stores configuration with lowest energy; only relevant for temperature */
    double best = FLT_MAX;
    int *bestLabels = calloc(pcount*lcount,sizeof(int));
    
    if(verbose){
        printf("Starting iterations...\n");
    }
    
    updateDistrictInformation(mbos);

    
    for(int i = 0; i < maxIters; i++){
        
        saveState(mbos,i);
        
        /* Calculate coefficients */
        if(verbose){
            printf("Calculating coefficients...\n");
        }
        
        calculateCoefficients(mbos, dual);
        
        /* For temperature: compute energy to store lowest configuration */
        /* NOTE: Very important that energy is calculated AFTER coefficients are calculated */
        double energy = calculateEnergy(mbos, dual);
        
        if(verbose){
            printf("Energy at iteration %d is: %f\n", i, energy);
        }
        
        energyAtStep[i] = energy;
        
        if(energy < best){
            best = energy;
            memcpy(bestLabels, mbos.labelMasses, pcount*lcount*sizeof(int));
        }
        
        /* Annealing step + temperature add step */
        mbos.temperature = mbos.temperature*mbos.annealing;
        addTemperature(mbos);
        
        if(verbose){
            printf("Auction step...\n");
        }
        
        /* Run auction */
        districtAuction(mbos, lowerBound);
        
        /* Update post-auction*/
        updateDistrictInformation(mbos);
        
        /* See if stopping criteria has been met */
        int changed = massChanged(labelMasses, oldLabelMasses, pcount, lcount);
        if(changed/(totalPop*1.0) <= stoppingCriteria && mbos.temperature < TEMP_THRESH){
            break;
        }
        
        if(verbose){
            printf("Auction complete.\n");
        }
        
        memcpy(oldLabelMasses,labelMasses,pcount*lcount*sizeof(int));
 
    }
    
    
    
    
    /* Replace labelMasses with best configuration (in case of temperature) */
    memcpy(labelMasses,bestLabels,pcount*lcount*sizeof(int));
    
    updateDistrictInformation(mbos);
    
    free(bestLabels);
    free(oldLabelMasses);
    destroyWeightedGraph(dual);
}

/* Wrapper for vlfeat functions */
void computeNearestNeighbors(double *unitCentroids, double *distances, int *neighbors, int pcount, int k){
    
    int numTrees = 1;
    int fullDim = 2;
   
    vl_size numN = k;
    vl_size tot = pcount;
    VlKDForest *forest = vl_kdforest_new(VL_TYPE_DOUBLE,fullDim,numTrees,VlDistanceL2);
    vl_kdforest_build(forest,tot,unitCentroids);
    vl_kdforest_set_max_num_comparisons(forest,NUM_COMPARISONS);
    
    vl_uint32 *indices=(vl_uint32 *) neighbors;
    vl_kdforest_query_with_array(forest,indices, numN, tot,distances,unitCentroids);
    
}


/* Sorts driving distances each time the code is called */
/* Safer alternative that uses driving distance to compute all neighbors */
/* Later change to make: Make driving_distance.h -- or weight matrix construction as its own plug-in <----  */
void computeNearestNeighborsDriving(double *drivingDistance, double *distances, int *neighbors, int pcount, int k){
    
    indexedDouble *sortedDrivingDistances = calloc(pcount*pcount,sizeof(indexedDouble));
    
    for(int i = 0; i < pcount; i++){
        for(int j = 0; j < pcount; j++){
            
            sortedDrivingDistances[pcount*i + j].dist = drivingDistance[pcount*i + j];
            sortedDrivingDistances[pcount*i + j].index = j;
        }
        
        qsort(&sortedDrivingDistances[pcount*i], pcount, sizeof(indexedDouble), cmpIndexedDouble);
    }
    
    for(int i = 0; i < pcount; i++){
        for(int j = 0; j < k; j++){
            distances[k*i + j] = sortedDrivingDistances[pcount*i + j].dist;
            neighbors[k*i + j] = sortedDrivingDistances[pcount*i +j].index;
        }
    }
    
    free(sortedDrivingDistances);
}


/* Create mbo struct from parameters and loaded files   */
mboStruct createMBO(unitGeoData *stateData, double msParameter, double stoppingCriteria, int pcount, int lcount, int k, int maxIters, unitDistrictData *unitData, double temp, double annealing, int initializeRandom){
        mboStruct mbos;
    
    mbos.convolutionCoefficients = calloc(pcount*lcount,sizeof(double));
    mbos.labelMasses = calloc(pcount*lcount,sizeof(int));
    mbos.districtCentroids = calloc(2*lcount,sizeof(double));
    mbos.unitCentroids = calloc(2*pcount,sizeof(double));
    mbos.unitPopulation = calloc(pcount,sizeof(int));
    mbos.districtPopulation = calloc(lcount,sizeof(int));
    mbos.temperature = temp;
    mbos.annealing = annealing;
    mbos.initializeRandom = initializeRandom;
    
    /* Scale by latitude and longtiude so max dimension is 1*/
    double maxLat = -FLT_MAX;
    double minLat = FLT_MAX;
    double maxLong = -FLT_MAX;
    double minLong = FLT_MAX;
    
    for(int i = 0; i < pcount; i++){
        if(stateData[i].latitude > maxLat){
            maxLat = stateData[i].latitude;
        }
        if(stateData[i].latitude < minLat){
            minLat = stateData[i].latitude;
        }
        if(stateData[i].longitude > maxLong){
            maxLong = stateData[i].longitude;
        }
        if(stateData[i].longitude < minLong){
            minLong = stateData[i].longitude;
        }
    }
    
    /* Get scalar to divide by*/
    double scaleDimension = (((maxLat - minLat) > (maxLong - minLong)) ? (maxLat - minLat) : (maxLong - minLong));
    
    /* Initialize labelMasses to zero */
    for(int d = 0; d < lcount; d++){
        for(int i = 0; i < pcount; i++){
            mbos.labelMasses[lcount*i + d] = 0;
        }
    }
    
    for(int i = 0; i < pcount; i++){
        mbos.unitPopulation[i] = stateData[i].pop + 1;
        mbos.unitCentroids[2*i] = (stateData[i].latitude - minLat)/scaleDimension;
        mbos.unitCentroids[2*i + 1] = (stateData[i].longitude - minLong)/scaleDimension;
        
        /* Initialize with current ID */
        int lbl = (int)unitData[i].district;
        
        if(mbos.initializeRandom == 1){
            lbl = rand() % lcount;
        }
        
        mbos.labelMasses[lcount*i + lbl] = mbos.unitPopulation[i];
    }
    
    
    mbos.msParameter = msParameter;
    mbos.stoppingCriteria = stoppingCriteria;
    mbos.pcount = pcount;
    mbos.lcount = lcount;
    mbos.maxIters = maxIters;
    
    
    double *distances = calloc(k*pcount,sizeof(double));
    int *neighbors = calloc(k*pcount,sizeof(int));
    
    
    /* Compute nearest neighbors, create and normalize W*/
    if(useDrivingDistance){
       computeNearestNeighborsDriving(drivingDistance, distances, neighbors, pcount, k);
    }else{
        computeNearestNeighbors(mbos.unitCentroids, distances, neighbors, pcount, k);
    }
    
    
    mbos.g = createSymmetricMatrix(kernel, distances, neighbors, pcount, k);
    
    normalizeMatrix(mbos.g,mbos.pcount); /* Normalizing is an important step that improves speed and performance */

    free(distances);
    free(neighbors);
    
    return mbos;
    
}

/* Clean up */
void destroyMBOStruct(mboStruct mbos){
    
    free(mbos.convolutionCoefficients);
    free(mbos.labelMasses);
    free(mbos.districtCentroids);
    free(mbos.unitCentroids);
    free(mbos.unitPopulation);
    free(mbos.districtPopulation);
    
    destroyWeightedGraph(mbos.g);
    
}




void runColorAuction(double *assignmentCoefficients, int *colorAssignments, int lcount){
    
    double epsilonMin=1e-6;
    double scalingFactor=4;
    double epsilon=1;
    
    double prices[lcount];
    memset(prices,0,lcount*sizeof(double));
    
    
    
    
    int done=0;
    
    while(epsilon>epsilonMin){
        
        
        memset(colorAssignments, -1, lcount*sizeof(int));
        done=0;
        
        
        while(!done){
            
            for(int i=0;i<lcount;i++){
                if(colorAssignments[i]<0){
                    int ml=0;
                    double max=-FLT_MAX;
                    double next=-FLT_MAX;
                    for(int j=0; j<lcount;j++){
                        if(assignmentCoefficients[lcount*i+j]-prices[j]>max){
                            next=max;
                            ml=j;
                            max=assignmentCoefficients[lcount*i+j]-prices[j];
                            
                            
                        }else if(assignmentCoefficients[lcount*i+j]-prices[j]>next){
                            next=assignmentCoefficients[lcount*i+j]-prices[j];
                        }
                    }
                    
                    double bid=max-next+epsilon;
                    prices[ml]+=bid;
                    for(int j=0;j<lcount;j++){
                        if(colorAssignments[j]==ml){
                            colorAssignments[j]=-1;
                        }
                        
                    }
                    colorAssignments[i]=ml;
                   
                    
                    
                }
            }
            
            
            done=1;
            
            for(int i=0;i<lcount;i++){
                if(colorAssignments[i]<0){
                    done=0;
                }
            }
        }
        
        epsilon/=scalingFactor;
        printf("%f\n",epsilon);
        
        
    }
    
   
    
    
}


void assignDistrictColors(double *colors, double *districtCentroids, int *colorAssignments, int lcount){
    
    double weights[lcount*lcount];
    double helper[lcount];
    double median[lcount];
    
   
    for(int i=0;i<lcount;i++){
        colorAssignments[i]=i;
    }
    
    for(int i = 0 ;i < lcount; i++){
        for(int j = 0 ;j < lcount; j++){
            
            weights[lcount*i + j] = pow(districtCentroids[2*i] - districtCentroids[2*j],2) +  pow(districtCentroids[2*i + 1] - districtCentroids[2*j + 1],2);
            helper[j] = weights[lcount*i +j];
            
        }
        
        qsort(helper, lcount, sizeof(double), cmp);
        median[i] = helper[lcount/2];
    }
    
    double sums[lcount];
    memset(sums,0,lcount*sizeof(double));
    
    for(int i = 0 ;i < lcount; i++){
        
        for(int j = 0 ;j < lcount; j++){
            
            double factor = sqrt(median[i]*median[j]);
            weights[lcount*i + j] = kernel(weights[lcount*i + j]/factor);
            sums[i] += weights[lcount*i + j];
            
        }
    }
    
    for(int i = 0 ;i < lcount; i++){
        for(int j = 0 ;j < lcount; j++){
            weights[lcount*i + j] /= sqrt(sums[i]*sums[j]);
            
        }
    }
    
    for(int i=0;i<lcount;i++){
        weights[lcount*i+i]=0;
    }
    
    
    double convolution[3*lcount];
    
    memset(convolution, 0, 3*lcount*sizeof(double));
    
    int iterations=50;
    
    for(int i=0; i<iterations;i++){
        for(int a = 0; a < lcount ; a++){
            for(int b = 0; b < lcount ; b++){
               
               
                
                convolution[3*a] += weights[lcount*a + b]*colors[3*colorAssignments[b]];
                convolution[3*a+1] += weights[lcount*a + b]*colors[3*colorAssignments[b]+1];
                convolution[3*a+2] += weights[lcount*a + b]*colors[3*colorAssignments[b]+2];
                
               
               
            }
            
        }
        
        double assignmentCoefficients[lcount*lcount];
        
        for(int j=0;j< lcount*lcount;j++){
            
            int a=j%lcount;
            int b=j/lcount;
            
            assignmentCoefficients[lcount*b+a]=pow(convolution[3*a]-convolution[3*b],2)+pow(convolution[3*a+1]-convolution[3*b+1],2)+pow(convolution[3*a+2]-convolution[3*b+2],2);
            
            
        }
        
        
        runColorAuction(assignmentCoefficients, colorAssignments, lcount);
        
        
    }
}




/* For flair: compute colors that are dissimilar for near districts for visualizing */
void computeDistrictColors(double *colors, double *districtCentroids, int lcount){
    
    double weights[lcount*lcount];
    double helper[lcount];
    double median[lcount];
    
    for(int i = 0; i < lcount ; i++){
        
        colors[3*i] = rand()/(RAND_MAX*1.0);
        colors[3*i + 1] = rand()/(RAND_MAX*1.0);
        colors[3*i + 2] = rand()/(RAND_MAX*1.0);
        
    }
    
    for(int i = 0; i < lcount; i++){
        for(int j = 0 ;j < lcount; j++){
            
            weights[lcount*i + j] = pow(districtCentroids[2*i] - districtCentroids[2*j],2) +  pow(districtCentroids[2*i + 1] - districtCentroids[2*j + 1],2);
            helper[j] = weights[lcount*i +j];
            
        }
        qsort(helper, lcount, sizeof(double), cmp);
        median[i] = helper[lcount/2];
    }
    
    double sums[lcount];
    memset(sums,0,lcount*sizeof(double));
    
    for(int i = 0 ;i < lcount; i++){
        
        for(int j = 0 ;j < lcount; j++){

            double factor = sqrt(median[i]*median[j]);
            weights[lcount*i + j] = kernel(weights[lcount*i + j]/factor);
            sums[i] += weights[lcount*i + j];
            
        }
    }
    
    for(int i = 0 ;i < lcount; i++){
        for(int j = 0 ;j < lcount; j++){
            weights[lcount*i + j] /= sqrt(sums[i]*sums[j]);
        }
    }
    
    int iterations = 500;

    for(int i = 0; i < iterations; i++){
        
        double gradient[3*lcount];
        memset(gradient,0,3*lcount*sizeof(double));
        
        for(int a = 0; a < lcount ; a++){
            for(int b = 0; b < lcount ; b++){
                
                gradient[3*a] += weights[lcount*a + b]*(colors[3*a] - colors[3*b])/2;
                gradient[3*a + 1] += weights[lcount*a + b]*(colors[3*a + 1] - colors[3*b + 1])/2;
                gradient[3*a + 2] += weights[lcount*a + b]*(colors[3*a + 2] - colors[3*b + 2])/2;
                
                gradient[3*b] += weights[lcount*a + b]*(colors[3*b] - colors[3*a])/2;
                gradient[3*b + 1] += weights[lcount*a + b]*(colors[3*b + 1] - colors[3*a + 1])/2;
                gradient[3*b + 2] += weights[lcount*a + b]*(colors[3*b + 2] - colors[3*a + 2])/2;
                
            }
        }
        
        for(int a = 0; a < lcount ; a++){
            
            colors[3*a] += .5*gradient[3*a];
            colors[3*a] = fmin(fmax(colors[3*a],0), 1);
            
            colors[3*a + 1] += .5*gradient[3*a + 1];
            colors[3*a + 1] = fmin(fmax(colors[3*a + 1],0), 1);
            
            colors[3*a + 2] += .5*gradient[3*a + 2];
            colors[3*a + 2] = fmin(fmax(colors[3*a + 2],0), 1);
            
        }
    }
}



int main(int argc, char *argv[]){
    
    srand(time(NULL));

    fileName = argv[1]; /* Filename of census unit input data file; format is (idnum, lat, long, pop) [doubles] */

    int pcount = atoi(argv[2]); /* Number of census units */
    int lcount = atoi(argv[3]); /* Number of districts */
    int k = atoi(argv[4]); /* Number of neighbors in the weighted graph */
    int maxIters = atoi(argv[5]); /* Max number of iterations to run */
    int initializeRandom = atoi(argv[6]); /* Random initialization if 1, 0 otherwise */
    double msParameter = atof(argv[7]); /* Mumford-Shah Penalty Weight: determines relative importance of centroid distances */
    double stoppingCriteria = atof(argv[8]);  /* Convergence halting parameter */
    int lowerBound = atoi(argv[9]);  /* Min population per district. */
    char *unitDistrictFile = argv[10]; /* Initial districts file for units */
    double temp = atof(argv[11]); /* Temperature */
    double annealing = atof(argv[12]); /* Multiplicative annealing rate */
    verbose = atoi(argv[13]); /* Multiplicative annealing rate */
    useDrivingDistance = atoi(argv[14]); /* 0: use default distance; 1: use driving distance */

    /* Read data */
    stateData = readDataFromFile(fileName); //read in the census unit data
    unitDistricts = readDataFromFile(unitDistrictFile);
    
    if(useDrivingDistance){
        printf("Loading MD driving data...\n");
        drivingDistance = readDataFromFile("data/MD_DrivingDistance");
    }
    
    /* Create holder for energy */
    energyAtStep = calloc(maxIters,sizeof(double));
    

    /* MBO struct created here */
    mboStruct mbos = createMBO(stateData, msParameter, stoppingCriteria, pcount, lcount, k,  maxIters, unitDistricts, temp, annealing, initializeRandom);
    
    /* Run MBO */
    runMBO(mbos, lowerBound);
    
   
    
    
    /* Sort by geography to match to color */
    indexedDouble *latitude = calloc(lcount,sizeof(indexedDouble));

    for(int d = 0; d < lcount; d++){
        latitude[d].dist = mbos.districtCentroids[2*d];
        latitude[d].index = d;
    }
    
    qsort(&latitude[0], lcount, sizeof(indexedDouble), cmpIndexedDouble);

    
    /* Print final populations to terminal */
    int *output = calloc(pcount*(lcount + 1),sizeof(int));
    
    int popCount[lcount];
    memset(popCount,0,lcount*sizeof(int));

    int relabelings[lcount];
    memset(relabelings,0,lcount*sizeof(int));

    
    for(int i = 0; i < pcount; i++){
        output[(lcount + 1)*i] = stateData[i].id;
        for(int d = 0; d < lcount; d++){
            output[(lcount + 1)*i + d + 1] = mbos.labelMasses[lcount*i + d];
            popCount[d] += mbos.labelMasses[lcount*i + d];
        }
    }
    
    for(int d = 0; d < lcount; d++){
        relabelings[d] = latitude[d].index;

        if(verbose){
            printf("Population %d: %d\n",d,popCount[d]);
        }
    }
    
    
    /* Write energy output to file */
    writeDataToFile(energyAtStep,"energy",maxIters*sizeof(double));

    /* Write labeling output to file */
    writeDataToFile(output, "mbo_output",pcount*(lcount+1)*sizeof(int));
    
    /* Write labeling output to file */
    writeDataToFile(relabelings, "geosorted_labels",lcount*sizeof(int));
    
    /* Write color file for districts */
    
    /*
    double colors[3*lcount];
    for(int i = 0; i < lcount ; i++){
        
        colors[3*i] = rand()/(RAND_MAX*1.0);
        colors[3*i + 1] = rand()/(RAND_MAX*1.0);
        colors[3*i + 2] = rand()/(RAND_MAX*1.0);
        
    }

    
    
    int colorAssignments[lcount];
    assignDistrictColors(colors, mbos.districtCentroids, colorAssignments, lcount);
    writeDataToFile(colorAssignments, "district_colors", lcount*sizeof(double));
*/
    /* Clean-up */
    destroyMBOStruct(mbos); //cleanup
    free(stateData);
    free(output);
    free(unitDistricts);
    free(energyAtStep);

}









