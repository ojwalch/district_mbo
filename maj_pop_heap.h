#include <stdlib.h>
#ifndef BRANCHING_NUMBER
#define BRANCHING_NUMBER 8
#endif
#ifndef MAJ_SIMPLE_IMPLICIT_HEAP
#define MAJ_SIMPLE_IMPLICIT_HEAP

typedef struct{
    k_type key;
    int originalIndex;
    int mass;
    int lbl;
    
}pop_node;

typedef struct{
    pop_node *root;
    int count;
    int totalMass;
    int limit;
}pop_heap;

static void pop_heap_swap_nodes(pop_heap *heap, int nodeIndex1, int nodeIndex2){
    
    int i1 = heap->root[nodeIndex1].originalIndex;
    int i2 = heap->root[nodeIndex2].originalIndex;
    k_type key1 = heap->root[nodeIndex1].key;
    k_type key2 = heap->root[nodeIndex2].key;
    int mass1 = heap->root[nodeIndex1].mass;
    int mass2 = heap->root[nodeIndex2].mass;
    int lbl1 = heap->root[nodeIndex1].lbl;
    int lbl2 = heap->root[nodeIndex2].lbl;
    heap->root[nodeIndex1].originalIndex = i2;
    heap->root[nodeIndex2].originalIndex = i1;
    heap->root[nodeIndex1].key = key2;
    heap->root[nodeIndex2].key = key1;
    heap->root[nodeIndex1].mass = mass2;
    heap->root[nodeIndex2].mass = mass1;
    heap->root[nodeIndex1].lbl=lbl2;
    heap->root[nodeIndex2].lbl=lbl1;
    
}

static void pop_heap_add_node_to_bottom(pop_heap *heap, int originalIndex, k_type key, int mass, int lbl){
    
    int count = heap->count;
    heap->root[count].originalIndex = originalIndex;
    heap->root[count].key = key;
    heap->root[count].mass = mass;
    heap->root[count].lbl = lbl;
    heap->count++;
    heap->totalMass += mass;
}

static void pop_heap_bubble_up(pop_heap *heap, int nodeIndex){
    
    while (nodeIndex > 0) {
        k_type myKey = heap->root[nodeIndex].key;
        int parentIndex = nodeIndex/BRANCHING_NUMBER - ((nodeIndex % BRANCHING_NUMBER) == 0);
        k_type parentKey = heap->root[parentIndex].key;
        if(myKey < parentKey){
            pop_heap_swap_nodes(heap,nodeIndex,parentIndex);
            nodeIndex = parentIndex;
        }else{
            break;
        }
    }
}

static void pop_heap_push_down(pop_heap *heap, int nodeIndex){
    
    int count = heap->count;
    int childIndex = nodeIndex*BRANCHING_NUMBER+1;
    while(childIndex < count){
        int minIndex = childIndex;
        k_type myKey = heap->root[nodeIndex].key;
        k_type min = myKey;
        for(int i = 0; i < BRANCHING_NUMBER; i++){
            if(childIndex + i < count){
                k_type childKey = heap->root[childIndex+i].key;
                if(childKey < min){
                    min = childKey;
                    minIndex = childIndex + i;
                }
            }
        }
        if(min < myKey){
            pop_heap_swap_nodes(heap,nodeIndex,minIndex);
            nodeIndex = minIndex;
            childIndex = nodeIndex*BRANCHING_NUMBER+1;
        }else{
            break;
        }
        
    }
}


static void pop_heap_delete_min(pop_heap *heap){
    
    int count = heap->count;
    int mass = heap->root[0].mass;
    pop_heap_swap_nodes(heap,count-1,0);
    heap->count--;
    heap->totalMass -= mass;
    pop_heap_push_down(heap,0);
}



pop_node extract_min(pop_heap *heap){
    
    pop_node n = heap->root[0];
    pop_heap_delete_min(heap);
    return n;
    
}

static void remove_excess_mass(pop_heap *heap, int gap){
    
    while(heap->totalMass > gap){
        int minMass = heap->root[0].mass;
        int tot = heap->totalMass;
        
        if(tot - minMass >= gap){
            pop_heap_delete_min(heap);
        }else{
            int excess = tot-gap;
            heap->root[0].mass -= excess;
            heap->totalMass -= excess;
            
        }
    }
}

/* */ 
static void pop_heap_consider_candidate(pop_heap *heap, int nodeIndex, k_type newKey, int mass, int lbl, int gap){
    int tot = heap->totalMass;
    k_type min = heap->root[0].key;
    
    if(tot < gap){
        pop_heap_add_node_to_bottom(heap,nodeIndex,newKey, mass, lbl);
        pop_heap_bubble_up(heap,heap->count-1);
    }else if(newKey>min){
        pop_heap_add_node_to_bottom(heap,nodeIndex,newKey, mass, lbl);
        pop_heap_bubble_up(heap,heap->count-1);
    }
    
    if(heap->totalMass>gap){
        remove_excess_mass(heap, gap);
    }
}

static pop_heap pop_heap_create_empty_heap(int pcount){
    pop_heap heap;
    heap.count = 0;
    heap.totalMass = 0;
    heap.root = malloc(pcount*sizeof(pop_node));
    heap.limit = pcount;
    return heap;
}


void pop_heap_reset(pop_heap *heap){
    heap->count=0;
    heap->totalMass=0;
}



static void pop_heap_destroy_heap(pop_heap *heap){
    free(heap->root);
}


#endif






