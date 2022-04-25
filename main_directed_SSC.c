//
//  CI_directed_SCC.c
//
//  Created by Gino Del Ferraro on 2017-08-12.
//  Copyright © 2017 Gino Del Ferraro. All rights reserved.
//

// HOW TO RUN: ./a.out NoN.txt L M N
//INPUT:
// 1 - NoN.txt: node i, node j, jij or whathever, mod i, mod j
// 2 - L = CI ball length
// 3 - M = number of rows in the file NoN.txt (# of links)
// 4 - N = number of nodes

#include "CI_directed_SSC.h"

#define dir_INPUT "/Volumes/LAB 1/BRAIN/LTP_Rats/2/NoN_dir_undir" //directory containing the input files
#define dir_OUTPUT " " // output directory

int main(int argc, const char * argv[]) {
    
    chdir(dir_INPUT); // Move to Input the directory
    
    int **SCC, **A,**Ain,**Aout,**Asym,**A_und, *queue, *length,i,*CI_score,k,j;
    int *interlink,node_maxCI,lapse,neigh, tot_cl, **size_comp,n_removed;
    double GC,size_SCC;
    struct node *v;
    struct vertex *w;
    struct stack *s;
    
    lapse = 1;      // # of runs
    node_maxCI = 1; // node with max CI, initialize to 1
    n_removed = 0;  // number of removed nodes
    count = 0;
    
    input_parameters(argc,argv); // get input parameters from stdin 
    allocate_memory(&v,&SCC,&w,&s,&queue,&length,&interlink,&CI_score);   // allocate memory
    read_Aij_LinkType(argv,&A,&Ain,&Aout,&Asym,&A_und,w);   // read Adjacency Aij and create the matrix that are needed
    initialize(v,w);        // initialize the nodes as UNVISITED
    
    compute_interlinks(A,Ain,interlink);  // compute how many in-going interlinks for each node
    s->top = -1;     // set the pointer to the stack in SCC to the initial value
 
   // chdir(dir_OUTPUT); // Move to Output the directory
    
    draw_graph(A,Aout,w,1);
    
    for (i=1; i<=N; i++){
        tarjan(Aout,v,w,i,SCC,s); //  compute the SCC component with a tarjan algorithm
    }
    empty_stack(s);
    reset(v);           // reset the initial values of the vector v to run another tarjan at next lapse
    
    size_SCC = bubble_sort(SCC); // sort the sizes of SCC and get the max value
    empty_SCC(SCC); // set to zeros the vector of SCC
    
    print_out(n_removed,size_SCC,0,w);  // print out results on stdin

    
    /* run till the SCC is larger than a certain value */
    while(node_maxCI != 0 && size_SCC > 0.05 ){
    
        draw_graph(A,Aout,w,lapse);
   //     printf("\n\n ------------- lapse =  % d ----------------------------\n",lapse);
        
        /* FIND CI SCORES */
        for (i=1;i<=N;i++){  // for all the nodes
            if(w[i].n == IN && v[i].stack == OUTSTACK){ // if the node has not been removed from the network and if it was part of the SCC at the previous lapse
                CI_score[i] = get_CI_NoN(i,w,Aout,A,interlink,queue,length); // compute CI of the node (only of nodes still part of the SCC
              //  printf("i  = %d, CIscore = %d\n\n",i,CI_score[i]);
            }
        }
        node_maxCI = get_max_CI_node(CI_score);  // get node with CI max
        remove_node(w,node_maxCI,A,Ain,Aout,CI_score,interlink); // remove node with max CI from the network
        printf("removed node = %d\n",node_maxCI);
        n_removed++;
        reset_CI(CI_score); // set the CI scores to zero (of all the nodes, specially of those not in the SCC anymore)
        /* ------------------ */
        
        /* FIND SCC - TARJAN */
        count = 0;
        for (i=1; i<=N; i++){
            tarjan(Aout,v,w,i,SCC,s); // compute the SCC component with a tarjan algorithm
        }
        empty_stack(s); // empty the stac
        reset(v); // reset the initial values of the vector v to run another tarjan at next lapse
        size_SCC = bubble_sort(SCC);
        empty_SCC(SCC);
        
        
        print_out(n_removed,size_SCC,node_maxCI,w); // print out results
        
        lapse++;
    }
    
    
    return 0;
}


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

/*
 
 for(i=1; i<= N; i++){
 printf("i = %d --> degree %d\n",i,LinkType[i][0]);
 for(k=1;k<=LinkType[i][0];k++){
 printf(" %d\t",LinkType[i][k]);
 }
 printf("\n\n");
 }
 
 // for (i=1;i<=N;i++){
 //      tarjan(A,v,i,SCC);
 //  }
 
 
 
 
 
 
 
 */

