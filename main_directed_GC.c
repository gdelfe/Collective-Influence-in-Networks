//
//  CI_directed.c
//
//  Created by Gino Del Ferraro on 2017-08-12.
//  Copyright © 2017 Gino Del Ferraro. All rights reserved.
//

// HOW TO RUN: ./a.out NoN.txt L M
//INPUT:
// 1 - NoN.txt: node i, node j, jij or whathever, mod i, mod j
// 2 - L = CI ball length
// 3 - M = number of rows in the file NoN.txt

// N,i.e. # of nodes, is defined in CI_directed_v3.h

#include "CI_directed_GC.h"

#define dir_INPUT "/Volumes/LAB/BRAIN/LTP_Rats/1/NoN" //directory containing the input files
#define dir_OUTPUT "/Volumes/LAB/BRAIN/LTP_Rats/1/CI" // output directory

int main(int argc, const char * argv[]) {
    
    chdir(dir_INPUT); // Move to Input the directory
    
    int **SCC, **A,**Ain,**Aout,**Asym,**A_und, *queue, *length,i,*CI_score,k,j;
    int *interlink,node_maxCI,lapse,neigh, tot_cl, **size_comp,n_removed;
    double GC;
    struct node *v;
    struct vertex *w;

    s.top = -1;     // for the stack in SCC
    lapse = 1;      // # of runs
    node_maxCI = 1; // node with max CI, initialize to 1
    n_removed = 0;  // number of removed nodes
    
    
    input_parameters(argc,argv); // get input parameters from stdin 
    allocate_memory(&v,&SCC,&w,&queue,&length,&interlink,&CI_score);   // allocate memory
    read_Aij_LinkType(argv,&A,&Ain,&Aout,&Asym,&A_und,w);   // read Adjacency Aij and create the matrix that are needed
    initialize(v,w);        // initialize the nodes as UNVISITED
    
    compute_interlinks(A,Ain,interlink);  // compute how many in-going interlinks for each node
    
   // chdir(dir_OUTPUT); // Move to Output the directory
    
    /* compute GC with zero nodes removed */
    tot_cl = compute_graph_components(w,A_und,queue);  // compute graph components
    GC = count_size_clusters(w,&size_comp,tot_cl,0,0);    // get the size of the GC
    
    /* run till the GC is larger than a certain value */
    while(node_maxCI != 0 && GC > 0.05 ){
        
//        draw_graph(A,Aout,w,lapse);
//        printf("\n\n ------------- lapse =  % d ----------------------------\n",lapse);
        
        
        for (i=1;i<=N;i++){  // for all the nodes
            if(w[i].n == IN){ // if the node has not been removed from the network
                CI_score[i] = get_CI_NoN(i,w,Aout,A,interlink,queue,length); // compute CI of the node
             //   printf("i  = %d, CIscore = %d\n\n",i,CI_score[i]);
            }
        }
        
        node_maxCI = get_max_CI_node(CI_score);  // get node with CI max
        remove_node(w,node_maxCI,A,Ain,Aout,CI_score,interlink); // remove node with max CI from the network
        n_removed++;
        
        tot_cl = compute_graph_components(w,A_und,queue);  // compute graph components
      //  printf("graph components = %d\n",tot_cl);
        GC = count_size_clusters(w,&size_comp,tot_cl,n_removed,node_maxCI);    // get the size of the GC
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

