//
//  stack.c
//  CI_directed
//
//  Created by gino on 2017-08-12.
//  Copyright © 2017 gino. All rights reserved.
//

// HOW TO RUN: 1. list_of_nodes (see below for format) 2. L (length CI ball) 3. M (numb rows in the list of node file)
// 1. list_of_nodes format: node i, node j, Jij, mod i, mod j
// CHECK: #define N in library is the number of total nodes

// OUTPUT on stdout:   number of node removed/ fraction of nodes removed / GC fraction / label of CI node / module of CI node 


#include "CI_undirected.h"

#define directory_INPUT "/Users/gino/Dropbox (Brain Conquistadores)/OPERATING_ROOM/SPEECH_STUDY/NO_ARREST/23.-35414460/CI" //directory containing the input files

int main(int argc, const char * argv[]) {
    
    chdir(directory_INPUT); // Move to Input the directory
    
    int **SCC, **A,**Ain,**Aout,**Asym,**A_und, *queue, *length,i,*CI_score,k,j;
    int *interlink,node_maxCI,lapse,neigh, tot_cl, **size_comp,n_removed;
    double GC;
    struct node *v;
    struct vertex *w;
    
    s.top = -1;
    lapse = 1;
    node_maxCI = 1;
    n_removed = 0;
    
    
    input_parameters(argc,argv);
    
    
    allocate_memory(&v,&SCC,&w,&queue,&length,&interlink,&CI_score);   // allocate memory
    read_Aij_LinkType(argv,&A,&Ain,&Aout,&Asym,&A_und,w);   // read Adjacency Aij and create the matrix that are needed

    initialize(v,w);        // initialize the nodes as UNVISITED
    
    compute_interlinks(Asym,A_und,interlink);  // compute how many in-going interlinks for each node
    
    
    /* compute GC with zero nodes removed */
    tot_cl = compute_graph_components(w,A_und,queue);  // compute graph components
    
   // printf("tot cl = %d\n",tot_cl);
    
    GC = count_size_clusters(w,&size_comp,tot_cl,0,0);    // get the size of the GC
    
    // printf("GC = %d\n",GC);
   // printf("N = %d\n",N);
    
    /* run till the GC is larger than a certain value */
    while(node_maxCI != 0 && GC > 0.03){
        
   //    draw_graph(Asym,A_und,w,lapse);
   //     printf("\n\n ------------- lapse =  % d ----------------------------\n",lapse);
    
        
        for (i=1;i<=N;i++){  // for all the nodes
            if(w[i].n == IN){ // if the node has not been removed from the network
                CI_score[i] = get_CI_NoN(i,w,A_und,Asym,interlink,queue,length); // compute CI of the node
//                printf("i  = %d, CIscore = %d\n\n",i,CI_score[i]);
            }
        }
        
        
        node_maxCI = get_max_CI_node(CI_score);  // get node with CI max
    //    printf("removed node is %d\n",node_maxCI);
        remove_node(w,node_maxCI,Asym,A_und,A_und,CI_score,interlink); // remove node with max CI from the network
        n_removed++;
//        printf("------------------- removed node %d\n",node_maxCI);

        
        tot_cl = compute_graph_components(w,A_und,queue);  // compute graph components
      //  printf("graph components = %d\n",tot_cl);
        GC = count_size_clusters(w,&size_comp,tot_cl,n_removed,node_maxCI);    // get the size of the GC
        lapse++;
    }
  
    return 0;
}


    


/*
 
 
 for(i=1; i<= N; i++){
 printf("i = %d -->",i);
 for(k=1;k<=A_und[i][0];k++){
 printf(" %d\t",A_und[i][k]);
 }
 printf("\n");
 }
 
 for(i=1; i<= N; i++){
 for(k=1;k<=N;k++){
 printf("%d\t",Asym[i][k]);
 }
 printf("\n");
 }
 
 for (i=1;i<=N;i++){
 printf("%d in-degree  = %d, out-degree = %d, interlink = %d \n",i,w[i].deg_in,w[i].deg_out,interlink[i]);
 }
 
 // for (i=1;i<=N;i++){
 //      tarjan(A,v,i,SCC);
 //  }
 
 
 */

