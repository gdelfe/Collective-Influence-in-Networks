//
//  CI_directed.h
//  CI_directed
//
// This code uses the GC = Giant Component of the network as a measure of network connectivity 
// The codes can only handle networks with bidirectional links
//
//  Created by gino on 2017-08-10.
//  Copyright © 2017 gino. All rights reserved.
//

#ifndef CI_directed_h
#define CI_directed_h


#endif /* CI_directed_h */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <unistd.h>

typedef enum{OFF, ON} off_on;
typedef enum{UNVISITED,VISITED} visited_unvisited;
typedef long int int_t;
typedef double dbl_t;

#define N 1590
#define MIN(a,b) (((a)<(b))?(a):(b))

/*------------*/

typedef struct stack
{
    int stk[N]; // array to store the queue
    int top;          // array index which refers to the top element in the queue
} STACK;

STACK s;

/*------------*/

// Structure node, used to run SCC, trajan algorithm, DFS
struct node{
    
    int ind;
    int lowlink;
    visited_unvisited status; // two states: the node has been visited/unvisited
    enum{OUTSTACK,INSTACK} stack;  // two states: the node is/is not on the stack
};

/*------------*/

// Structure vertex, used to run BFS to rank CI nodes
struct vertex{
    
    enum{OUT,IN} n;
    visited_unvisited status; // two states: the node has been visited/unvisited
    int deg_in, deg_out;  // current kin/kout which is update during the nodes removal
    int LabCluster;
    int mod;
};

//int N; // tot numb of nodes
int index1, count, L, M;


/*                                 **************                              */
/*                          MODULE FUNCTIONS DEFINITIONS                       */
/*                                 **************                              */


/* =========================================================================== */
/*                              GET INPUT PARAMETERS                           */
/* =========================================================================== */

/*  get parameter from stdout and convert them in the right format */
/*INPUT MUST BE: ./a.out NoN.txt CI_length numb_of_nodes_in_NoN.txt */

void input_parameters(int argc, const char **argv){
    
    if(argc != 4 ){
        printf("\nWrong number of input parameters! - EXIT \n\n");
        exit(errno);
    }
    if(argc == 4){
        L = atoi(argv[2]);  // Length of the CI ball
        M = atoi(argv[3]);  // number of rows in the NoN.txt file
    }
    
    return;
}

/* =========================================================================== */
/*                                  INPUT                                      */
/* =========================================================================== */

// allocate memory for J[i][j] and A[i][j] and assign value to Aij
void read_Aij_LinkType(const char **argv, int ***A, int ***Ain, int ***Aout, int ***Asym, int ***A_und, struct vertex *w){
    
    int i,j,node_i,node_j,mod_i,mod_j,ind,ind_sym,degree,degree_tot;
    int J;
    
    char filename[101];
    
    FILE *fpList;; // pointer to input file
    
    sprintf(filename,"%s",argv[1]); // read name_file of the out-going A matrix from stdin
    fpList = fopen(filename,"r");
    
    
    if(fpList == NULL){
        fprintf(stderr,"PROBLEM OPENING List of Node INPUT FILE\n");
        exit(errno);
    }
    
    
    // allocate memory for the different adjacency matrices
    *A = (int **)malloc( (N + 1) * sizeof(int *));      // Adj matrix in a NxN form
    *Ain = (int **)malloc( (N + 1) * sizeof(int *));    // Adj matrix-list for in-links in the form i, j1, j2,...
    *Aout = (int **)malloc( (N + 1) * sizeof(int *));   // Adj matrix-list for out-links in the form i, j1, j2,...
    
    *Asym = (int **)malloc( (N + 1) * sizeof(int *));  // Adj matrix symmetrized to store the undirected network
    *A_und = (int **)malloc( (N + 1) * sizeof(int *));  // A undirected. Adjacency matrix-list regardless link direction in the form  i, j1,j2,...
    
    for( i = 0; i < N + 1; i++){
        (*A)[i] = (int *)malloc((N + 1) * sizeof(int));
        (*Ain)[i] = (int *)malloc((N + 1) * sizeof(int));
        (*Aout)[i] = (int *)malloc((N + 1) * sizeof(int));
        
        (*Asym)[i] = (int *)malloc((N + 1) * sizeof(int));
        (*A_und)[i] = (int *)malloc((N + 1) * sizeof(int));
    }
    
    // ** read List of Nodes. File tipe: node_i, node_j, Jij, module_i, module_j
    // A[node_i][node_j] is a N+1 x N+1 matrix. It's non-zero only between 1 and N
    
    for (i = 0; i < M; i++){
        
        fscanf(fpList, "%d%d%d%d%d",&node_i,&node_j,&J,&mod_i,&mod_j); // read in-degree
        
        if(node_i > N || node_j > N){
            printf("ERROR ALLOCATING \n");
            exit(1);
        }
            
        if(mod_i == mod_j){ // if node_i and node_j are in the same cluster
            (*A)[node_i][node_j] = 1;
            
            (*Asym)[node_i][node_j] = 1;
            (*Asym)[node_j][node_i] = 1; // symmetrize
            
            w[node_i].mod = mod_i; // assign the modules
            w[node_j].mod = mod_j;
            
            
        }
        else{         // if node_i and node_j are in a different cluster
            (*A)[node_i][node_j] = 2;
           
            (*Asym)[node_i][node_j] = 2;
            (*Asym)[node_j][node_i] = 2; // symmetrize
            
            w[node_i].mod = mod_i; // assign the modules
            w[node_j].mod = mod_j;
            

        }
    }
    
    // ** Construct Adjacency list for the out-links
    for (i = 1; i <= N; i++){
        ind = 1;
        degree = 0;
        ind_sym = 1;
        degree_tot = 0;
        for(j = 1; j <= N; j++){ // for all the i-out-neighbours
            
            // out-link Adjacency list
            if((*Asym)[i][j]!=0){
                (*Aout)[i][ind] = j;
                
                degree++;
                ind++;
            }
            (*Aout)[i][0] = degree; // assign out-degree
            w[i].deg_out = degree;
            
            // Symmetrized Adjacency matrix
            if((*Asym)[i][j]!=0){
                (*A_und)[i][ind_sym] = j;
                
                degree_tot++;
                ind_sym++;
            }
            (*A_und)[i][0] = degree_tot; // assign degree
        }
    }
    
    
    // ** Construct Adjacency matrix for the in-links
    for (j = 1; j <= N; j++){
        ind = 1;
        degree = 0;
        for(i = 1; i <= N; i++){
            if((*Asym)[i][j]!=0){
                (*Ain)[j][ind] = i;
                degree++;
                ind++;
            }
            (*Ain)[j][0] = degree; // assign in-degree
            w[j].deg_in = degree;
        }
    }
    
    
    return;
    
}

/* =========================================================================== */
/*                    Compute NUMB of INTERLINKS per NODE                      */
/* =========================================================================== */


// count how many in-going interlinks each node has
// LinkType is the matrix which stores the type of links.
// LinkType = 2 is an interlink

void compute_interlinks(int **A,int **Ain,int *interlink){
    
    int i,j,neigh;
    
    for(i = 1; i <= N; i++){
        for(j = 1; j <= Ain[i][0]; j++){ // for all the in-going neighbours
            neigh = Ain[i][j];
            if(A[neigh][i] == 2)  // check if it's a control link (interlink)
                interlink[i] += 1;
        }
    }
    return;
}




/* =========================================================================== */
/*                             ALLOCATE MEMORY                                 */
/* =========================================================================== */

void allocate_memory(struct node **v, int ***SCC, struct vertex **w, int **queue, int **length, int **interlink, int **CI_score){
    
    int i;
    *v = (struct node *)malloc((N + 1) * sizeof(struct node));  // allocate memory for N struct variable
    *w = (struct vertex *)malloc((N + 1) * sizeof(struct vertex));  // allocate memory for N struct variable
    *queue = (int *)calloc(N + 1,sizeof(int));
    *length = (int *)calloc(N + 1,sizeof(int));
    *interlink = (int *)calloc(N + 1,sizeof(int));
    *CI_score = (int *)calloc(N + 1,sizeof(int));
    
    *SCC = (int **)calloc(N,sizeof(int *));
    for (i=0; i<N; i++)
        (*SCC)[i] = (int *)calloc(2,sizeof(int));
    
    
    return;
    
}

/* =========================================================================== */
/*                             STACK FUNCTIONS                                 */
/* =========================================================================== */

/*  Function to add an element to the stack */
void push (int vertex_lab){
    
    //  int num;
    if (s.top == (N - 1)){
        
        printf ("Stack is Full\n");
        return;
    }
    else{
        
        //  printf ("Enter the element to be pushed\n");
        // scanf ("%d", &num);
        s.top = s.top + 1;
        s.stk[s.top] = vertex_lab;
        //      printf ("pushed element is = %d\n", s.stk[s.top]);
        
    }
    return;
}

/* =========================================================================== */

/*  Function to delete an element from the stack */
int pop (void){
    
    int node;
    if (s.top == - 1){
        
        printf ("Pop function: Stack is Empty\n");
        return (s.top);
    }
    else{
        
        node = s.stk[s.top];
        //   printf ("*********** popped element is = %d\n", s.stk[s.top]);
        s.top = s.top - 1;
    }
    return (node);
}

/* =========================================================================== */

/*  Function to display the status of the stack */
void display (void){
    
    int i;
    if (s.top == -1){
        printf ("Display function: Stack is empty\n");
        return;
    }
    else{
        printf ("\n The status of the stack is \n");
        for (i = s.top; i >= 0; i--)
        {
            printf ("%d\n", s.stk[i]);
        }
    }
    printf ("\n");
}

/* =========================================================================== */
/*                             TARJAN ALGORITHM                                */
/* =========================================================================== */


void tarjan(int **Aout, struct node *v, int i, int ** SCC){
    
    int k,neigh,n;
    
    if (v[i].status == UNVISITED){ // if the node was unvisited add it to the stack
        
        v[i].ind = index1;
        v[i].lowlink = index1;
        index1++;
        
        push(i);   // add node i to the stack
        v[i].stack = INSTACK;  // label node as in stack
        v[i].status = VISITED;   // label node as visited
        
        for (k =1; k<= Aout[i][0]; k++){       // for all the descendant
            neigh = Aout[i][k];
            
            // if the descendent was unvisited
            if (v[neigh].status == UNVISITED){
                tarjan(Aout,v,neigh,SCC);          // start spanning from that descendent (DFS)
                v[i].lowlink = MIN(v[i].lowlink,v[neigh].lowlink);
                
            }
            // if the neighbours has been visited and it's still on the stack
            else if (v[neigh].status == VISITED && v[neigh].stack == INSTACK){
                v[i].lowlink = MIN(v[i].lowlink,v[neigh].ind);
            }
            
        }
        
        // if the root and the head of the path are the same
        if (v[i].lowlink == v[i].ind){
            
            // Start a new SCC //
            SCC[count][0] = count + 1;          // SCC[0][count] keeps track of the label of the component
            do{               // pop out all the descendent of i and i itself from the stack
                
                n = pop();          // pop out the top node from the stack
                v[n].stack = OUTSTACK;  // label the node as out of the stack
                SCC[count][1] += 1;    // increment the size of the component: SCC[1][count] keeps track of the SCC component size
                
                
            }while(n != i); // while we are not in the root
            count++;
        }
    }
    
    return;
}

/* =========================================================================== */
/*                         INITIALIZE STRUCTURE of NODES                       */
/* =========================================================================== */

void initialize(struct node *v, struct vertex *w){
    
    int i;
    for (i=1; i<=N; i++){
        v[i].ind = 1;
        v[i].lowlink = i;
        v[i].status = UNVISITED;
        
        w[i].n = IN;
        w[i].status = UNVISITED;
        
    }
    return;
}

/* =========================================================================== */
/*                            BUBBLE SORT the SCC                              */
/* =========================================================================== */

// SCC[i][0] indicates the label of the i-th component
// SCC[i][1] indicate the size of the i-th component
void bubble_sort(int **SCC){
    
    int i,j,temp0,temp1;
    
    for (i = 0; i < count; i++){
        for (j = i+1; j < count;j++){
            
            if (SCC[i][1] < SCC[j][1])
            {
                temp1 = SCC[i][1];
                SCC[i][1] = SCC[j][1];
                SCC[j][1] = temp1;
                
                temp0 = SCC[i][0];
                SCC[i][0] = SCC[j][0];
                SCC[j][0] = temp0;
            }
        }
    }
    
    return;
}

/* =========================================================================== */
/*                            GET CI SCORE  one network                        */
/* =========================================================================== */

//get CI, one single network, for node i

int get_CI(int  i, struct vertex *w, int **Aout, int *queue, int *length) {
    
    
    int  s, degree, temp, delta, cnt, k, neigh, ind, CI;
    int *current, *end;
    
    degree = Aout[i][0];
    
    if( degree == 0 || degree == 1)
        return 0;
    
    else {
        queue[0] = i;
        w[i].status = VISITED;
        current = queue;  // current
        end = queue + 1; // end of the queue
        temp = 1;
        delta = 1;
        length[0] = 1;
        s = 1;
        cnt = 0;
        
        
        while(current != end) { // while queue is not empty
            if(s <= L) {
                degree = Aout[*current][0];
                for(k = 1; k <= degree; k++) {
                    neigh = Aout[*current][k];
                    if( (w[neigh].n == IN) && (w[neigh].status == UNVISITED)) {
                        
                        queue[temp++] = neigh;   //assign the neighbour to the queue
                        w[neigh].status = VISITED;
                        length[s] += 1;
                    }
                }
            }
            current += 1;
            end += temp - delta;
            delta = temp;
            cnt += 1;
            
            if(cnt == length[s-1]) {
                s++;
                cnt = 0;
            }
        }
        
        ind = 0;
        for(s = 0; s < L; s++) // count how many node there are in s = [0,L)
            ind += length[s];
        
        CI = 0;
//        printf("CI -- Node i = %d\n",i);
        for(k = ind; k < (ind + length[L]); k++) {
            CI += (Aout[queue[k]][0] - 1);    //(degree_k - 1)
//            printf("k = %d, queue[k] = %d, k_j = %d\n",k,queue[k],Aout[queue[k]][0]);
        }
        CI *= (Aout[i][0] - 1);       // (degree_i - 1)
//        printf("k_i = %d\n",Aout[i][0]);
        
        for(k = 0; k < temp; k++){
            w[queue[k]].status = UNVISITED;
            queue[k] = 0;           // set the queue to zero
        }
        for(s = 0; s <= L; s++)
            length[s] = 0;
        
        return CI;
    }
}


/* =========================================================================== */
/*                                CI BALL                                      */
/* =========================================================================== */

int CI_ball(int i, struct vertex *w, int **Aout, int *queue, int *length){
    
    int j,k,degree,temp,delta,cnt,ind,ktot,CI;
    int *current,*end,s,neigh;
    
    degree = w[i].deg_in + w[i].deg_out; // total degree
    
    CI = 1; // initialize First (Second) as True
    if( degree == 0 || degree == 1)
        CI = 0;
    
    else if(CI == 1){ // If the First (Second) term is not zero
        queue[0] = i;
        w[i].status = VISITED;
        current = queue;  // current
        end = queue + 1; // end of the queue
        temp = 1;
        delta = 1;
        length[0] = 1;
        s = 1;
        cnt = 0;
        
        while(current != end) { // while queue is not empty
            if(s <= L) {
                degree = Aout[*current][0];
                for(k = 1; k <= degree; k++) {
                        neigh = Aout[*current][k];
                    if( (w[neigh].n == IN) && (w[neigh].status == UNVISITED)) {
                        
                        queue[temp++] = neigh;   //assign the neighbour to the queue
                        w[neigh].status = VISITED;
                        length[s] += 1;
                    }
                }
            }
            current += 1;
            end += temp - delta;
            delta = temp;
            cnt += 1;
            
            if(cnt == length[s-1]) {
                s++;
                cnt = 0;
            }
        }
        
        ind = 0;
        for(s = 0; s < L; s++) // count how many node there are in s = [0,L)
            ind += length[s];
        
        CI = 0; // set CI to zero before starting compute it
//        printf("CI -- Node i = %d\n",i);
        for(k = ind; k < (ind + length[L]); k++) {
            ktot = w[queue[k]].deg_in + w[queue[k]].deg_out;
            CI += (ktot - 1);    //(degree_k - 1)
//            printf("    j = %d, queue[j] = %d, ktot_j = %d, ktot - 1 = %d \n",k,queue[k],ktot,ktot-1);
        }
        CI *= (w[i].deg_in + w[i].deg_out - 1);       // (degree_i - 1)
//        printf("    ktot_i = %d  ktot_i -1 = %d\n",w[i].deg_in + w[i].deg_out,w[i].deg_in + w[i].deg_out-1);
//        printf("-- CI  = %d\n",CI);
        
        
        for(k = 0; k < temp; k++){
            w[queue[k]].status = UNVISITED;  // set the previously visited node as unvisited
            queue[k] = 0;           // set the queue to zero
        }
        for(s = 0; s <= L; s++)
            length[s] = 0;   // remove all the elements from the frontier at each distance 0 <= s <= L
        
    }
    
    return (CI);
}

/* =========================================================================== */
/*                               GET CI SCORE  on NoN                          */
/* =========================================================================== */

//get CI, NoN, for node i
// This function makes use of the function CI_ball()


int get_CI_NoN(int  i, struct vertex *w, int **Aout, int **A, int *interlink, int *queue, int *length) {
    
    
    int  degree, temp, delta, cnt, k, ind, sum, CI_first, CI_second, ktot;
    int *current, *end, j, neigh_j;
    
    /* -------------------------------------- */
    /* First term in CI NoN formula on the RHS */
    
    CI_first = CI_ball(i,w,Aout,queue,length);
    
    /* -------------------------------------- */
    /* Second term in CI NoN formula on the RHS */
    sum = 0;
    for(j = 1; j<= Aout[i][0]; j++){ // for all the i-neighbours
        neigh_j = Aout[i][j];
        if(A[i][neigh_j] == 2 && interlink[neigh_j] == 1){ // check if the link is an interlink and if it's the unique interlink
            
//            printf("-- CI second partial  -- node = %d\n",neigh_j);
            CI_second = CI_ball(neigh_j,w,Aout,queue,length);
//            printf("-- CI second partial  = %d\n",CI_second);
            sum += CI_second; // sum over the j neighbours with only one outgoing interlink
        }
        
    }
//    printf("first =  %d, second = %d, CI totale  = %d\n\n",CI_first,sum,CI_first + sum);
    
    
    return (CI_first + sum);
}

/* =========================================================================== */
/*                         get max CI SCORE                                    */
/* =========================================================================== */

int get_max_CI_node(int *CI_score){
    
    int max,node_maxCI, i;
    
    max = 0;
    node_maxCI = 0;
    for(i = 1; i<=N; i++){
        if(CI_score[i] > max){
            max = CI_score[i];
            node_maxCI = i;
        }
    }
    return node_maxCI;
}

/* =========================================================================== */
/*                         get max CI SCORE                                    */
/* =========================================================================== */

/*int apply_rule(struct vertex *w, int **Aout, int node_maxCI){
    
    int j,node_maxCI;
    
    if(interlink[node_maxCI] == 1){
        for(j = 1; j <= Aout[node_maxCI][j]; j++){ // for all the out-going neighbours
            neigh = Aout[node_maxCI][j];
            if(w[neigh].n == IN && A[][] ==2) // if the neigh is still in the network
                
        }
    
    
    }

    return;
} */

/* =========================================================================== */
/*                         REMOVE NODE from network                            */
/* =========================================================================== */

void remove_node(struct vertex *w, int node, int **A, int **Ain, int **Aout, int *CI_score, int *interlink){
    
    int j,neigh;
    w[node].n = OUT; // remove the node from the network
    CI_score[node] = 0;   // set its CI score to zero
    
    for(j = 1; j <= Aout[node][j]; j++){ // for all the out-going neighbours
        neigh = Aout[node][j];
        if(w[neigh].n == IN) // if the neigh is still in the network
            w[neigh].deg_in--; // decrese the in-degree of the neighbour by a factor one
        
        // apply the rule
        if(w[neigh].n == IN && interlink[neigh] == 1 && A[node][neigh] == 2){  // if the neigh is IN, node i point to the neigh with an interlink and this is the only control link that the neighbour has left
//            printf("CASE --- neigh %d OUT\n",neigh);
            remove_node(w,neigh,A,Ain,Aout,CI_score,interlink); // iterate the rule on neighbours of the neighbour
        }
    }
    
    for(j = 1; j <= Ain[node][j]; j++){ // for all the in-going neighbours
        neigh = Ain[node][j];
        if(w[neigh].n == IN){ // if the neigh is still in the network
            w[neigh].deg_out--; // decrese the out-degree of the neighbour by a factor one
            
            if(A[neigh][node] == 2) // if the in-link is a control/interlink
                interlink[neigh]-= 1;  // reduce the control link of the neighbours by a factor one
        }
    }
    
    return;
}



/* =========================================================================== */
/*                         PRINT and DRAW GRAPH                                */
/* =========================================================================== */

/* Print a graph in gephViz format. Read the graph and plot it in a pdf file */
/* then remove the data file containing the graph */

void draw_graph(int **A, int **Aout, struct vertex *w, int lapse){
    
    int i,j,neigh;
    
    FILE *fp;
    char filename[101];
    char color[8][100]; // color array for the nodes
    
 /*
     sprintf(color[0],"grey");
     sprintf(color[1],"gold");
     sprintf(color[2],"brown1");
     sprintf(color[3],"cyan");
     sprintf(color[4],"ghostwhite");
     sprintf(color[5],"chocolate1");
     sprintf(color[6],"chartreuse1");
     sprintf(color[7],"deeppink");
   */
    
    sprintf(filename,"graph_%d.dat",lapse); // read name_file from stdin, file 1 or 2 depending on var
    fp = fopen(filename,"w");
    
    
    fprintf(fp,"digraph G {\n");
    
  /*     fprintf(fp,"    node[style=filled]\n");
     for (i=0; i<N; i++) { //print different colors for different clusters
     fprintf(fp,"    node[fillcolor=\"%s\"] %d\n",color[list_tot[i][1]],i);
     }
    */
    
    for (i=1; i<=N; i++) { // print edges
        
        if(w[i].n == IN){ // if the node is still in the network
            if(w[i].deg_out==0) // if degree is zero
                fprintf(fp," %d\n",i);  // print just the node
            for (j=1; j<= Aout[i][0]; j++) {  // for all the neighbours
                
                neigh = Aout[i][j];
                if(w[neigh].n == IN){ // if the neighbour is still in the network
                fprintf(fp," %d -> %d ",i,Aout[i][j]);  // i -> j
                
                neigh = Aout[i][j];
                if(A[i][neigh] == 2)
                    fprintf(fp,"[color = \"blue\",style = \"dashed\", arrowhead = \"empty\"]\n");
                else fprintf(fp,"\n");
                }
            }
        }
        
        
    }
    fprintf(fp,"}\n");
    fclose(fp);
    
    sprintf(filename,"dot -Tpdf graph_%d.dat -o graph_%d.pdf",lapse,lapse);
    system(filename);
    sprintf(filename,"rm graph_%d.dat",lapse);
    system(filename);
    
    return ;
    
}



/* =========================================================================== */
/*         COMPUTE THE GIANT COMPONENT - BREADTH FIRST SEARCH                  */
/* =========================================================================== */

/* Breadth First Search algorithm. Label each node with the cluster label of belonging */

 int compute_graph_components(struct vertex *w, int **A_und, int *queue){
    
    int i,j,temp,delta,neigh,Label_Cluster;
    int *current,*end;
    
    for (i = 1; i<= N; i++) {
        w[i].LabCluster = 0;
    }
    
    Label_Cluster = 0;
    
    for(i = 1; i <= N; i++){
        
        if(w[i].LabCluster == 0 && (w[i].n == IN)){ // if the node has not been already visited and it's still on the network
            
            Label_Cluster++ ;      // increment the total number of clusters
            
            w[i].LabCluster = Label_Cluster;
            
            queue[0] = i;
            w[i].status = VISITED;
            
            current = queue;    // pointer to the current element in the queue
            end = queue + 1;    // pointer to the last element added in the queue
            
            temp = 1;
            delta = 1;
            
            while (current != end) {  // while queue is not empty
                for (j=1; j<= A_und[*current][0]; j++) { // for all the neighbourds of i
                    
                    neigh = A_und[*current][j];
                    if( (w[neigh].n == IN) && (w[neigh].status == UNVISITED)){
                        
                        w[neigh].LabCluster = Label_Cluster;
                        w[neigh].status = VISITED;
                        
                        queue[temp] = neigh; //assign the neighbour to the queue
                        temp++;
                        
                    }
                }
                
                current += 1; // point to next node in the queue
                end += temp - delta; // point to the end of the queue
                delta = temp;
            }
        }
    }
     
     for(i = 1; i <= N; i++){
         queue[i] = 0; // set the queue to zero
         w[i].status = UNVISITED;
     }
     
    
    return Label_Cluster; // tot number of clusters
}


/* =========================================================================== */
/*                      COMPUTE SIZE OF EACH COMPONENT                         */
/* =========================================================================== */


/* Count how many nodes there are in each cluster. Put in an array[M][2]
 the label of the cluster and the size of it sorted descending based on the size */

double count_size_clusters(struct vertex *w, int ***size_comp, int tot_num_cl, int n_removed, int node_maxCI){
    
    int i,j,LabClust,temp0,temp1;
    
    
    /* Allocate memory for a mat[M][2] where M is the total number of different clusters */
    /* size_comp[i][0]is the label of the cluster */
    /* size_comp[i][1]is the size of the cluster */
    /* Note: labeling of the cluster start from 1 to M = tot_num_cl */
    
    *size_comp = (int **)calloc((tot_num_cl+1),sizeof(int *));
    for (i=0; i<(tot_num_cl+1); i++)
        (*size_comp)[i] = (int *)calloc(2,sizeof(int));
    
    
    for (i = 1; i <= N; i++){
        
        LabClust = w[i].LabCluster;
        
        (*size_comp)[LabClust][0] = w[i].LabCluster;  // Label of the cluster
        (*size_comp)[LabClust][1] = (*size_comp)[LabClust][1] + 1;  // count how many node are in each cluster
    }
    
    
    /* sort the vector size_comp in descending order based on the size of the clusters, second column  */
    /* BUBBLE SORT */
    for (i = 1; i < tot_num_cl+1; i++){
        for (j = i+1; j < tot_num_cl+1;j++){
            if ((*size_comp)[i][1] < (*size_comp)[j][1])
            {
                temp1 = (*size_comp)[i][1];
                (*size_comp)[i][1] = (*size_comp)[j][1];
                (*size_comp)[j][1] = temp1;
                
                temp0 = (*size_comp)[i][0];
                (*size_comp)[i][0] = (*size_comp)[j][0];
                (*size_comp)[j][0] = temp0;
            }
        }
    }
    
    
  //  if(tot_num_cl == 1){ /* if there is more than one component */
        /* number of node removed/ fraction of nodes removed / GC fraction / label of CI node / module of CI node */
        printf("%d\t%lf\t%lf\t%d\t%d\n",n_removed,(double)n_removed/N,(double)(*size_comp)[1][1]/N,node_maxCI,w[node_maxCI].mod);
  //  }
    
  //  if(tot_num_cl == 2){ /* if there is more than one component */
        /* Print Label GC, Size GC, Label 2nd largest component, Size 2nd largest component, 3rd ... */
   //     printf("%d\t%lf\t%d\t%lf\t0\t0\n",(*size_comp)[1][0],(double)(*size_comp)[1][1]/N,(*size_comp)[2][0],(double)(*size_comp)[2][1]/N);
  //  }
    
    return ((double)(*size_comp)[1][1]/N);
    
}



