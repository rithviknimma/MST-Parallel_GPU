#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include <cuda.h>

#include <time.h>

#define THREADSPERBLOCK 64

// graph
struct edge{
	int v;
	int u;
	int weight;
};

struct graph{
	int num_edges;
	int num_vertices;
	struct edge* edges;
};

// bipartite graph
struct b_vertex_a{
	int v;
	int small_edge; // don't know what this is for
};

struct b_vertex_b{
	int e; // edge number
};

struct b_edge{
	int v;
	int u;
	int cv;
	int weight;
};

struct b_graph{
	int num_vertex_a;
	int num_vertex_b;
	int num_bipartite_edges;
	struct b_vertex_a* vertices_a;
	struct b_vertex_b* vertices_b; 
	struct b_edge* edges;
};

// strut
struct strut_edge{
    int v;
    int u; // edge number index
    int cv; // correspondent vertex
};

struct strut_u_vertex{ 
    int degree; // degree in struts
    int v1;
    int v2;
    int weight;
};

struct strut{
    int num_v; // num of bipartite vertices
    int num_u; // num of u vertices adjacent to strut edge
    int num_strut_edges; // same number as num_v
    struct strut_edge* edges; 
    struct strut_u_vertex* vertices_u; // u vertices - 0 value indicates not in strut, value > 0 indicates how many strut edges it is connected to
};

void get_graph(struct graph* og_graph, char* input);
__global__ void get_bipartite_graph(int num_edges, int num_vertices, struct edge* graphEdges, struct b_vertex_a* vetices_a, struct b_vertex_b* vetices_b, struct b_edge* bg_graphEdges) ;


__global__ void get_smallest_edges(int bp_num_edges, int num_smallest_edges, struct b_edge* bg_graphEdges, int* smallest_weights, int* smallest_edges);
__global__ void mst_edges_init(int og_num_edges, bool *mst_edges);
__global__ void get_mst_edges(int num_smallest_edges, int* smallest_edges, struct b_edge* bg_graphEdges, bool *mst_edges);
__global__ void get_num_mst(int og_num_edges, bool *mst_edges, int* num_mst);

// strut stuff
__global__ void get_strut_edges(int bg_num_vertices, int* smallest_edges, struct b_edge* bg_graphEdges, strut_edge* strut_edges);
__global__ void strut_u_init(int bg_num_vertex_b, struct strut_u_vertex* vertices_u);
__global__ void get_strut_u_degree(int num_strut_edges, strut_edge* strut_edges, struct strut_u_vertex* vertices_u);
__global__ void get_strut_u_vertices(int bg_num_edges, struct b_edge* bg_graphEdges, struct strut_u_vertex* vertices_u);
__global__ void get_zero_diff_num(int bg_num_vertex_b, struct strut_u_vertex* vertices_u, int* zero_diff_edges);

__global__ void super_vertices_init(int num_strut_vertices, int* super_vertices);
__global__ void get_new_bg_vertex_b(int num_bg_vertexb, int* super_vertices, struct strut_u_vertex* vertices_u, int* new_vertex_b, int* num_newbg_vertexb);

__global__ void prefixCopy(int prefixNum, int* old_prefix, int *new_prefix);
__global__ void getPrefixSum(int* entries, int* entriesC, int d);
__global__ void get_super_vertices(int num_strut_vertices, strut_edge* strut_edges, struct strut_u_vertex* vertices_u, int* super_vertices);
__global__ void get_new_bg_edges(int num_bg_vertex_b, int* new_bg_edges, int* prefixSum, struct strut_u_vertex* vertices_u, int* super_vertices, struct b_edge* bg_graphEdges, int * max_super_vertex);

__global__ void init_smallest_edges_weights(int num_edges, int *smallest_weights, int* smallest_edges);
/* NOTES: 
    - Remember to free all malloced and cuda malloced variables 
    - Comment out debugging statements
    - Output to file 
    - find sequential algorithm that outputs result in same way
    - compare results with that
    - Time the algorithm - put in output
    - documentation
        - read piazza and term project info for documentation
    - prep for presentation
    - submit on github - ask chonyang and email garg by thursday morning
*/

// driver
int main(int argc, char** argv){
	if(argc != 3){
		printf("mst: incorrect formatting\n");
		printf("Valid input: mst.out <Input file name> <Output file name>\n");
		return 0;
	}

	//***** ACQUIRE INPUT GRAPH *****//
	struct graph og_graph; // input
	get_graph(&og_graph, argv[1]);

	//debugging
	// printf("Graph:\n");
	// printf("vertices:%d edges:%d\n",og_graph.num_vertices, og_graph.num_edges);
	// for(int i = 0; i < og_graph.num_edges; i++){
	// 	printf("index:%d - %d   %d   %d\n", i, og_graph.edges[i].v, og_graph.edges[i].u, og_graph.edges[i].weight);
	// }

	//***** CREATE BIPARTITE GRAPH *****//
	struct b_graph bg_graph;
	bg_graph.num_vertex_a = og_graph.num_vertices;
	bg_graph.num_vertex_b = og_graph.num_edges;
	bg_graph.num_bipartite_edges = og_graph.num_edges * 2;

	// allocate GPU array
	cudaMalloc((void**) &(bg_graph.vertices_a), bg_graph.num_vertex_a * sizeof(struct b_vertex_a));
	cudaMalloc((void**) &(bg_graph.vertices_b), bg_graph.num_vertex_b * sizeof(struct b_vertex_b));
	cudaMalloc((void**) &(bg_graph.edges), bg_graph.num_bipartite_edges * sizeof(struct b_edge));

	struct edge* d_og_edges = NULL;
	cudaMalloc((void**) &(d_og_edges), og_graph.num_edges * sizeof(struct edge));
	cudaMemcpy(d_og_edges, og_graph.edges, og_graph.num_edges*sizeof(struct edge), cudaMemcpyHostToDevice);

	get_bipartite_graph<<<(og_graph.num_edges + THREADSPERBLOCK-1)/THREADSPERBLOCK, THREADSPERBLOCK>>>(og_graph.num_edges, og_graph.num_vertices, d_og_edges, bg_graph.vertices_a, bg_graph.vertices_b, bg_graph.edges);

	//    debugging
	struct b_graph debugging;
	debugging.num_vertex_a = og_graph.num_vertices;
	debugging.num_vertex_b = og_graph.num_edges;
	debugging.num_bipartite_edges = og_graph.num_edges * 2;

	debugging.vertices_a = (struct b_vertex_a*) malloc(debugging.num_vertex_a * sizeof(struct b_vertex_a));
	debugging.vertices_b = (struct b_vertex_b*) malloc(debugging.num_vertex_b * sizeof(struct b_vertex_b));
	debugging.edges = (struct b_edge*) malloc(debugging.num_bipartite_edges * sizeof(struct b_edge));

	cudaMemcpy(debugging.vertices_a, bg_graph.vertices_a, debugging.num_vertex_a * sizeof(struct b_vertex_a), cudaMemcpyDeviceToHost);
	cudaMemcpy(debugging.vertices_b, bg_graph.vertices_b, debugging.num_vertex_b * sizeof(struct b_vertex_b), cudaMemcpyDeviceToHost);
	cudaMemcpy(debugging.edges, bg_graph.edges, debugging.num_bipartite_edges * sizeof(struct b_edge), cudaMemcpyDeviceToHost);

	// printf("Bipartite Graph:\n");
	// printf("verticesA: %d, verticesB: %d, edges: %d\n", debugging.num_vertex_a, debugging.num_vertex_b, debugging.num_bipartite_edges);
	// for(int i = 0; i < debugging.num_bipartite_edges; i++){
	// 	printf("index: %d - %d   %d   %d   %d\n", i, debugging.edges[i].v, debugging.edges[i].u, debugging.edges[i].cv, debugging.edges[i].weight);
	// }
	free(debugging.vertices_a);
	free(debugging.vertices_b);
    free(debugging.edges);
    
	//***** SMALLEST EDGE WEIGHT EDGE FOR EACH VERTEX IN BG_GRAPH *****//
	int* smallest_weights = NULL;
	int* smallest_edges = NULL;

	
    //***** GET SOLUTION *****//
    bool* d_mst_edges = NULL;
    bool* mst_edges = NULL;
    mst_edges = (bool*) malloc(og_graph.num_edges * sizeof(bool));

    // don't malloc again for this variable
    cudaMalloc((void**) &(d_mst_edges), og_graph.num_edges* sizeof(bool));
    
    mst_edges_init<<<(og_graph.num_edges + THREADSPERBLOCK-1)/THREADSPERBLOCK, THREADSPERBLOCK>>>(og_graph.num_edges, d_mst_edges);

    int * solution_size = (int*) malloc (sizeof(int));
    *solution_size = 0;
    int* d_solutionSize = NULL;
    cudaMalloc((void**) &(d_solutionSize),sizeof(int));
    cudaMemcpy(d_solutionSize, solution_size, sizeof(int), cudaMemcpyHostToDevice);

    int* max_super_vertex = (int*) malloc (sizeof(int));
    *max_super_vertex = bg_graph.num_vertex_a;
    
    while(*solution_size <  (og_graph.num_vertices - 1)){
        // cudaMalloc((void**) &(smallest_weights), bg_graph.num_vertex_a * sizeof(int));
        // cudaMalloc((void**) &(smallest_edges), bg_graph.num_vertex_a * sizeof(int));
        //get_smallest_edges<<<(bg_graph.num_bipartite_edges + THREADSPERBLOCK-1)/THREADSPERBLOCK, THREADSPERBLOCK>>>(bg_graph.num_bipartite_edges, bg_graph.num_vertex_a, bg_graph.edges, smallest_weights, smallest_edges);
        cudaMalloc((void**) &(smallest_weights), *max_super_vertex * sizeof(int));
        cudaMalloc((void**) &(smallest_edges), *max_super_vertex * sizeof(int));
        init_smallest_edges_weights<<<(*max_super_vertex + THREADSPERBLOCK-1)/THREADSPERBLOCK, THREADSPERBLOCK>>>(*max_super_vertex, smallest_weights, smallest_edges);
        get_smallest_edges<<<(bg_graph.num_bipartite_edges + THREADSPERBLOCK-1)/THREADSPERBLOCK, THREADSPERBLOCK>>>(bg_graph.num_bipartite_edges,*max_super_vertex, bg_graph.edges, smallest_weights, smallest_edges);
    
        // debugging
        int* debug_smallest_weights = NULL;
        debug_smallest_weights = (int*) malloc(*max_super_vertex * sizeof(int));
        cudaMemcpy(debug_smallest_weights, smallest_weights, *max_super_vertex * sizeof(int), cudaMemcpyDeviceToHost);
    
        // for(int i = 0; i < *max_super_vertex; i++){
        //     printf("bg index of smallest weight: %d\n", debug_smallest_weights[i]);
        // }

        int* debug_smallest_edges = NULL;
        debug_smallest_edges = (int*) malloc(*max_super_vertex * sizeof(int));
        cudaMemcpy(debug_smallest_edges, smallest_edges, *max_super_vertex * sizeof(int), cudaMemcpyDeviceToHost);
    
        // for(int i = 0; i < *max_super_vertex; i++){
        //     printf("bg index of smallest edge: %d\n", debug_smallest_edges[i]);
        // }

        get_mst_edges<<<(*max_super_vertex + THREADSPERBLOCK-1)/THREADSPERBLOCK, THREADSPERBLOCK>>>(*max_super_vertex , smallest_edges, bg_graph.edges, d_mst_edges);
        
        get_num_mst<<<(og_graph.num_edges  + THREADSPERBLOCK-1)/THREADSPERBLOCK, THREADSPERBLOCK>>>(og_graph.num_edges , d_mst_edges, d_solutionSize);

        // debugging
        // printf("MST:\n");
        cudaMemcpy(mst_edges, d_mst_edges, og_graph.num_edges * sizeof(bool), cudaMemcpyDeviceToHost);
        // for(int i = 0; i < og_graph.num_edges; i++){
        //     if(mst_edges[i] == true){
        //         //printf("mst edges index: %d\n", i);
        //         printf("index: %d - %d   %d   %d\n", i, og_graph.edges[i].v, og_graph.edges[i].u, og_graph.edges[i].weight);
        //     }
        // }
        
        cudaMemcpy(solution_size, d_solutionSize,  sizeof(int), cudaMemcpyDeviceToHost);
        // printf("Num MST edges found: %d\n",*solution_size);
        
        if(*solution_size <  (og_graph.num_vertices - 1)){
            //***** GET STRUT *****//
            struct strut new_strut;
            new_strut.num_v = bg_graph.num_vertex_a;
            new_strut.num_u = bg_graph.num_vertex_b;
            new_strut.num_strut_edges = bg_graph.num_vertex_a;

            struct strut_edge* d_strut_edges = NULL; 
            
            cudaMalloc((void**) &(d_strut_edges), new_strut.num_v * sizeof(struct strut_edge));
            get_strut_edges<<<((bg_graph.num_vertex_a) + THREADSPERBLOCK-1)/THREADSPERBLOCK, THREADSPERBLOCK>>>(bg_graph.num_vertex_a, smallest_edges, bg_graph.edges, d_strut_edges);

            // debugging
            struct strut_edge* strut_edges = (struct strut_edge* ) malloc(new_strut.num_v * sizeof(struct strut_edge));
            cudaMemcpy(strut_edges, d_strut_edges, new_strut.num_v * sizeof(struct strut_edge), cudaMemcpyDeviceToHost);
            // printf("STRUT EDGES:\n");
            // for(int i = 0; i < new_strut.num_v ; i++){
            //     printf("%d   %d   %d\n", strut_edges[i].v,strut_edges[i].u, strut_edges[i].cv);
            // }

            // getting strut_u
            struct strut_u_vertex* d_vertices_u = NULL;
            cudaMalloc((void**) &(d_vertices_u), new_strut.num_u * sizeof(struct strut_u_vertex));
            strut_u_init<<<((new_strut.num_u) + THREADSPERBLOCK-1)/THREADSPERBLOCK, THREADSPERBLOCK>>>(new_strut.num_u, d_vertices_u);
            get_strut_u_degree<<<((new_strut.num_v) + THREADSPERBLOCK-1)/THREADSPERBLOCK, THREADSPERBLOCK>>>(new_strut.num_v, d_strut_edges, d_vertices_u);
        
            get_strut_u_vertices<<<((bg_graph.num_bipartite_edges) + THREADSPERBLOCK-1)/THREADSPERBLOCK, THREADSPERBLOCK>>>(bg_graph.num_bipartite_edges, bg_graph.edges, d_vertices_u);

            // debugging
            struct strut_u_vertex* vertices_u = (struct strut_u_vertex* ) malloc(new_strut.num_u * sizeof(struct strut_u_vertex));
            cudaMemcpy(vertices_u, d_vertices_u, new_strut.num_u * sizeof(struct strut_u_vertex), cudaMemcpyDeviceToHost);
            // printf("STRUT U VERTICES DEGREE:\n");
            // for(int i = 0; i < new_strut.num_u ; i++){
            //     printf("index: %d degree: %d v1: %d v2: %d\n",i, vertices_u[i].degree,  vertices_u[i].v1,  vertices_u[i].v2);
            // }

            /* ZERO DIFF */
            int* d_zero_diff_edges = NULL;
            cudaMalloc((void**) &(d_zero_diff_edges), new_strut.num_u * sizeof(int));
            get_zero_diff_num<<<((new_strut.num_u) + THREADSPERBLOCK-1)/THREADSPERBLOCK, THREADSPERBLOCK>>>(new_strut.num_u, d_vertices_u, d_zero_diff_edges);

            // debugging
            int*zero_diff_edges = (int*) malloc (sizeof(int));
            cudaMemcpy(zero_diff_edges, d_zero_diff_edges, sizeof(int), cudaMemcpyDeviceToHost);
            // printf("zero diff edges: %d\n", *zero_diff_edges);

            // /*SUPER VERTEX*/
            int* d_super_vertices = NULL;
            cudaMalloc((void**) &(d_super_vertices), bg_graph.num_vertex_a* sizeof(int));
            super_vertices_init<<<((bg_graph.num_vertex_a) + THREADSPERBLOCK-1)/THREADSPERBLOCK, THREADSPERBLOCK>>>(bg_graph.num_vertex_a, d_super_vertices);

            int* super_vertices = (int*) malloc(bg_graph.num_vertex_a* sizeof(int));
            cudaMemcpy(super_vertices, d_super_vertices,bg_graph.num_vertex_a* sizeof(int), cudaMemcpyDeviceToHost);

            for(int i = 0; i < new_strut.num_u ; i++){
                int super_vertex;
                if(vertices_u[i].degree > 0){ // if incident to strut edge
                    if(super_vertices[vertices_u[i].v1 - 1] < vertices_u[i].v1){
                        super_vertex = super_vertices[vertices_u[i].v1 - 1];
                        super_vertices[vertices_u[i].v1 - 1] = super_vertex;
                        super_vertices[vertices_u[i].v2 - 1] = super_vertex;
                    }
                    else{
                        super_vertex = vertices_u[i].v1;
                        super_vertices[vertices_u[i].v1 - 1] = super_vertex;
                        super_vertices[vertices_u[i].v2 - 1] = super_vertex;
                    }
                }
            }

            cudaMemcpy(d_super_vertices, super_vertices,bg_graph.num_vertex_a* sizeof(int), cudaMemcpyHostToDevice);

            // debugging
            // printf("Supervertices\n:");
            // for(int i = 0; i < bg_graph.num_vertex_a ; i++){
            //     printf("vertex: %d supervertex: %d\n", i+1, super_vertices[i]);
            // }

            /******** CREATING NEW BIPARTITE GRAPH **********/
            struct b_graph new_bg_graph;
            new_bg_graph.num_vertex_a = *zero_diff_edges;
            
            int* new_num_vertex_b = NULL;
            cudaMalloc((void**) &(new_num_vertex_b), sizeof(int));
            int* new_vertex_b = NULL;
            cudaMalloc((void**) &(new_vertex_b), bg_graph.num_vertex_b * sizeof(int));
            get_new_bg_vertex_b<<<((bg_graph.num_vertex_b) + THREADSPERBLOCK-1)/THREADSPERBLOCK, THREADSPERBLOCK>>>(bg_graph.num_vertex_b, d_super_vertices,d_vertices_u, new_vertex_b, new_num_vertex_b);

            // num_vertex_b and num_bipartite edges
            cudaMemcpy(&new_bg_graph.num_vertex_b, new_num_vertex_b, sizeof(int), cudaMemcpyDeviceToHost);
            new_bg_graph.num_bipartite_edges = new_bg_graph.num_vertex_b * 2;
            
            // debugging 
            int* new_vertex_b_debug = (int*)malloc( bg_graph.num_vertex_b * sizeof(int));
            cudaMemcpy(new_vertex_b_debug, new_vertex_b,  bg_graph.num_vertex_b * sizeof(int), cudaMemcpyDeviceToHost);
            // printf("New Bipartie edges to choose:\n");
            // for(int i =0 ; i < bg_graph.num_vertex_b ; i++){
            //     printf("index: %d value: %d\n", i, new_vertex_b_debug[i]);
            // }
            // printf("New vertex b num: %d\n", new_bg_graph.num_vertex_b);

            int* prefixSum = NULL;
            cudaMalloc((void**) &(prefixSum), bg_graph.num_vertex_b * sizeof(int));
            prefixCopy<<<((bg_graph.num_vertex_b) + THREADSPERBLOCK-1)/THREADSPERBLOCK, THREADSPERBLOCK>>>(bg_graph.num_vertex_b, new_vertex_b, prefixSum);

            //get index of bipartite edges
            int* d_prefix_helper = NULL;
            cudaMalloc((void**) &(d_prefix_helper), bg_graph.num_vertex_b * sizeof(int));
            /* prefix sum belloch scan */
            int d = 1;
            while(d<bg_graph.num_vertex_b){
                getPrefixSum<<<(bg_graph.num_vertex_b + THREADSPERBLOCK-1)/THREADSPERBLOCK, THREADSPERBLOCK>>>(prefixSum, d_prefix_helper, d);
                d = 2*d;    
            }

            //debugging
            // printf("prefix sum:\n");
            // int* vertex_b_print = (int*) malloc( bg_graph.num_vertex_b * sizeof(int));
            // cudaMemcpy(vertex_b_print, prefixSum,  bg_graph.num_vertex_b * sizeof(int), cudaMemcpyDeviceToHost);
            // for(int i = 0; i <  bg_graph.num_vertex_b ; i++){
            //     printf("vertex: %d index: %d\n", i, vertex_b_print[i]);
            // }

            cudaMalloc((void**) &(new_bg_graph.edges), new_bg_graph.num_bipartite_edges * sizeof(struct b_edge));

            int* d_max_super_vertex = NULL;
            cudaMalloc((void**) &(d_max_super_vertex),sizeof(int));
            get_new_bg_edges<<<(bg_graph.num_vertex_b + THREADSPERBLOCK-1)/THREADSPERBLOCK, THREADSPERBLOCK>>>(bg_graph.num_vertex_b , new_vertex_b, prefixSum, d_vertices_u, d_super_vertices, new_bg_graph.edges, d_max_super_vertex);
            cudaMemcpy(max_super_vertex, d_max_super_vertex, sizeof(int), cudaMemcpyDeviceToHost);

            // debugging
            debugging.num_vertex_a = new_bg_graph.num_vertex_a;
            debugging.num_vertex_b = new_bg_graph.num_vertex_b;
            debugging.num_bipartite_edges = new_bg_graph.num_bipartite_edges ;
            
            debugging.edges = (struct b_edge*) malloc(new_bg_graph.num_bipartite_edges * sizeof(struct b_edge));
            cudaMemcpy(debugging.edges, new_bg_graph.edges, new_bg_graph.num_bipartite_edges * sizeof(struct b_edge), cudaMemcpyDeviceToHost);

            // printf("New Bipartite Graph:\n");
            // printf("verticesA: %d, verticesB: %d, edges: %d\n", debugging.num_vertex_a, debugging.num_vertex_b, debugging.num_bipartite_edges);
            // for(int i = 0; i < debugging.num_bipartite_edges; i++){
            //     printf("index: %d - %d   %d   %d   %d\n", i, debugging.edges[i].v, debugging.edges[i].u, debugging.edges[i].cv, debugging.edges[i].weight);
            // }

            bg_graph.num_vertex_a = new_bg_graph.num_vertex_a;
            bg_graph.num_vertex_b = new_bg_graph.num_vertex_b;
            bg_graph.num_bipartite_edges = new_bg_graph.num_bipartite_edges;

            
            cudaFree(bg_graph.edges);
            cudaMalloc((void**) &(bg_graph.edges), bg_graph.num_bipartite_edges * sizeof(struct b_edge));
            cudaMemcpy(bg_graph.edges, debugging.edges, bg_graph.num_bipartite_edges * sizeof(struct b_edge), cudaMemcpyHostToDevice);

            free(debug_smallest_weights);
            free(debug_smallest_edges);
            free(strut_edges);
            free(vertices_u);
            free(zero_diff_edges);
            free(super_vertices);
            free(new_vertex_b_debug);
            
            cudaFree(d_prefix_helper);
            cudaFree(prefixSum);
            cudaFree(new_vertex_b);
            cudaFree(new_num_vertex_b);
            cudaFree(d_max_super_vertex);
            cudaFree(d_super_vertices);
            cudaFree(d_zero_diff_edges);
            cudaFree(d_vertices_u);
            cudaFree(new_num_vertex_b);
            cudaFree(new_vertex_b);
            cudaFree(d_strut_edges);
            cudaFree(smallest_weights);
            cudaFree(smallest_edges);
        }
    }

    //printf("done with loop\n");
    /*end of while loop*/
    FILE *file;
    file = fopen(argv[argc-1],"w+");
    fprintf(file,"Input Graph\nVertices: %d Edges: %d\n", og_graph.num_vertices, og_graph.num_edges);
    fprintf(file, "MST Edges:\n");
    for(int i = 0; i < og_graph.num_edges; i++){
        if(mst_edges[i] == true){
            fprintf(file, "index: %d - v: %d  u: %d  weight: %d\n", i, og_graph.edges[i].v, og_graph.edges[i].u, og_graph.edges[i].weight);
        }
    }
    fclose(file);

    // malloc frees
    free(max_super_vertex);
    free(og_graph.edges);
    free(mst_edges);
    free(solution_size);
    free(max_super_vertex);
    

    // cuda malloc frees
    cudaFree(d_solutionSize);
    cudaFree(d_mst_edges);
    cudaFree(mst_edges);
	cudaFree(d_og_edges);
	cudaFree(bg_graph.vertices_a);
	cudaFree(bg_graph.vertices_b);
	cudaFree(bg_graph.edges);
}

void get_graph(struct graph* og_graph, char* input){
	FILE *file;
	char buff[255];
	int num_vertices;
	int num_edges;

	file = fopen(input , "r");
    if(file == NULL){
        perror(input);
        exit(1);
    }
    else{
    	fscanf(file, "%s", buff);
    	num_vertices = atoi(buff);

    	fscanf(file, "%s", buff);
    	num_edges = atoi(buff);

    	(*og_graph).num_edges = num_edges;
    	(*og_graph).num_vertices = num_vertices;
    	(*og_graph).edges = (struct edge*) malloc(sizeof(struct edge) * num_edges);

    	for(int i = 0; i < num_edges; i++){
			fscanf(file, "%s", buff);
			(*og_graph).edges[i].v = atoi(buff);

			fscanf(file, "%s", buff);
			(*og_graph).edges[i].u = atoi(buff);

			fscanf(file, "%s", buff);
			(*og_graph).edges[i].weight = atoi(buff);
    	}
    }
    fclose(file);
}


__global__ void get_bipartite_graph(int num_edges, int num_vertices, struct edge* graphEdges, struct b_vertex_a* vertices_a, struct b_vertex_b* vertices_b, struct b_edge* bg_graphEdges) {
    int edge = threadIdx.x + blockIdx.x * blockDim.x;

    if(edge < num_edges){
    	// acquire two bipartite edges for each orginal graph edge
    	bg_graphEdges[2*edge].v = graphEdges[edge].v;
    	bg_graphEdges[2*edge].u = edge;
    	bg_graphEdges[2*edge].cv = 2*edge+1; // corresponding edge/vertex
    	bg_graphEdges[2*edge].weight = graphEdges[edge].weight;

    	bg_graphEdges[2*edge+1].v = graphEdges[edge].u;
    	bg_graphEdges[2*edge+1].u = edge;
    	bg_graphEdges[2*edge+1].cv = 2*edge; // corresponding edge/vertex
    	bg_graphEdges[2*edge+1].weight = graphEdges[edge].weight;

    	vertices_b[edge].e = edge;
    	if(edge < num_vertices)
    		vertices_a[edge].v = edge;
    }
}

__global__ void init_smallest_edges_weights(int num_edges, int *smallest_weights, int* smallest_edges){
    int edge = threadIdx.x + blockIdx.x * blockDim.x;
    if(edge < num_edges){
        smallest_weights[edge] = -1;
        smallest_edges[edge] = -1;
    }
}

// fills in smallest edges array with the index of smallest bipartite edges for each vertex (index of smallest_edges corresponds to vertex number) in graph
__global__ void get_smallest_edges(int bp_num_edges, int num_smallest_edges, struct b_edge* bg_graphEdges, int* smallest_weights, int* smallest_edges){
    int edge = threadIdx.x + blockIdx.x * blockDim.x;
    if(edge< bp_num_edges){
        int index = bg_graphEdges[edge].v - 1;
        // smallest_weights[index] = bg_graphEdges[edge].weight; // filler weight to compare with
        smallest_weights[index] = INT_MAX;
        __syncthreads(); // acquire all smallest weights
        atomicMin(&(smallest_weights[index]), bg_graphEdges[edge].weight); // save actual smallest weight

        __syncthreads(); // acquire all smallest weights

        smallest_edges[index] = bp_num_edges - 1; // filler edge number to comapre with, max edge
        // if(edge < num_smallest_edges){
        //     if(edge != index){
        //         smallest_weights[index] = -1;
        //         smallest_edges[index] = -1;
        //     }
        // }
        __syncthreads(); // acquire all smallest edges

        if(bg_graphEdges[edge].weight == smallest_weights[index]) // save smallest edge if the the bg edge has same weight as smallest weight
			//atomicMin(&(smallest_edges[index]), bg_graphEdges[edge].u);
			atomicMin(&(smallest_edges[index]), edge);
    }
}

// flags all edges to false
__global__ void mst_edges_init(int og_num_edges, bool *mst_edges){
    int edge = threadIdx.x + blockIdx.x * blockDim.x;
    if(edge < og_num_edges){
        mst_edges[edge] = false;
    }
}

// sets which edges go in mst
__global__ void get_mst_edges(int num_smallest_edges, int* smallest_edges, struct b_edge* bg_graphEdges, bool *mst_edges){
    int edge = threadIdx.x + blockIdx.x * blockDim.x;
    int bg_index;
    int vertex;
    if(edge < num_smallest_edges){
        bg_index = smallest_edges[edge];
        if(bg_index != -1){
            vertex = bg_graphEdges[bg_index].u;
            mst_edges[vertex] = true;
        }
    }
}

// gets num of mst edges in solution set
__global__ void get_num_mst(int og_num_edges, bool *mst_edges, int* num_mst){
    int edge = threadIdx.x + blockIdx.x * blockDim.x;
    if(edge < og_num_edges){
        *num_mst = 0; // reset
        __syncthreads();

        if(mst_edges[edge] == true)
            atomicAdd(num_mst, 1);
    }
}

// makes the strut edges
__global__ void get_strut_edges(int bg_num_vertices, int* smallest_edges, struct b_edge* bg_graphEdges, strut_edge* strut_edges){
    int bg_vertex = threadIdx.x + blockIdx.x * blockDim.x;
    
    if(bg_vertex < bg_num_vertices){
        strut_edges[bg_vertex].v = bg_vertex + 1; // vertex
        strut_edges[bg_vertex].u = bg_graphEdges[smallest_edges[bg_vertex]].u; // edge index (u vertex)
        strut_edges[bg_vertex].cv = bg_graphEdges[bg_graphEdges[smallest_edges[bg_vertex]].cv].v; // save vertex that is connected to same edge index (u vertex);
    }
}

// init strut u vertices degree
__global__ void strut_u_init(int bg_num_vertex_b, struct strut_u_vertex* vertices_u){
    int vertex_b = threadIdx.x + blockIdx.x * blockDim.x;
    if(vertex_b < bg_num_vertex_b)
        vertices_u[vertex_b].degree = 0;
}


// fill in degree of strut u vertices
__global__ void get_strut_u_degree(int num_strut_vertices, strut_edge* strut_edges, struct strut_u_vertex* vertices_u){
    int strut_edge = threadIdx.x + blockIdx.x * blockDim.x;

    if(strut_edge < num_strut_vertices){
        atomicAdd(&(vertices_u[strut_edges[strut_edge].u]).degree, 1);
    }
}

// fill in what vertices the vertices_u from the strut is connected
__global__ void get_strut_u_vertices(int bg_num_edges, struct b_edge* bg_graphEdges, struct strut_u_vertex* vertices_u){
    int bg_edge = threadIdx.x + blockIdx.x * blockDim.x;
    if(bg_edge < bg_num_edges){
        if(bg_edge%2 == 0){ // only even edges
            vertices_u[bg_graphEdges[bg_edge].u].v1 = bg_graphEdges[bg_edge].v;
            vertices_u[bg_graphEdges[bg_edge].u].v2 = bg_graphEdges[bg_graphEdges[bg_edge].cv].v;
            vertices_u[bg_graphEdges[bg_edge].u].weight =  bg_graphEdges[bg_edge].weight;
        }
    }

}

// get number of zero difference vertrices u in strut
__global__ void get_zero_diff_num(int bg_num_vertex_b, struct strut_u_vertex* vertices_u, int* zero_diff_edges){
    int vertex_b = threadIdx.x + blockIdx.x * blockDim.x;
    if(vertex_b < bg_num_vertex_b){
        if(vertices_u[vertex_b].degree == 2)
            atomicAdd(zero_diff_edges, 1);
    }
}

// initialize super vertices
__global__ void super_vertices_init(int num_strut_vertices, int* super_vertices){
    int vertex = threadIdx.x + blockIdx.x * blockDim.x; 
    if(vertex < num_strut_vertices){
        super_vertices[vertex] = vertex + 1;
    }
}

// set which verticies_u will be in new bipartitie graph and get how many there are 
__global__ void get_new_bg_vertex_b(int num_bg_vertexb, int* super_vertices, struct strut_u_vertex* vertices_u, int* new_vertex_b, int* num_newbg_vertexb){
    int vertex = threadIdx.x + blockIdx.x * blockDim.x; 
    if(vertex < num_bg_vertexb){
        new_vertex_b[vertex] = 0; // setting all to false
        __syncthreads();
        if(super_vertices[vertices_u[vertex].v1 - 1] != super_vertices[vertices_u[vertex].v2 - 1]){
            new_vertex_b[vertex] = 1; 
        }
        *num_newbg_vertexb = 0;
        __syncthreads();
        if(new_vertex_b[vertex] == 1)
            atomicAdd(num_newbg_vertexb, 1);
    }
}

// copy maker
__global__ void prefixCopy(int prefixNum, int* old_prefix, int *new_prefix){
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    if(index < prefixNum){
        new_prefix[index] = old_prefix[index];
    }
}
__global__ void getPrefixSum(int* entries, int* entriesC, int d) {
    int index = threadIdx.x + blockIdx.x * blockDim.x;

    if(index >= d)
        entriesC[index] = entries[index - d];
    else
        entriesC[index] = 0;

    __syncthreads();

    entries[index] = entries[index] + entriesC[index];
}

// makes new bipartite edges
__global__ void get_new_bg_edges(int num_bg_vertex_b, int* new_bg_edges, int* prefixSum, struct strut_u_vertex* vertices_u, int* super_vertices, struct b_edge* bg_graphEdges, int * max_super_vertex){
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    int edge1;
    int edge2;

    if(index < num_bg_vertex_b){
        if(new_bg_edges[index] == 1){
            edge1 = (prefixSum[index] - 1) * 2;
            edge2 = edge1+1;

            bg_graphEdges[edge1].v = super_vertices[vertices_u[index].v1-1];
            bg_graphEdges[edge1].u = index;
            bg_graphEdges[edge1].cv = edge2;
            bg_graphEdges[edge1].weight = vertices_u[index].weight;

            bg_graphEdges[edge2].v = super_vertices[vertices_u[index].v2-1];
            bg_graphEdges[edge2].u = index;
            bg_graphEdges[edge2].cv = edge1;
            bg_graphEdges[edge2].weight = vertices_u[index].weight;

            atomicMax(max_super_vertex, super_vertices[vertices_u[index].v1-1]);
            atomicMax(max_super_vertex, super_vertices[vertices_u[index].v2-1]);
        }
    }
}

// get what vertex each vertex is compacted to during compression of bipartite graph
__global__ void get_super_vertices(int num_strut_vertices, strut_edge* strut_edges, struct strut_u_vertex* vertices_u, int* super_vertices){
	int strut_edge = threadIdx.x + blockIdx.x * blockDim.x; 
    int cv;
    int min_cv;
    if(strut_edge < num_strut_vertices){
        min_cv = strut_edges[strut_edge].v;
        // cv = strut_edges[strut_edge].cv;
        // while(cv < min_cv){
		// 	min_cv = cv;
        //     cv = strut_edges[cv-1].cv;
        // }
        if(vertices_u[strut_edges[strut_edge].u].v1 == min_cv) 
            cv = vertices_u[strut_edges[strut_edge].u].v2;
        else
            cv = vertices_u[strut_edges[strut_edge].u].v1;
        while(cv < min_cv){
            min_cv = cv;
            if(vertices_u[strut_edges[strut_edge].u].v1 == min_cv)
                cv = vertices_u[strut_edges[cv-1].u].v2;
            else
                cv = vertices_u[strut_edges[cv-1].u].v1;
        }
        super_vertices[strut_edge] = min_cv;
    }
}