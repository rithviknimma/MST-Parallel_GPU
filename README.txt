Parallel computation of minimum spanning trees using GPU in CUDA C/C++

/*****************************************************************/

Compile with:
nvcc -o mst.out mst.cu

To run:
mst.out <Input file> <Output file>

/*****************************************************************/

- Input file has to be a text file that contains the graph

- Provided some sample input files

- Output will be all the MST edges written to the output file name specified by the arguments

/*****************************************************************/

This algorithm creates bipartite graphs and uses a decomposition called "strut" to find MST edges
in O(log p) iterations, where p is for the processors. 

This algorithm was found in a research paper that we have also put in this folder titled "RESEARCH_PAPER"

It contains the detailed description for the algorithm we have implemented.

/*****************************************************************/

This is the pseudocode for the algorithm:

Input: A connected graph G = (V, E), where V = {v1, v2, . . . , vn} is a set of n vertices and E is a
set of m edges (vi, vj ), where vi and vj are vertices of V.

Each edge (vi, vj ) has a weight denoted by wij .
Output: A minimum spanning tree of G whose edges are in SolutionEdgeSet.

1: SolutionEdgeSet := empty.
2: Transform the given graph G into a bipartite graph H =(V, U, E0).
3: condition := true.
4: while condition do
5:   Obtain a strut.
6: 	for every strut-edge (vi, uj ) do
7: 		Add the edge (vi, vk) to the SolutionEdgeSet, where vi and vk are adjacent to uj . In other words,
		add original edge(uj ) to the SolutionEdgeSet.
8: 	end for
9:   Compute the number of zero-difference vertices.
10:   if number of zero-difference vertices = 1 then
11: 	condition := false
12:   else
13: 	for every strut-edge (vi, uj ) do
14: 	  Compact the two vertices adjacent to uj therebyproducing a new bipartite graph.
15: 	end for
16:   end if
17: end while


/*****************************************************************/