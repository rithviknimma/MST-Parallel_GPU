/*               A Graph Generation Package
            Richard Johnsonbaugh and Martin Kalin
   Department of Computer Science and Information Systems
                      DePaul University
                     Chicago, IL  60604
          johnsonbaugh@cs.depaul.edu, kalin@cs.depaul.edu

     This program generates graphs of specified  sizes  and
properties. Among the types of graphs are

     * random graphs
     * random connected graphs
     * random directed acyclic graphs
     * random complete weighted graphs
     * random pairs of isomorphic regular graphs
     * random graphs with Hamiltonian cycles
     * random networks

Graphs  may be specified further with respect to one or more
of these properties:

     * weighted or unweighted
     * directed or undirected
     * simple or nonsimple

We donate the program freely and hope others find it useful. We only
request attribution whenever this program is (successfully) shared. In
the spirit of the million dollar commercial vendors, we assume absolutely
no responsibility whatsoever for whatever mayhem this program
may inflict on your own system. Good luck!

     The program is menu-driven. The top-level menu is

     0. Exit Program
     1. Random Graph
     2. Random Digraph
     3. Random Connected Graph
     4. Random DAG (Directed Acyclic Graph)
     5. Random Complete Weighted Graph
     6. Random Isomorphic Graph Pair
     7. Random Hamiltonian Cycle Graph
     8. Random Network
     Choice?

Selecting  any  choice  except  0  brings  a  sub-menu.  For
example,  if  the  user  selects  1,  a  series  of  prompts
results:

          Number of Vertices:

          Number of Edges:

          1. Not Simple
          2. Simple
          Choice?

          1. Unweighted Edges
          2. Weighted Edges
          Choice?

          File for graph:

A  nonsimple  graph  is  allowed  to  contain parallel edges
and/or  loops.  If  the  user  selects  the  Weighted  Edges
option, the prompt

     Maximum Edge Weight:

allows  the  user  to  enter a maximum edge weight. A random
graph  meeting  the  chosen specifications is then generated
and  written to the file named by the user. Output takes the
form

n     m
v1    w1    c1
v2    w2    c2
      .
      .
      .

where  n  is  the  number of vertices and m is the number of
edges  in  the  graph. Edge i is incident on vertices vi and
wi.  The weight of edge i is ci. If the graph is unweighted,
the value of ci is 1 for all i.

     Other  options work similarly. For example, if the user
selects 2, Random Digraph, an output line

     vi    wi    ci

denotes  an  edge  directed  from  vi  to wi with weight ci.
Random  graphs  and digraphs are useful as input to programs
that  find components, that test for connectedness, and that
find   spanning   forests. Programs  that  find
spanning  trees  or minimal spanning trees can use
as  input  random  connected graphs, provided as option 3 on
the  top-level  menu.  Option 3 allows the user to designate
the connected graph's edges as weighted or unweighted.

     Directed  acyclic graphs, suitable as input to programs
that  find  topological orderings, are provided by
option  4  in  the  top-level  menu.  The program writes the
generated  graph to one file and and a topological sort to a
second  file.  The  user  may specify the number of vertices
and  edges,  the  maximum  weight per edge, and the two file
names.

     Complete   weighted  graphs  are  useful  as  input  to
programs   that   either   solve   exactly   the   traveling
salesperson  problem  or  approximate  an  exact solution to
this  problem. Option  5 generates such graphs,
allowing  the  user  to  specify  the  number  of  vertices,
whether  the graph is directed, the maximum weight per edge,
and the file name.

     Pairs  of  isomorphic  graphs,  suitable  as  input  to
programs  that  test  for  graph  isomorphism, are
provided  by  option  6  in  the top-level menu. The program
writes  two  isomorphic,  regular,  undirected graphs (i.e.,
isomorphic,  undirected  graphs in which all of the vertices
have   the   same   degree)  to  two  files  and  writes  an
isomorphism  to  a  third  file. Pairs of isomorphic regular
graphs  make  challenging  input  to  programs that test for
graph  isomorphism. Under this
option,  the  user  specifies  the  number  of vertices, the
common vertex degree, and the file names.

     Graphs  with  Hamiltonian cycles are useful as input to
programs  that  find  such  cycles or determine that they do
not  exist. Option  7  on  the  top-level  menu
provides  graphs with Hamiltonian cycles. The program writes
the  graph  to  one  file and a Hamiltonian cycle to another
file.  The  user may choose the number of vertices and edges
as well as the file names for the graph and the cycle.

     Networks,  provided  as option 8 on the top-level menu,
are  suitable  as  input to programs that find maximal flows.
Under  this  option,  the program generates a
simple, weighted, directed graph satisfying the following:

(a)  Vertex 1 has no incoming edges.

(b)  Vertex  n,  where n is equal to the number of vertices,
     has no outgoing edges.

The  user  selects  the  number  of  vertices and edges, the
maximum weight per edge, and the file name.

     This program does some error and integrity checking. As
examples,  the  user  may  not  specify  a  graph  with zero
vertices  or  edges.  If  a  simple, undirected graph with n
vertices  is  to  be  generated,  the user should specify at
most  n(n - 1)/2  edges. If a Hamiltonian cycle graph with n
vertices  is  to  be  generated,  the user should specify at
least  n  edges.  If  a  network  with  n  vertices is to be
generated,  the  user  should  specify  at  most n*n - 3n + 3
edges.  The program corrects erroneous user input wherever
feasible.  In the examples just cited, the program would set
the edge count correctly if the user did not.

   For a full discussion, see the article by Johnsonbaugh
and Kalin, "A Graph Generation Software Package,"
in Proceedings of the 22nd SIGCSE Technical Symposium.
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

/*** record for qsorting a graph to provide another isomorphic to it ***/
#define RecordSize 24
typedef struct {
   char record[ RecordSize ];
} Record;

int zero_vertex_or_edge_count( void );
void fix_imbalanced_graph( void );
void print_graph( int v,
                  int e,
                  char* out_file,
                  int* adj_matrix,
                  int dir_flag );
int unpicked_vertices_p( int* closed, int n );
int pick_to( int* closed, int vertices );
void network_flow_graph( int vertices,
                         int edges,
                         int max_weight,
                         char* outfile );
void hamiltonian_cycle_graph( int v,
                              int e,
                              int dir_flag,
                              char* out_file,
                              char* ham_file );
void complete_graph( int v,
                     int max_wgt,
                     int dir_flag,
                     char* out_file );
void directed_acyclic_graph( int v,
                             int e,
                             int max_wgt,
                             int weight_flag,
                             char* out_file,
                             char* dag_file );
void random_digraph( int v,
                     int e,
                     int max_wgt,
                     int weight_flag,
                     int simple_flag,
                     char* out_file );
void r_graph( int v,
              int e,
              int max_wgt, int weight_flag,
              char* out_file );
void random_graph( int v,
                   int e,
                   int max_wgt,
                   int weight_flag,
                   int simple_flag,
                   char* out_file );
void random_connected_graph( int v,
                             int e,
                             int max_wgt,
                             int weight_flag,
                             char* out_file );
void rename_edges( Record* ptr, int* aliases, int count );
void rename_vertices( int* aliases, int count );
int compare_strings( const void* a, const void* b );
void isomorphic_graph_pair( int vertices,
                            int spec_degree,
                            char* first,
                            char* second,
                            char* third );
void warning( char* message );
int bogus_file_name( int count );
void build_graphs( void );
void build_next_graph( void );
void display_main_menu( void );
void process_choice( void );
void display_simple_or_not_menu( void );
void display_dimensions_menu( int which );
void display_weighted_menu( void );
void display_max_weight_menu( void );
void display_directed_menu( void );
void display_max_degree_menu( void );
void display_outfile_menu( int count, int index1, int index2 );
void seed_ran( void );
int illegal_parms( int files );
int ran( int k );     /* customized random number generator */
void permute( int* a, int n ); /* gives a random permutation */

/* initialize array a so that a[ i ] = i */
void init_array( int* a, int end );

void swap( int* a, int *b ); /* swap two ints */
int get_int( char* indent ); /* get user input & check for valid input */
/*** miscellany ***/
#define True         1
#define False        0
#define None         -999
#define Undirected   0
#define Directed     1
#define ExitProgram  0

/*** graph generators ***/
#define RandomGraph            1
#define RandomDigraph          2
#define RandomConnectedGraph   3
#define DirectedAcyclicGraph   4
#define CompleteGraph          5
#define IsomorphicGraphPair    6
#define HamiltonianCycleGraph  7
#define NetworkFlowGraph       8

/*** warnings ***/
#define IllegalDimensions   "\n\n\tNeither edge nor vertex count may be zero.\n\n"
#define IllegalFileName     "\n\n\tOutput file has illegal name.\n\n"
#define InsufficientStorage "\n\n\tThere is not enough space for this size graph.\n\n"
#define SingleVertexNetwork "\n\n\tA single vertex network is not allowed.\n\n"

/*** prompts ***/
#define ToplevelPrompts  9
static char  *toplevel_prompts[] =
   { " 0. Exit Program ",
     " 1. Random Graph ",
     " 2. Random Digraph ",
     " 3. Random Connected Graph ",
     " 4. Random DAG (Directed Acyclic Graph) ",
     " 5. Random Complete Weighted Graph ",
     " 6. Random Isomorphic Graph Pair ",
     " 7. Random Hamiltonian Cycle Graph ",
     " 8. Random Network "
   };

#define SimpleGraphAns  1
#define SimpleOrNotPrompts  2
static char *simple_or_not_prompts[] =
   { " 1. Simple ",
     " 2. Not Simple "
   };

#define WeightedAns  2
#define WeightedOrNotPrompts  2
static char *weighted_or_not_prompts[] =
   { " 1. Unweighted Edges ",
     " 2. Weighted Edges "
   };

#define EdgeIndex         0
#define VertexIndex       1
#define VerticesOnly      -111
#define EdgesAndVertices  -222
static char *dimension_prompts[] =
   { " Number of Edges:    ",
     " Number of Vertices: "
   };

#define DirectedGraphAns 1
#define DirectedPrompts 2
static char *directed_prompts[] =
   { " 1. Directed Graph ",
     " 2. Undirected Graph "
   };

static char *degree_prompt =  " Degree Per Edge:  ";

#define ConnectedAns  1
#define ConnectedPrompts  2
static char *connected_prompts[] =
   { " 1. Connected Graph ",
     " 2. Unconnected Graph "
   };

static char *main_file_prompt = " File for graph:  ";

#define TopologicalSort      0
#define IsomorphicGraph      1
#define Isomorphism          2
#define HamiltonianCycle     3
static char *file_prompts[] =
   { " File for topological sort:  ",
     " File for isomorphic graph:  ",
     " File for isomorphism:  ",
     " File for Hamiltonian cycle:  "
   };


/*** output files ***/
#define MaxFileName   128
typedef struct out_files {
   char  outfile1[ MaxFileName + 1 ];
   char  outfile2[ MaxFileName + 1 ];
   char  outfile3[ MaxFileName + 1 ];
} Outfile;

/*** graphs ***/
typedef struct graph_parms {
   int      edge_count;
   int      vertex_count;
   int      max_weight;
   int      max_degree;
} GraphParms;

typedef struct graph_props {
   int      simple_p;
   int      weighted_p;
   int      directed_p;
   int      dag_p;
   int      isomorphic_p;
   int      network_p;
} GraphProps;

/*** menu ***/
typedef struct menu {
   GraphProps   props;
   GraphParms   parms;
   Outfile      outfiles;
} MenuStruct;

typedef MenuStruct *Menu;

static MenuStruct  menu_structure;
static Menu        menu = &menu_structure;

static int   build_more_graphs,
             menu_choice;

/* min, max, odd */
#undef min
#define min( x, y )   ((( x ) < ( y )) ? ( x ) : ( y ))

#undef max
#define max( x, y )   ((( x ) > ( y )) ? ( x ) : ( y ))

#define odd( num )    ( ( num ) % 2 )


int main()
{
   build_more_graphs = True;
   build_graphs();
   return 0;
}

void build_graphs( void )
{
   while ( build_more_graphs )
      build_next_graph();
}

void build_next_graph( void )
{
   display_main_menu();
   process_choice();
}


#define ChoicePrompt  " Choice? "
void display_main_menu( void )
{
   int  i;

   printf( "\n\n" );
   for ( i = 0; i < ToplevelPrompts; i++ )
      printf( "\n\t%s", toplevel_prompts[ i ] );
   printf( "\n\t" );
   printf( ChoicePrompt ); menu_choice = get_int( "\n\t" );
}

void process_choice( void )
{

   seed_ran();  /* seed system's random number generator */

   menu -> props.weighted_p =
      menu -> props.dag_p =
         menu -> props.isomorphic_p =
            menu -> props.simple_p =
               menu -> props.directed_p =
                  menu -> props.network_p = False;

   switch ( menu_choice ) {
      case ExitProgram:
         build_more_graphs = False;
         return;
      case RandomGraph:
         display_dimensions_menu( EdgesAndVertices );
         display_simple_or_not_menu();
         display_weighted_menu();
         if ( menu -> props.weighted_p )
            display_max_weight_menu();
         display_outfile_menu( 1, None, None );
         if ( illegal_parms( 1 ) )
            return;
         else
            random_graph( menu -> parms.vertex_count,
                          menu -> parms.edge_count,
                          menu -> parms.max_weight,
                          menu -> props.weighted_p,
                          menu -> props.simple_p,
                          menu -> outfiles.outfile1 );
      break;
      case RandomDigraph:
         menu -> props.directed_p = True;
         display_dimensions_menu( EdgesAndVertices );
         display_simple_or_not_menu();
         display_weighted_menu();
         if ( menu -> props.weighted_p )
            display_max_weight_menu();
         display_outfile_menu( 1, None, None );
         if ( illegal_parms( 1 ) )
            return;
         else
            random_digraph( menu -> parms.vertex_count,
                            menu -> parms.edge_count,
                            menu -> parms.max_weight,
                            menu -> props.weighted_p,
                            menu -> props.simple_p,
                            menu -> outfiles.outfile1 );
      break;
      case RandomConnectedGraph:
         menu -> props.simple_p = True;
         display_dimensions_menu( EdgesAndVertices );
         display_weighted_menu();
         if ( menu -> props.weighted_p )
            display_max_weight_menu();
         display_outfile_menu( 1, None, None );
         if ( illegal_parms( 1 ) )
            return;
         else
            random_connected_graph( menu -> parms.vertex_count,
                                    menu -> parms.edge_count,
                                    menu -> parms.max_weight,
                                    menu -> props.weighted_p,
                                    menu -> outfiles.outfile1 );
      break;
      case DirectedAcyclicGraph:
         menu -> props.dag_p = True;
         display_dimensions_menu( EdgesAndVertices );
         display_weighted_menu();
         if ( menu -> props.weighted_p )
            display_max_weight_menu();
         display_outfile_menu( 2, TopologicalSort, None );
         if ( illegal_parms( 2 ) )
            return;
         else
            directed_acyclic_graph( menu -> parms.vertex_count,
                                    menu -> parms.edge_count,
                                    menu -> parms.max_weight,
                                    menu -> props.weighted_p,
                                    menu -> outfiles.outfile1,
                                    menu -> outfiles.outfile2 );
      break;
      case CompleteGraph:
         display_dimensions_menu( VerticesOnly );
         display_max_weight_menu();
         display_directed_menu();
         display_outfile_menu( 1, None, None );
         complete_graph( menu -> parms.vertex_count,
                         menu -> parms.max_weight,
                         menu -> props.directed_p,
                         menu -> outfiles.outfile1 );
      break;
      case IsomorphicGraphPair:
         menu -> props.isomorphic_p = True;
         display_dimensions_menu( VerticesOnly );
         display_max_degree_menu();
         display_outfile_menu( 3, IsomorphicGraph, Isomorphism );
         if ( illegal_parms( 3 ) )
            return;
         else
            isomorphic_graph_pair( menu -> parms.vertex_count,
                                   menu -> parms.max_degree,
                                   menu -> outfiles.outfile1,
                                   menu -> outfiles.outfile2,
                                   menu -> outfiles.outfile3 );
      break;
      case HamiltonianCycleGraph:
         menu -> props.simple_p = True;
         display_dimensions_menu( EdgesAndVertices );
         display_directed_menu();
         display_outfile_menu( 2, HamiltonianCycle, None );
         if ( illegal_parms( 2 ) )
            return;
         else
            hamiltonian_cycle_graph( menu -> parms.vertex_count,
                                     menu -> parms.edge_count,
                                     menu -> props.directed_p,
                                     menu -> outfiles.outfile1,
                                     menu -> outfiles.outfile2 );
      break;
      case NetworkFlowGraph:
         menu -> props.network_p = True;
         menu -> props.directed_p = True;
         display_dimensions_menu( EdgesAndVertices );
         display_max_weight_menu();
         display_outfile_menu( 1, None, None );
         if ( illegal_parms( 1 ) )
            return;
         else
            network_flow_graph( menu -> parms.vertex_count,
                                menu -> parms.edge_count,
                                menu -> parms.max_weight,
                                menu -> outfiles.outfile1 );
      break;
   }
}

void display_simple_or_not_menu( void )
{
   int i;

   for ( i = 0; i < SimpleOrNotPrompts; i++ )
      printf( "\n\t\t%s", simple_or_not_prompts[ i ] );
   printf( "\n\t\t" );
   printf( ChoicePrompt );
   if ( get_int( "\n\t\t" ) == SimpleGraphAns )
      menu -> props.simple_p = True;
   else
      menu -> props.simple_p = False;
}

void display_directed_menu( void )
{
   int i;

   for ( i = 0; i < DirectedPrompts; i++ )
      printf( "\n\t\t%s", directed_prompts[ i ] );
   printf( "\n\t\t" );
   printf( ChoicePrompt );
   if ( get_int( "\n\t\t" ) == DirectedGraphAns )
      menu -> props.directed_p = True;
   else
      menu -> props.directed_p = False;
}

void display_dimensions_menu( int which )
{
   printf( "\n\t\t%s", dimension_prompts[ VertexIndex ] );
   menu -> parms.vertex_count = get_int( "\n\t\t" );
   if ( which == EdgesAndVertices ) {
      printf( "\t\t%s", dimension_prompts[ EdgeIndex ] );
      menu -> parms.edge_count = get_int( "\n\t\t" );
   }
}

void display_max_degree_menu( void )
{
   int  degree;

   printf( "\n\t\t%s", degree_prompt );
   degree = get_int( "\n\t\t" );
   menu -> parms.max_degree = degree;
   menu -> parms.edge_count = 1;       /* to avoid error message */
}

#define DefaultMaxWeight 1
void display_weighted_menu( void )
{
   int  i;

   for ( i = 0; i < WeightedOrNotPrompts; i++ )
      printf( "\n\t\t%s", weighted_or_not_prompts[ i ] );
   printf( "\n\t\t" );
   printf( ChoicePrompt );
   if ( get_int( "\n\t\t" ) == WeightedAns )
      menu -> props.weighted_p = True;
   else {
      menu -> props.weighted_p = False;
      menu -> parms.max_weight = DefaultMaxWeight;
   }
}

void display_max_weight_menu( void )
{
   printf( "\n\t\t\tMaximum Edge Weight:  " );
   menu -> parms.max_weight = get_int( "\n\t\t\t" );
}

void display_outfile_menu( int count, int index1, int index2 )
{
   int  i;

   for ( i = 0; i < count; i++ ) {
      if ( i == 0 ) {
         printf( "\n\t\t%s", main_file_prompt );
         scanf( "%s", menu -> outfiles.outfile1 );
      }
      else if ( i == 1 ) {
         printf( "\n\t\t%s", file_prompts[ index1 ] );
         scanf( "%s", menu -> outfiles.outfile2 );
      }
      else {
         printf( "\n\t\t%s", file_prompts[ index2 ] );
         scanf( "%s", menu -> outfiles.outfile3 );
      }
   }
}

int zero_vertex_or_edge_count( void )
{
   return ( menu -> parms.vertex_count == 0  ||
            menu -> parms.edge_count == 0 );
}

int bogus_file_name( int count )
{
   if ( count == 1 )
      return ( !strlen( menu -> outfiles.outfile1 ) );
   else
      return ( !strlen( menu -> outfiles.outfile1 ) &&
               !strlen( menu -> outfiles.outfile2 ) );
}

void fix_imbalanced_graph( void )
{
   int  max_edges;

   if ( menu -> props.simple_p ) {
      max_edges = menu -> parms.vertex_count
                  * ( menu -> parms.vertex_count - 1 );
      if ( !menu -> props.directed_p )
         max_edges /= 2;
      if ( menu -> parms.edge_count > max_edges )
         menu -> parms.edge_count = max_edges;
   }
   else if ( menu -> props.dag_p ) {
      max_edges = ( menu -> parms.vertex_count
                  * ( menu -> parms.vertex_count - 1 ) ) / 2;
      if ( menu -> parms.edge_count > max_edges )
         menu -> parms.edge_count = max_edges;
   }
   else if ( menu -> props.isomorphic_p ) {
      if ( odd( menu -> parms.max_degree ) )
         if ( odd( menu -> parms.vertex_count ) )
            menu -> parms.vertex_count++;
         if ( menu -> parms.vertex_count <= menu -> parms.max_degree )
            menu -> parms.vertex_count = menu -> parms.max_degree + 1;
   }
   else if ( menu -> props.network_p ) {
      max_edges =
         ( menu -> parms.vertex_count * menu -> parms.vertex_count ) -
         ( menu -> parms.vertex_count * 3 ) + 3;
      if ( menu -> parms.edge_count >= max_edges )
         menu -> parms.edge_count = max_edges;
   }
}

int illegal_parms( int files )
{
   if ( zero_vertex_or_edge_count() ) {
      warning( IllegalDimensions );
      return ( True );
   }
   else if ( bogus_file_name( files ) ) {
      warning( IllegalFileName );
      return ( True );
   }
   else {
      fix_imbalanced_graph();
      return ( False );
   }
}

void warning( char* message )
{
   printf("%s",message );
}


/******************** graph generators **************************/

/* This function writes a pair of isomorphic, regular,
   undirected graphs to two files and writes an isomorphism to a third
   file.

   To  generate  a  pair  of isomorphic simple and regular
   graphs,  we  first  generate  a  random  simple  and regular
   graph.  We then construct an isomorphic copy by renaming its
   vertices using a random permutation.

   To generate a random simple graph
   in   which  every  vertex  has  degree  r, we  execute  the
   following,  which  we  subsequently  refer  to  as  the main
   procedure:

     1.   Choose  vertices  v and w each of degree less than
          r  such  that there is no edge between v and w. If
          no such pair exists, execute the fixup procedure.

     2.   Add an edge between v and w.

     3.   If some vertex has degree less than r, go to 1.

   The  fixup  procedure  works as follows. There are two cases
   to consider:

     1.   All vertices except one have degree equal to r.

     2.   There  are  pairs  of vertices each of degree less
          than  r,  but  each  such  pair is connected by an
          edge.

   In  the  first case, we increase the edge count in the graph
   in  the following way. Suppose that vertex v has degree less
   than  r.  Choose an edge (w1,w2) such that (v,w1) and (v,w2)
   are  not edges. Delete the edge (w1,w2) and add edges (v,w1)
   and  (v,w2).  In the second case, we increase the edge count
   in  the graph in a similar way. Suppose that vertices v1 and
   v2  have  degree  less than r and (v1,v2) is an edge. Choose
   an  edge  (w1,w2)  such  that  (v1,w1)  and  (v2,w2) are not
   edges.  Delete  the  edge  (w1,w2) and add edges (v1,w1) and
   (v2,w2).  If  some  vertex in the resulting graph has degree
   less than r, we resume by executing the main procedure.
 */
void isomorphic_graph_pair( int vertices,
                            int spec_degree,
                            char* first,
                            char* second,
                            char* third )
{
   Record  *records, *ptr;
   int     i, j,
           size = vertices * vertices,
           edge_count = 0,
           edges = ( vertices * spec_degree ) / 2,
           base, offset, index1, index2, index3,
           v, v1, v2,
           w, w1, w2,
           *adj_matrix,
           *vertex_degree,
           *aux,
           *aliases;
   FILE    *fptr;

   /* adjacency matrix to represent graph */
   if ( ( adj_matrix = ( int * ) calloc( size, sizeof( int ) ) ) == NULL ) {
      warning( InsufficientStorage );
      return;
   }

   /* vectors to hold degree count per vertex */
   if ( ( vertex_degree = ( int * ) calloc( vertices, sizeof( int ) ) ) == NULL ) {
      warning( InsufficientStorage );
      free( adj_matrix );
      return;
   }
   if ( ( aux = ( int * ) calloc( vertices, sizeof( int ) ) ) == NULL ) {
      warning( InsufficientStorage );
      free( adj_matrix );
      free( vertex_degree );
      return;
   }

   /* vector to hold renamed vertices */
   if ( ( aliases = ( int * ) calloc( vertices, sizeof( int ) ) ) == NULL ) {
      warning( InsufficientStorage );
      free( adj_matrix );
      free( vertex_degree );
      free( aux );
      return;
   }
   
   /* file entries to be sorted to get isomorphic graph */
   if ( ( records = ( Record * ) calloc( edges, sizeof( Record ) ) ) == NULL ) {
      warning( InsufficientStorage );
      free( adj_matrix );
      free( vertex_degree );
      free( aux );
      free( aliases );
      return;
   }
   else
      ptr = records;
   
   /* Initialize adjacency matrix for random permutation. */
   init_array( aux, vertices );

   /* loop until every vertex has spec_degree */
   while ( edge_count < edges ) {
      permute( aux, vertices );
      w = -1;
      /* Look for 2 vertices whose degree is less than r.
         If found & no edge, put in an edge & repeat.
         If we can find 2 vertices whose degree 
         is less than r, but each such pair is joined by an edge,
         set w != -1 and record one such pair as v1, w1.
         If there is no pair of vertices whose degree is less than r,
         w will still be -1 and there will be
         exactly one vertex v whose degree is less than r.  */
      for ( i = 0; i < vertices; i++ ) {
         if ( vertex_degree[ aux[ i ] ] < spec_degree ) {
            v = aux[ i ];
            for ( j = i + 1; j < vertices; j++ ) {
               if ( vertex_degree[ aux[ j ] ] < spec_degree ) {
                  w = aux[ j ];
                  base = min( v, w );
                  offset = max( v, w );
                  index1 = ( base * vertices ) + offset;
                  if ( !adj_matrix[ index1 ] ) { /* no edge so add one */
                     adj_matrix[ index1 ] = 1;
                     vertex_degree[ v ]++;
                     vertex_degree[ w ]++;
                     goto next;
                  }
                  else { /* save the vertices where there is an edge */
                     v1 = v;
                     w1 = w;
                  }
               }
            }
         }
      }
      if ( w == -1 ) { /* only one vertex v with degree < spec_degree */
         v1 = v;
         /* Find an edge (v2,w2) such that

               v2 != v1 && w2 != v1 && there are no edges (v1,v2)
               and (v1,w2).

            Add edges (v1,v2) and (v1,w2) and delete edge (v2,w2).  */
         for ( i = 0; i < vertices - 1; i++ ) {
            if ( v1 != ( v2 = aux[ i ] ) ) {
               base = min( v1, v2 );
               offset = max( v1, v2 );
               index1 = ( base * vertices ) + offset;
               if ( !adj_matrix[ index1 ] ) {
                  for ( j = i + 1; j < vertices; j++ ) {
                     if ( ( w2 = aux[ j ] ) != v1 ) {
                        base = min( v2, w2 );
                        offset = max( v2, w2 );
                        index2 = ( base * vertices ) + offset;
                        base = min( v1, w2 );
                        offset = max( v1, w2 );
                        index3 = ( base * vertices ) + offset;
                        if ( adj_matrix[ index2 ]  &&  !adj_matrix[ index3 ] ) {
                           adj_matrix[ index1 ] =
                              adj_matrix[ index3 ] = 1;
                           adj_matrix[ index2 ] = 0;
                           vertex_degree[ v1 ]++;
                           vertex_degree[ v1 ]++;
                           goto next;
                        }
                     }
                  }
               }
            }
         }
      }
      else {  /* two vertices v1 and w1 with degree < spec_degree */
         /* Find an edge (v2,w2) such that

               v2 != v1 && w1 != v2 && w2 != v1 && w2 != w1
               and there are no edges (v1,v2) and (w1,w2).

            Add edges (v1,v2) and (w1,w2) and delete edge (v2,w2).

            [Note that both loops must run from 0 to n-1 since
            there is an asymmetry in the conditions. For example,
            there might be an edge (v2,w2) such that there is an
            edge (v1,v2) yet there are no edges (v1,w2) and (w1,v2).]      */
         for ( i = 0; i < vertices; i++ ) {
            if ( v1 != ( v2 = aux[ i ] ) && v2 != w1 ) {
               base = min( v1, v2 );
               offset = max( v1, v2 );
               index1 = ( base * vertices ) + offset;
               if ( !adj_matrix[ index1 ] ) {
                  for ( j = 0; j < vertices; j++ ) {
                     if ( ( w2 = aux[ j ] ) != v1 && w2 != w1 ) {
                        base = min( v2, w2 );
                        offset = max( v2, w2 );
                        index2 = ( base * vertices ) + offset;
                        base = min( w1, w2 );
                        offset = max( w1, w2 );
                        index3 = ( base * vertices ) + offset;
                        if ( adj_matrix[ index2 ]  &&  !adj_matrix[ index3 ] ) {
                           adj_matrix[ index2 ] = 0;
                           adj_matrix[ index1 ] =
                              adj_matrix[ index3 ] = 1;
                           vertex_degree[ v1 ]++;
                           vertex_degree[ w1 ]++;
                           goto next;
                        }
                     }
                  }
               }
            }
         }
      }
      next: edge_count++;
   }

   /* print original graph */
   print_graph( vertices, edges, first, adj_matrix, Undirected );
   
   /* print an isomorphism */
   rename_vertices( aliases, vertices );
   if ( ( fptr = fopen( third, "w" ) ) == NULL ) {
      printf( "\n\t\tCould not open file %s.\n\n", third );
      free( adj_matrix );
      free( vertex_degree );
      free( aux );
      free( aliases );
      free( records );
      return;
   }
   else 
      for ( i = 0; i < vertices; i++ )
         fprintf( fptr, "%5d  %5d\n", i + 1, aliases[ i ] + 1 );
   fclose( fptr );

   /* print a graph isomorphic to the original:
        -- change edge names to aliases
        -- sort file to help disguise isomorphism
   */
   if ( ( fptr = fopen( first, "r" ) ) == NULL ) {
      printf( "\n\t\tCould not open file %s.\n\n", first );
      free( adj_matrix );
      free( vertex_degree );
      free( aux );
      free( aliases );
      free( records );
      return;
   }
   else {
      fgets( ptr -> record, RecordSize, fptr ); /* 1st record, vertex and edge count, is ignored */
      while ( ( fgets( ptr -> record, RecordSize, fptr ) ) != NULL )
         ptr++;
      fclose( fptr );
      rename_edges( records, aliases, edges );
      qsort( ( void * ) records, edges, sizeof( Record ), compare_strings );
      if ( ( fptr = fopen( second, "w" ) ) == NULL ) {
         printf( "\n\t\tCould not open file %s.\n", second );
         free( adj_matrix );
         free( vertex_degree );
         free( aux );
         free( aliases );
         free( records );
         return;
      }
      else {
         fprintf( fptr, "%5d   %5d\n", vertices, edges );
         ptr = records;
         for ( i = 0; i < edges; i++ ) {
            fprintf( fptr, "%s\n", ptr -> record );
            ptr++;
         }
      }
   }
   fclose( fptr );
   
   free( adj_matrix );
   free( vertex_degree );
   free( aux );
   free( aliases );
   free( records );
}

/* This function serves the same purpose as strcmp.  However,
   the arguments in compare_strings are of type const void*, thus
   providing compatibility with the ANSI function prototype header for qsort.
   If strcmp is substituted for compare_strings in the qsort
   invocation, in ANSI C a warning will be issued.
   (Notice that in ANSI C, the const void* types are automatically converted to
   const char* when strcmp is called.)
 */
int compare_strings( const void* a, const void* b )
{
   return strcmp( a, b );
}

void rename_vertices( int* aliases, int count )
{
   init_array( aliases, count );
   permute( aliases, count );
}

void rename_edges( Record* ptr, int* aliases, int count )
{
   int   i, vertex1, vertex2, v1, v2;
   char  buffer[ RecordSize ];

   for ( i = 0; i < count; i++ ) {
      sscanf( ptr -> record, "%5d  %5d",
              &vertex1, &vertex2 );
      if ( ( v1 = aliases[ vertex1 - 1 ] + 1 ) <
           ( v2 = aliases[ vertex2 - 1 ] + 1 ) )
         sprintf( buffer, "%5d   %5d   %5d", v1, v2, 1 );
      else
         sprintf( buffer, "%5d   %5d   %5d", v2, v1, 1 );
      strcpy( ptr -> record, buffer );
      ptr++;
   }
}

/* This function generates a random connected simple graph with
   v vertices and max(v-1,e) edges.  The graph can be weighted
   (weight_flag == 1) or unweighted (weight_flag != 1). If
   it is weighted, the weights are in the range 1 to max_wgt.
   It is assumed that e <= v(v-1)/2. (In this program, this assured
   because of the call to fix_imbalanced_graph.)

   To  generate  a  random  connected  graph,  we begin by
   generating  a  random  spanning  tree.  To generate a random
   spanning  tree,  we  first  generate  a  random  permutation
   tree[0],...,tree[v-1]. (v = number of vertices.)
   We  then  iteratively  add edges  to form a
   tree.  We  begin with the tree consisting of vertex tree[0] and
   no   edges.   At   the   iterative   step,  we  assume  that
   tree[0],tree[1],...,tree[i-1]  are  in  the  tree.  We  then add vertex
   tree[i]  to     the    tree    by    adding    the    edge
   (tree[i],tree[rand(i)]).  (This  construction  is similar to
   that  of  Prim's algorithm.) Finally, we add random edges to
   produce the desired number of edges.
 */
void random_connected_graph( int v,
                             int e,
                             int max_wgt,
                             int weight_flag,
                             char* out_file )
{
   int i, j, count, index, *adj_matrix, *tree;

   if ( ( adj_matrix = ( int * ) calloc( v * v, sizeof( int ) ) )
        == NULL ) {
      printf( "Not enough room for this size graph\n" );
      return;
   }

   if ( ( tree = ( int * ) calloc( v, sizeof( int ) ) ) == NULL ) {
      printf( "Not enough room for this size graph\n" );
      free( adj_matrix );
      return;
   }

   printf( "\n\tBeginning construction of graph.\n" );

   /*  Generate a random permutation in the array tree. */
   init_array( tree, v );
   permute( tree, v );

   /*  Next generate a random spanning tree.
       The algorithm is:

         Assume that vertices tree[ 0 ],...,tree[ i - 1 ] are in
         the tree.  Add an edge incident on tree[ i ]
         and a random vertex in the set {tree[ 0 ],...,tree[ i - 1 ]}.
    */

   for ( i = 1; i < v; i++ ) {
      j = ran( i );
      adj_matrix[ tree[ i ] * v + tree[ j ] ] =
         adj_matrix[ tree[ j ] * v + tree[ i ] ] =
         weight_flag ? 1 + ran( max_wgt ) : 1;
   }

   /* Add additional random edges until achieving at least desired number */

   for ( count = v - 1; count < e; ) {
      i = ran( v );
      j = ran( v );

      if ( i == j )
         continue;

      if ( i > j )
         swap( &i, &j );

      index = i * v + j;
      if ( !adj_matrix[ index ] ) {
         adj_matrix[ index ] = weight_flag ? 1 + ran( max_wgt ) : 1;
         count++;
      }
   }

   print_graph( v, count, out_file, adj_matrix, Undirected );

   free( tree );
   free( adj_matrix );
}

/* This function generates a random undirected graph with
   v vertices and e edges.  The graph can be weighted
   (weight_flag == 1) or unweighted (weight_flag != 1). If
   it is weighted, the weights are in the range 1 to max_wgt.
   The graph can be simple (simple_flag == 1) or not simple
   (simple_flag != 1). If the graph is simple, it is assumed
   that e <= v(v-1)/2. (In this program, this assured because
   of the call to fix_imbalanced_graph.)
 */
void random_graph( int v,
                   int e,
                   int max_wgt,
                   int weight_flag,
                   int simple_flag,
                   char* out_file )
{
   int i, j, count, index, *adj_matrix;

   /* generate a not necessarily simple random graph */
   if ( !simple_flag ) {
      r_graph( v, e, max_wgt, weight_flag, out_file );
      return;
   }

   /* generate a simple random graph */
   if ( ( adj_matrix = ( int * ) calloc( v * v, sizeof( int ) ) )
        == NULL ) {
      printf( "Not enough room for this size graph\n" );
      return;
   }

   for ( count = 0; count < e; ) {
      i = ran( v );
      j = ran( v );

      if ( i == j )
         continue;

      if ( i > j )
         swap( &i, &j );

      index = i * v + j;
      if ( !adj_matrix[ index ] ) {
         adj_matrix[ index ] = weight_flag ? 1 + ran( max_wgt ) : 1;
         count++;
      }
   }

   print_graph( v, e, out_file, adj_matrix, Undirected );

   free( adj_matrix );
}

/* This function generates a random, not necessarily simple
   graph. The graph can be weighted (weight_flag == 1)
   or unweighted (weight_flag != 1). It can be interpreted
   as directed or undirected.

   Random edges are generated and written directly to the
   output file.
 */
void r_graph( int v,
              int e,
              int max_wgt, int weight_flag,
              char* out_file )
{
   FILE *fp;
   int i;

   if ( ( fp = fopen( out_file, "w" ) ) == NULL ) {
      printf( "Unable to open file %s for writing\n", out_file );
      return;
   }

   printf( "\n\tWriting graph to file %s.\n", out_file );

   fprintf( fp, "%5d   %5d\n", v, e );

   for ( i = 0; i < e; i++ )
      fprintf( fp, "%5d   %5d   %5d\n", 1 + ran( v ), 1 + ran( v ),
               weight_flag ? 1 + ran( max_wgt ) : 1 );

   printf( "\tGraph is written to file %s.\n", out_file );
   fclose( fp );
}

/* This function generates a random directed graph with
   v vertices and e edges.  The graph can be weighted
   (weight_flag == 1) or unweighted (weight_flag != 1). If
   it is weighted, the weights are in the range 1 to max_wgt.
   The graph can be simple (simple_flag == 1) or not simple
   (simple_flag != 1). If the graph is simple, it is assumed
   that e <= v(v-1). (In this program, this assured because
   of the call to fix_imbalanced_graph.)
 */
void random_digraph( int v,
                     int e,
                     int max_wgt,
                     int weight_flag,
                     int simple_flag,
                     char* out_file )
{
   int i, j, count, index, *adj_matrix;

   /* generate a not necessarily simple random digraph */
   if ( !simple_flag ) {
      r_graph( v, e, max_wgt, weight_flag, out_file );
      return;
   }

   /* generate a simple random digraph */
   if ( ( adj_matrix = ( int * ) calloc( v * v, sizeof( int ) ) )
        == NULL ) {
      printf( "Not enough room for this size graph\n" );
      return;
   }

   for ( count = 0; count < e; ) {
      i = ran( v );
      j = ran( v );

      if ( i == j )
         continue;

      index = i * v + j;
      if ( !adj_matrix[ index ] ) {
         adj_matrix[ index ] = weight_flag ? 1 + ran( max_wgt ) : 1;
         count++;
      }
   }

   print_graph( v, e, out_file, adj_matrix, Directed );
   free( adj_matrix );
}

/* This function writes a directed acyclic graph to one file and a
   topological sort to a second file. The graph can be weighted or
   unweighted. It is assumed that e <= v(v-1)/2. (In this program,
   this assured because of the call to fix_imbalanced_graph.)

   To   generate   a  directed  acyclic  graph,  we  first
   generate  a  random  permutation  dag[0],...,dag[v-1].
   (v = number of vertices.)
   This random permutation serves as a topological
   sort  of  the  graph. We then generate random edges of the
   form (dag[i],dag[j]) with i < j.
 */
void directed_acyclic_graph( int v,
                             int e,
                             int max_wgt,
                             int weight_flag,
                             char* out_file,
                             char* dag_file )
{
   int i, j, count, index, *adj_matrix, *dag;
   FILE *fptr;

   if ( ( adj_matrix = ( int * ) calloc( v * v, sizeof( int ) ) )
        == NULL ) {
      printf( "Not enough room for this size graph\n" );
      return;
   }

   if ( ( dag = ( int * ) calloc( v, sizeof( int ) ) ) == NULL ) {
      printf( "Not enough room for this size graph\n" );
      free( adj_matrix );
      return;
   }

   printf( "\n\tBeginning construction of graph.\n" );

   /*  Generate a random permutation in the array dag. */
   init_array( dag, v );
   permute( dag, v );

   for ( count = 0; count < e; ) {
      if ( ( i = ran( v ) ) == ( j = ran( v ) ) )
         continue;
      if ( i > j )
         swap( &i, &j );
      i = dag[ i ];
      j = dag[ j ];
      index = i * v + j;
      if ( !adj_matrix[ index ] ) {
         adj_matrix[ index ] = weight_flag ? 1 + ran( max_wgt ) : 1;
         count++;
      }
   }

   print_graph( v, e, out_file, adj_matrix, Directed );

   if ( ( fptr = fopen( dag_file, "w" ) ) == NULL )
      printf( "\n\t\t Could not open file %s.\n", dag_file );
   else {
      for ( i = 0; i < v; i++ )
         fprintf( fptr, "%5d\n", 1 + dag[ i ] );
      fclose( fptr );
   }

   free( adj_matrix );
   free( dag );
}

/* This function generates a weighted complete graph.  The
   graph can be directed or not.
 */
void complete_graph( int v,
                     int max_wgt,
                     int dir_flag,
                     char* out_file )
{
   int i, j;
   FILE *fp;

   if ( ( fp = fopen( out_file, "w" ) ) == NULL ) {
      printf( "Unable to open file %s for writing.\n", out_file );
      return;
   }
   printf( "\n\tWriting graph to file %s.\n", out_file );

   fprintf( fp, "%5d   %5d\n", v,
            dir_flag ? v * ( v - 1 ) : v * ( v - 1 ) / 2 );

   for ( i = 1; i < v; i++ )
      for ( j = i + 1; j <= v; j++ ) {
         fprintf( fp, "%5d   %5d   %5d\n", i, j, 1 + ran( max_wgt ) );
         if ( dir_flag )
            fprintf( fp, "%5d   %5d   %5d\n", j, i, 1 + ran( max_wgt ) );
      }

   fclose( fp );

   printf( "\tGraph is written to file %s.\n", out_file );
}

/* This function writes a simple graph with a Hamiltonian cycle to one
   file and a Hamiltonian cycle to another file. The graph will
   have max(e,v) edges. The graph can be directed or undirected. It is
   assumed that e <= v(v-1)/2 if the graph is undirected, and that
   e <= v(v-1) if the grpah is directed. (In this program,
   this assured because of the call to fix_imbalanced_graph.)

   To generate a random graph with a
   Hamiltonian cycle, we begin with a random permutation
   ham[0],...,ham[v-1] (v = number of vertices).  We then generate edges

   (ham[0],ham[1]),(ham[1],ham[2]),...,(ham[v-2],ham[v-1]),(ham[v-1],ham[0]),

   so that the graph has the Hamiltonian cycle

                     ham[0],...,ham[v-1],ham[0].

   Finally,  we  add random edges to produce the desired number
   of edges.
 */
void hamiltonian_cycle_graph( int v,
                              int e,
                              int dir_flag,
                              char* out_file,
                              char* ham_file )
{
   int i, j, k, l, count, index, *adj_matrix, *ham;
   FILE *fptr;

   if ( ( adj_matrix = ( int * ) calloc( v * v, sizeof( int ) ) )
        == NULL ) {
      printf( "Not enough room for this size graph\n" );
      return;
   }

   if ( ( ham = ( int * ) calloc( v, sizeof( int ) ) ) == NULL ) {
      printf( "Not enough room for this size graph\n" );
      free( adj_matrix );
      return;
   }

   printf( "\n\tBeginning construction of graph.\n" );

   /*  Generate a random permutation in the array ham. */
   init_array( ham, v );
   permute( ham, v );

   if ( ( fptr = fopen( ham_file, "w" ) ) == NULL ) {
      printf( "\n\t\t Could not open file %s.\n", ham_file );
      free( adj_matrix );
      free( ham );
      return;
   }

   /* print Hamiltonian cycle and store required edges */
   for ( i = 0; i < v; i++ ) {
      fprintf( fptr, "%5d\n", 1 + ham[ i ] );
      k = ham[ i ];
      l = ham[ ( i + 1 ) % v ];
      if ( k > l && !dir_flag )
         swap( &k, &l );
      adj_matrix[ k * v + l ] = 1;
   }

   fprintf( fptr, "%5d\n", 1 + ham[ 0 ] );
   fclose( fptr );
   free( ham );

   for ( count = v; count < e; ) {
      if ( ( i = ran( v ) ) == ( j = ran( v ) ) )
         continue;
      if ( i > j && !dir_flag )
         swap( &i, &j );
      index = i * v + j;
      if ( !adj_matrix[ index ] ) {
         adj_matrix[ index ] = 1;
         count++;
      }
   }

   print_graph( v, count, out_file, adj_matrix, dir_flag );
}

/* This function generates a random network.  The user selects
   the number of vertices (v) and edges (e), the maximum weight per edge, and the
   the output file name. It is assumed that e <= v**2 - 3v + 3. (In this
   program, this assured because of the call to fix_imbalanced_graph.) The
   function will produce a network with at least e edges. In some rare,
   cases, the actual number of edges may be slightly more than e because
   of the following:

      Each vertex is guaranteed to be on a directed path from the source
      to the sink. This is achieved by generating random directed paths
      from the source to the sink until every vertex is on such a path.
 */
void network_flow_graph( int vertices,
                         int edges,
                         int max_weight,
                         char* outfile )
{
   int  start = 0,
        finish = vertices - 1,
        edge_count = 0,
        from, to,
        *adj_matrix,
        *closed_list;

   if ( vertices == 1 ) {
      warning( SingleVertexNetwork );
      return;
   }

   /* adjacency matrix to represent graph */
   if ( ( adj_matrix = ( int * ) calloc( vertices * vertices, sizeof( int ) ) ) == NULL ) {
      warning( InsufficientStorage );
      return;
   }

   /* closed_list[ i ] == True if vertices[ i ] is on some path from start to finish */
   if ( ( closed_list = ( int * ) calloc( vertices, sizeof( int ) ) ) == NULL ) {
      warning( InsufficientStorage );
      free( adj_matrix );
      return;
   }

   /* put each vertex on some path from start to finish */
   closed_list[ start ] = True;
   do {
      /* from start to some vertex not already picked */
      to = pick_to( closed_list, vertices );
      adj_matrix[ to ] = ran( max_weight ) + 1;
      edge_count++;

      /* from current vertex to vertex != finish */
      if ( to != finish ) {
         from = to;
         closed_list[ finish ] = False;
         while ( ( to = pick_to( closed_list, vertices ) ) != finish ) {
            adj_matrix[ from * vertices + to ] = ran( max_weight ) + 1;
            edge_count++;
            from = to;
         }

         /* from vertex to finish */
         adj_matrix[ from * vertices + finish ] = ran( max_weight ) + 1;
         edge_count++;
      }
   } while ( unpicked_vertices_p( closed_list, vertices ) );

   /* add additional edges as needed */
   while ( edge_count < edges ) {
      from = ran( vertices - 1 );
      to = pick_to( NULL, vertices );
      if ( !adj_matrix[ from * vertices + to ]  &&  from != to ) {
         adj_matrix[ from * vertices + to ] = ran( max_weight ) + 1;
         edge_count++;
      }
   }
   print_graph( vertices, edge_count, outfile, adj_matrix, Directed );

   /* free storage */
   free( adj_matrix );
   free( closed_list );
}

int pick_to( int* closed, int vertices )
{
   int  vertex;

   if ( closed == NULL )        /* pick any vertex, even if picked before... */
      do {
         vertex = ran( vertices );
      } while ( vertex == 0 );  /* ...except for start */
   else {
      do {
         vertex = ran( vertices );
      } while ( closed[ vertex ] );
      closed[ vertex ] = True;
   }
   return ( vertex );
}

int unpicked_vertices_p( int* closed, int n )
{
   int  i = 1;

   while ( i < n )
      if ( !closed[ i ] )
         return ( True );
      else
         i++;
   return ( False );
}

/*** ran, etc. ***/
void seed_ran( void )
{
   srand( ( unsigned short ) time( NULL ) );
}

/* Return a random integer between 0 and k-1 inclusive. */
int ran( int k )
{
   return rand() % k;
}

void print_graph( int v,
                  int e,
                  char* out_file,
                  int* adj_matrix,
                  int dir_flag )
{
   int i, j, index;
   FILE *fp;

   if ( ( fp = fopen( out_file, "w" ) ) == NULL ) {
      printf( "Unable to open file %s for writing.\n", out_file );
      return;
   }
   printf( "\n\tWriting graph to file %s.\n", out_file );

   fprintf( fp, "%5d   %5d\n", v, e );

   if ( !dir_flag )
      for ( i = 1; i < v; i++ )
         for ( j = i + 1; j <= v; j++ ) {
            index = ( i - 1 ) * v + j - 1;
            if ( adj_matrix[ index ] )
               fprintf( fp, "%5d   %5d   %5d\n", i, j, adj_matrix[ index ] );
         }
   else
      for ( i = 1; i <= v; i++ )
         for ( j = 1; j <= v; j++ ) {
            index = ( i - 1 ) * v + j - 1;
            if ( adj_matrix[ index ] )
               fprintf( fp, "%5d   %5d   %5d\n", i, j, adj_matrix[ index ] );
         }
   fclose( fp );
   printf( "\tGraph is written to file %s.\n", out_file );
}

/* randomly permute a[ 0 ],...,a[ n - 1 ] */
void permute( int* a, int n )
{
   int i;

   for ( i = 0; i < n - 1; i++ )
      swap( a + i + ran( n - i ), a + i );
}

void swap( int* a, int *b )
{
   int temp;

   temp = *a;
   *a = *b;
   *b = temp;
}

/* set a[ i ] = i, for i = 0,...,end - 1 */
void init_array( int* a, int end )
{
   int i;

   for ( i = 0; i < end; i++ )
      *a++ = i;
}

/* Get integer input from user. If not an int, prompt for correct
   value. indent gives the proper indentation for error message.
 */
int get_int( char* indent )
{
   char buff[ 30 ];
   int val;

   for ( ; ; ) {
      scanf( "%s", buff );
      if ( sscanf( buff, "%d", &val ) != 1 )
         printf( "%s%s", indent, "***** Illegal input.  Reenter value: " );
      else
         return val;
   }
}