# CSM-Experiments
Source code and data for paper "An Experimental Study of Indexes in Continuous Subgraph Matching".

# Envoriement
To compile the code, users need c++ version 11 or higher.  

# Data and Data Organization

The datasets and queries that we use in our experiments can be found at .

The file of initial graph starts with the total vertex count, as we need this count to allocate vertex tables in the graph storage.

Then each line describes an edge, with the form

source vertex, label of source vertex, destination vertex, label of destination vertex, label of the edge.

The file of update stream has the same format except that it does not have the vertex count.

The query graphs are organized in the following format:

The first line is the number of vertcies, followed by the neighborhood information of each vertex.

For each vertex u, it starts with 

vertex ID, vertex label

Then followed by the number of lable groups. For each group, it starts with a neighbor vertex label and an edge label, u is connected to neighbors with this vertex label by edges with the given  



We also provide the c++ code that we use to extract the query graphs.


# Code Compile
To select between differnt algorihtms, users should edit the "setting" file and set the following variables:

### index_type 
control the indexing method of the algorithm, can be set to the following value:

0:   NI, no index, 

1:   NLF, NLF filter without materialization, proposed by NewSp,

2:   NLFI, materialized NLF index, proposed by RapidFlow

3:   DCG, index proposed by Turbofux 

4:   DCS, index proposed by Symbi

5:   CaliG, index proposed by CaliG

### enumeration type
control the enumeration method of the algorithm, can be set to the following value:

0:   GraphFlow

1:   NewSP

2:   RapidFlow

3:   TurboFlux

4:   Symbi

5:   CaliG

When setting the index type and the enumeration type to the same number, the combination will be the original form of a CSM algorithm, like 0 for GraphFlow, 1 for NewSP. But we also encourage to try different combinations.

### index_test
When setting this value to true, we will only update the index but do not enumerate matches. This option is used to test index update cost.

### index_NLF
When setting this value to false, DCG and DCS will not use NLF filter.

### use_edge_view
When setting this value to false, indexes will not build edge indexes, but only vertex indexes.

After setting the above values, users can compile the code with command

#### g++ -O3 -std=c++11 -o process processing.cpp

An executable file "process" will be generated.


# Execution

Users can 



# Result and Result Porcessing.

After execution, there will be two files generated in the output path. One is record.txt and the other is result.txt.

In result.txt, each query start with "query : [query number]". In the following,  we list the index of the edges upon whose insertion there are incremental matches generated, and the number of incremental matches.

In record.txt, we record the excecution result of each query in the following form:




。
