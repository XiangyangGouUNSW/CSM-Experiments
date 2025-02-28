# CSM-Experiments
Source code and data for paper "An Experimental Study of Indexes in Continuous Subgraph Matching".

# Envoriement
To compile the code, users need c++ version 11 or higher.  

# Data and Data Organization

The datasets and queries that we use in our experiments can be found at https://unsw-my.sharepoint.com/:u:/g/personal/z3544923_ad_unsw_edu_au/EQVwnHa9OONGgy6UqW1we0QB73QF0wi_nIqzdKKcA_fujQ?e=aw9uSb .

The file of initial graph starts with the total vertex count, as we need this count to allocate vertex tables in the graph storage.

Then each line describes an edge, with the form

[source vertex] [label of source vertex] [destination vertex] [label of destination vertex] [label of the edge]

The file of update stream has the same format except that it does not have the vertex count.

The query graphs are organized in the following format:

The first line is the number of vertcies, followed by the neighborhood information of each vertex.

For each vertex u, it starts with 

[vertex ID] [vertex label]

Then followed by the number of lable groups. For each group, it starts with 

[neighbor vertex label]  [edge label]  [number of neighbors in this group]

u is connected to neighbors with this vertex label by edges with the given edge label. Then we will list the IDs of the neighbors in this group.

For example, a neighborhood information of vertex 1 maybe as follows:

1 0

1 0 2

2 3

It means that vertex 1 has label 0, and is connected to 2 label 1 neighbors via edges with label 0. These neighbors are vertex 2 and 3.


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

Users can execute the experiments with conmand

#### ./process [path of initial graph file] [path of update stream file] [path of query set file] [output path]

The output path should be a folder and the other paths are paths to files. For example


#### ./process LiveJournal/data/initial_data.txt LiveJournal/data/update_data.txt LiveJournal/query/size6/dense-query.txt Output/


# Result and Result Porcessing.

After execution, there will be two files generated in the output path. One is 

#### record.txt 

and the other is

#### result.txt.

In the result file, each query start with "query : [query number]". In the following,  we list the indexes of the edges upon whose insertion there are incremental matches generated, and the number of incremental matches.

In the record file, we record the excecution result of each query in the following form:

#### query [query number]

#### [offline index building time] [Whether the query is solved] [total count of match results] [processing time]  [average intermediate result in the matching process induced by each update]

#### [total number vertex candidates in the index] [total number edge candidates in the index] [total memory usage of the index]

Note that for unsolved queries the total count of match results may be different for different algorithms. Because they have different processing speed and will find different number of matches in the time limit.


We also provide two c++ files for result processing.

The first is "result-process.cpp", used to coumpute the average values of all queries. Including the average processing time, number of solved and unsolved queries, average index memory usage and so on.

It can be complied with 

#### g++ -O3 -std=c++11 -o result-process result-processing.cpp

Then executed with 
#### result-process [path to the record file]

Note that when we compute the average match result count, we will only include the match result count of solved queries. This may result into different average match result count for differen algorithms with the same group of queries. Because theses algorithms may have different solved and unsolved query set. 

The second c++ file is "result-compare.cpp", which compares the processing time of two algorithms.

It can be complied with 

#### g++ -O3 -std=c++11 -o result-compare result-compare.cpp

Then executed with 
#### result-process [path to the record file of alg1]  [path to the record file of alg2]

This program will find the queries where at least one of the algrorithms have solved, and compute the improvement of alg2 against alg1 in this query. Them improvement is computed as 

####  1 - (prcossing time of alg2) / (prcossing time of alg1). 







。
