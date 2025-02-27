#pragma once
#include<time.h>
#include "../indexes/index-list.h"
#include "enumeration/enumeration-list.h"

#define index_type 0
// control the indexing method of the algorithm, can be set to the following value:
//0:   NI, no index, 
//1:   NLF, NLF filter without materialization, proposed by NewSp,
//2:   NLFI, materialized NLF index, proposed by RapidFlow
//3:   DCG, index proposed by Turbofux 
//4:   DCS, index proposed by Symbi
//5:   CaliG, index proposed by CaliG

#define enumeration type 0
// control the enumeration method of the algorithm, can be set to the following value:
//0:   GraphFlow
//1:   NewSP
//2:   RapidFlow
//3:   TurboFlux
//4:   Symbi
//5:   CaliG

// when setting the index type and the enumeration type to the same number, the combination will be the original form of a CSM algorithm, like 0 for GraphFlow, 1 for NewSP. But we also encourage to try different combinations.


#define index_test false
// when set this value to true, we will only update the index but do not enumerate matches. This option is used to test index update cost.
#define index_NLF true
// when set this value to false, DCG and DCS will not use NLF filter.
#define use_edge_view true
// when set this value to false, indexes will not build edge indexes, but only vertex indexes.


 

typedef symbi_index my_index;
//typedef TF_index my_index;
//typedef NLF_index my_index;
//typedef CaliG_index my_index; 

// change this type to change among different indexes. User can choose between the following indexes;
// NLF_index : the index of RapidFlow. namely NLFI in the paper
//TF_index: index of Turboflux, namely DCG;
//symbi_index: the index of symbi, namely DCS
//CaliG_index ; index of CaliG






typedef NewSP solution; // change this type to change among different enumberation method. User can choose between the following indexes;
// GraphFlow
// NewSP
//RapidFlow
//Symbi
//TurboFlux
//CaliG

#define use_NLF false // if set this flag to false, No-index enumberation algorithms will not perform NLF check


clock_t start_clock;

unsigned int time_limit = 3600;

unsigned long long inter_result;
unsigned int filtered_edge;

unsigned int label_filtered;


unsigned int expanding_num;
unsigned int neighbor_sum;
unsigned int filtered_expanding;

// Note that indexes other than NLF_index also uses the NLF filtering rule. This is requested by a recent work, we also think it is better in this way 
// to evaluate the influence of indexes in this way.


