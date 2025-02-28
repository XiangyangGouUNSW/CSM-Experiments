#pragma once
#include<time.h>


#define index_type 0
// control the indexing method of the algorithm, can be set to the following value:
//0:   NI, no index, 
//1:   NLF, NLF filter without materialization, proposed by NewSp,
//2:   NLFI, materialized NLF index, proposed by RapidFlow
//3:   DCG, index proposed by Turbofux 
//4:   DCS, index proposed by Symbi
//5:   CaliG, index proposed by CaliG

#define enumeration_type 1
// control the enumeration method of the algorithm, can be set to the following value:
//0:   GraphFlow
//1:   NewSP
//2:   RapidFlow
//3:   TurboFlux
//4:   Symbi
//5:   CaliG

// when setting the index type and the enumeration type to the same number, the combination will be the original form of a CSM algorithm, like 0 for GraphFlow, 1 for NewSP. But we also encourage to try different combinations.


#define indexing_test false
// when set this value to true, we will only update the index but do not enumerate matches. This option is used to test index update cost.
#define index_NLF true
// when set this value to false, DCG and DCS will not use NLF filter.
#define enable_edge_view true
// when set this value to false, indexes will not build edge indexes, but only vertex indexes.






clock_t start_clock;

unsigned int time_limit = 3600;

unsigned long long inter_result;
unsigned int filtered_edge;

unsigned int label_filtered;


unsigned int expanding_num;
unsigned int neighbor_sum;
unsigned int filtered_expanding;


#include "indexes/index-list.h"
 
#if index_type==0
	 typedef NLF_index my_index;
	 #define use_NLF false
	 #define enable_global_index false
#elif index_type==1
	typedef NLF_index my_index;
	 #define use_NLF true
	 #define enable_global_index false
#elif index_type==2
	typedef NLF_index my_index;
	 #define use_NLF true	
	 #define enable_global_index true
#elif index_type==3
	typedef TF_index my_index;
	 #define use_NLF true	
	 #define enable_global_index true
#elif index_type==4
	typedef symbi_index my_index;
	 #define use_NLF true	
	 #define enable_global_index true
#elif index_type==5
	typedef CaliG_index my_index;
	 #define use_NLF true	
	 #define enable_global_index true
#else
	typedef NLF_index my_index;
	 #define use_NLF true	
	 #define enable_global_index true
#endif



#include "enumeration/enumeration-list.h"

#if enumeration_type==0
	 typedef GraphFlow solution;
#elif enumeration_type==1
	 typedef NewSP solution;
#elif enumeration_type==2
	 typedef RapidFlow solution;
#elif enumeration_type==3
	 typedef TurboFlux solution;
#elif enumeration_type==4
	 typedef Symbi solution;
#elif enumeration_type==5
	 typedef CaliG solution;
#else
	 typedef NewSP solution;
#endif




