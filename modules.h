#ifndef HEADER_FILE
#define HEADER_FILE

struct graph
{
  struct spmat_ll A;
  long* vectorDegrees;
  long M;
};

/*Represents a group of divisions*/
struct division
{
	int len;  //number of divisions
	struct divisionGroup* divisions; //"P" in algorithm
};

/*represents a single group in divisions*
 * TODO: change short to bit? or perhaps boolean? need to see who's smaller*/
struct divisionGroup
{
	int groupSize;
	short* groupVector;
	/* For n number of nodes in graph , group vector is size n
	 * groupVector[i]=0 if node i not in group. 1 if in group.
	 * #{i | groupVector[i]=1 } = groupSize
	*/
};

#endif
