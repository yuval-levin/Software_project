// @CR Too generic - both file name and ifndef. Common convention is _MODULES_H_ (based on file name)
#ifndef HEADER_FILE
#define HEADER_FILE

struct graph
{
  struct spmat_ll A;
  long* vectorDegrees;
  long M;
};

// @CR the comment is a bit non-informative. Either you recognize the need for previous knowledge
//   from the assignment (such as A, and M) so its not necessary, or be more verbose and maybe even
//   copy the definitions from the assignment as comments.
/*Represents a group of divisions*/
struct division
{
	int len;  //number of divisions // @CR why not just call it numDivisions and remove the comment?
	struct node* divisions; //"P" in algorithm // @CR why not just call it P (if you called the graph A and M).
};

// @CR not valid - must be declared before division which uses it (or use forward declaration).
/*represents a single group in divisions*
 * TODO: change short to bit? or perhaps boolean? need to see who's smaller*/
// @CR about to do - char is the shortest, if you know the max length is 32 you can use int and use bits operations.
// @CR need to think if this is the correct data structure to use, we talked about an alternative, but could be other
//    options and I'm not fully into the exercise to give meaningfull feedback here.
struct divisionGroup
{
	int groupSize;
	struct spmat_ll* groupSubmatrix;
	int* sumOfRows; //(computeS populates it).  0 by default
	/* For n number of nodes in graph , group vector is size n
	 * groupVector[i]=0 if node i not in group. 1 if in group.
	 * #{i | groupVector[i]=1 } = groupSize
	*/
};

struct node
{
	union data data;
	struct node *next;
};

union data {
	int num;
	struct divisionGroup* group;
};

#endif
