/*@CR Too generic - both file name and ifndef. Common convention is _MODULES_H_ (based on file name)*/
#ifndef HEADER_FILE
#define HEADER_FILE

static const double epsilon = 0.00001;

struct graph
{
  struct _spmat* A;
  long* vectorDegrees;
  long M;
};

/* @CR the comment is a bit non-informative. Either you recognize the need for previous knowledge
   from the assignment (such as A, and M) so its not necessary, or be more verbose and maybe even
  copy the definitions from the assignment as comments.*/
/*Represents a group of divisions*/
struct division
{
	int len;  /*number of divisions // @CR why not just call it numDivisions and remove the comment?*/
	struct node* divisions; /*"P" in algorithm // @CR why not just call it P (if you called the graph A and M).*/
};

/*@CR not valid - must be declared before division which uses it (or use forward declaration).
*represents a single group in divisions*
 * TODO: change short to bit? or perhaps boolean? need to see who's smaller
 *@CR about to do - char is the shortest, if you know the max length is 32 you can use int and use bits operations.
*@CR need to think if this is the correct data structure to use, we talked about an alternative, but could be other
    *options and I'm not fully into the exercise to give meaningfull feedback here.*/
struct divisionGroup
{
	int groupSize;
	struct _spmat* groupSubmatrix;
	int* sumOfRows; /*(computeS populates it).  0 by default*/
	int* groupMembers;

};

struct shiftedDivisionGroup
{
	struct divisionGroup* group;
	double norm;
};


struct node
{
	union data {
		int num;
		struct divisionGroup* group;
	} data;
	struct node* next;
};



#endif
