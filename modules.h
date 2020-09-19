
#ifndef _MODULES_H
#define _MODULES_H

static const double epsilon = 0.00001;

struct graph
{
  struct _spmat* A;
  double* vectorDegrees;
  double M;
  double* degreesDividedByM; /* Ki/M is used many times, so we'll store*/
  int numOfNodes;
};

/*Represents a group of divisions*/
struct division
{
	int len;  /*number of divisions // @CR why not just call it numDivisions and remove the comment?*/
	struct node* divisions; /*"P" in algorithm // @CR why not just call it P (if you called the graph A and M).*/
};


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
	struct node* next;
	union  _data {
		int num;
		struct divisionGroup* group;
	}data;

};


#endif
