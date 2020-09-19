
#ifndef _MODULES_H
#define _MODULES_H

static const double epsilon = 0.00001;

struct graph
{
  struct _spmat* A;
  double* vectorDegrees;
  double M;
  double* degreesDividedByM; /* Ki/M is used many times, so it is stored*/
  int numOfNodes;
};

/*Represents a group of divisions*/
struct division
{
	int len;  /*number of divisions*/
	struct node* divisions; /*"P" in algorithm */
};

/*represents a division*/
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
