#include<stdio.h>
#include <stdlib.h>
#include "modules.h"

// @CR as discussed, should probably come in its own file.
int main(int args, char** argv)
{

}

/* PARAMS : n - number of nodes
 * The division is empty by default*/

//TODO SHOULD WE USE TYPEDEF FOR DIVISION?
//TODO: add function called "free division" to make sure there are no memoryleaks
// @CR define private functions which are not exported as static (it means that the linking scope is module only so it is not exported).
struct division createDivision(int n)
{
	struct division div;

	div=(struct division*)malloc(sizeof(struct division));
	div->len = 0;
	div->divisions = (struct divisionGroup*)malloc(n*sizeof(struct divisionGroup));
	//div[0] = createDefaultDivision(n);

	return div;
}

/*PARAMS : n is number of nodes in graph
 * DESC : Returns a divisonGroup with all nodes in graph
 */
struct divisionGroup createDefaultDivision(int n)
{
	int i;
	struct divisonGroup group = (struct divisionGroup)malloc(sizeof(struct divisionGroup));

	group->groupSize = n;
	group->groupVector = (short*)malloc(n*sizeof(short));
	for (i = 0; i < n; i++)
	{
		group->groupVector[i] = 1;
	}
	return group;
}

void addGroupToDivision(struct division* D,struct divisionGroup* group)
{
	D->divisions[D->len] = group; //TODO: make sure this is okay
	D->len = (D->len)+1;
}



/*
 *PARAMS: P - division POINTER (or else struct division would be copied .. )
 *PARAMS: g - the index of the group from P we're dividing
 *
 *DESCRIPTION : Implementation of Algorithm 2. We edit the division P.
 *DESCRIPTION : The new groups g1,g2 are added to the list of groups AT THE END
 *DESCRIPTION : If g is unidivisble we change nothing (g is STILL FIRST)
 *DESCRIPTION : if g is undivisible we return 0
 *DESCRIPTION : if g is divisble we return 1
 * TODO: if we always call Algorithm 2 with group 0 , consider erasing g param.
 */
short Algorithm2(struct division * P,int g)
{

}

