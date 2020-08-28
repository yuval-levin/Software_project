#include<stdio.h>
#include <stdlib.h>
#include "modules.h"

int main(int args, char** argv)
{


}

/* PARAMS : n - number of nodes
 * The division is empty by default*/

//TODO SHOULD WE USE TYPEDEF FOR DIVISION?
//TODO: add function called "free division" to make sure there are no memoryleaks
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

/*PARAMS: P,O are the divisions from the algorithm.
		They are created and initialized outside of the function,
		to enable memory freeing easily and outside of the function */
struct division* Algorithm3(struct division* P, struct division* O)
{

	struct divisionGroup* g;
	short  gWasDivided; //1 for "true", 0 for "false";
	while (P->len > 0)
	{
		g = (P->divisions)[0];
		gWasDivided = Algorithm2(&P,0);
		if (gWasDivided == 0) updateOAndP(&P,&O);
		else addToOSize1Groups(&O,&P);
		// STEP 3.4.2 is automatic, since Algorithm2 edits P.
		// So g1,g2 still remain in P is their size is > 0.
	}

	return &O;
}

void addToOSize1Groups(struct division* P,struct divison* O)
{
	int indexG1 = (P->len)-1;
	int indexG2 = (P->len)-2;
	if (P->divisions[indexG1]->groupSize == 1) updateOAndP(&P,&O,indexG1);
	if (P->divisions[indexG2]->groupSize == 1) updateOAndP(&P,&O,indexG2);
}

/*
 * PARAMS : P - division to remove group from
 * PARAMS : O - division to add group to
 * PARAMS : groupIndex - groupIndex in P
 *
 * DESC: Removes from P group in groupIndex.
 * 		 Adds that group to O
 */
void updateOAndP(struct division* P,struct divison* O,int groupIndex)
{
	addGroupToDivision(&O,&(P->divisions[groupIndex]));
	removeLastGroupFromDivision(&P);
}

void removeLastGroupFromDivision(struct division* D)
{
	D->len = (D->len)-1;
	D->divisions[D->len] = NULL; //todo: perhaps delete this row, since len-1 is enough?
	//next time someone will add or iterate it will be up to new len..
	//make sure this row doesn't delete the last group division, only the reference to it in D
}
/*
 * PARAMS : D - division Pointer
 * PARAMS : group - group to be added to D
 *
 * DESC : Ddd group to D at the end
 */
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

