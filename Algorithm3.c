#include <stdio.h>
#include <stdlib.h>
#include "modules.h"
#include "Algorithm2.h"
#include "error_codes.h"
#include "modularity_maximization.h"
#include "spmat.h"

/* free a linked list of type spmat_node*/
static void free_list(struct spmat_node* cur) {
	struct spmat_node* next;
	while (cur != NULL) {
		next = cur->next;
		free(cur);
		cur = next;
	}
}

/* free the division groups*/
static void free_group(struct divisionGroup* group) {
	int i,length;
	struct spmat_node** rows;

	rows = group->groupSubmatrix->private;
	length = group->groupSubmatrix->n;

	for (i = 0; i < length; i++) {
		free_list(rows[i]);
	}

	free(group->sumOfRows);
	free(group->groupMembers);
	free_A(group->groupSubmatrix);
	free(group); /*TODO: delete?*/
}

/*free the  divisions given to us by Algorithm 3
 * Which are inside of O */
void free_division_group(struct division* O) {
	int i, length;
	struct node* currDiv;
	struct node* nextDiv;
	length = O->len;

	currDiv = O->divisions;
	for (i = 0; i < length; i++) {
		nextDiv = currDiv->next;
		free_group(currDiv->data.group);
		free(currDiv);
		currDiv = nextDiv;
	}
	free(O);
}

/*helper function for handling groups:
 * adds a groupdivision g in in start of division D*/
static void add_groupDivision(struct division* D, struct divisionGroup* g) {
	struct node* add = (struct node*) malloc(sizeof(struct node));
	if (add == NULL)
		panic(ERROR_MALLOC_FAILED);

	add->data.group = g;
	add->next = NULL;
	if (D->len == 0)
		D->divisions = add;
	else {
		add->next = D->divisions;
		D->divisions = add;
	}

	D->len = (D->len) + 1;
}

/* for handling groups:
 * removes  first divisiongroup from division D and returns it*/
static struct divisionGroup* remove_first_group(struct division* D) {
	struct node* group;
	struct node* nextGroup;
	struct divisionGroup* returnGroup;
	group = D->divisions;
	nextGroup = group->next;
	group->next = NULL;
	D->divisions = nextGroup;
	D->len = (D->len) - 1;
	returnGroup = group->data.group;
	free(group);
	return returnGroup;
}

/* for Algorithm 2.
 * After split,if group (here represented by g as input argument, but is  g1 or g2 in algo2) larger than 1 are put in P,
 * Any groups of size 1 are put in O */
static void update_division_post_split(struct divisionGroup* g, struct division* P,
		struct division* O) {
	if (g->groupSize == 1)
		add_groupDivision(O, g);
	else
		add_groupDivision(P, g);

}

/* splitByS helper, returns the size of g1*/
static int calc_size(double* vectorS, int length) {
	int i, size;
	size = 0;

	for (i = 0; i < length; i++) {
		if (vectorS[i] == 1) {
			size++;
		}
	}
	return size;
}

/* splitByS helper, updates the mat rows,
 * removes irrelevant nodes and updates sum of rows accordingly*/
static void update_mat_rows(double* vectorS, int numMembers, int groupIndicator,
		struct spmat_node** gtRows, int* gtSumOfRows) {
	int i_row;
	struct spmat_node *cur, *prev, *next;

	/* delete irrelevant nodes*/
	for (i_row = 0; i_row < numMembers; i_row++) {
		prev = NULL;
		cur = gtRows[i_row];
		while (cur != NULL) {
			next = cur->next;

			/* remove cur, prev remains the same*/
			if (vectorS[cur->index] != groupIndicator) {
				gtSumOfRows[i_row]--;
				if (prev == NULL) {
					gtRows[i_row] = cur->next;
				} else {
					prev->next = cur->next;
				}
				free(cur);
				/* don't remove cur, prev should advance*/
			} else {
				prev = cur;
			}
			cur = next;
		}
	}
}

/* splitByS helper, updates nodes index according to the new mat */
static void update_mat_rows_index(int* gGroupMembers, int nuMembers,
		struct spmat_node** gRows) {
	int i_row, count, i;
	int* mapNameToColIndex;
	struct spmat_node *cur;

	/* create a map, will be used to update index field to fit the new column index in the new spmat*/
	count = 0;
	mapNameToColIndex = (int*) malloc(
			(gGroupMembers[nuMembers - 1] + 1) * sizeof(int));
	if (mapNameToColIndex == NULL)
		panic(ERROR_MALLOC_FAILED);

	for (i = 0; i < nuMembers; i++) {
		mapNameToColIndex[gGroupMembers[i]] = count;
		count++;
	}

	/* delete irrelevant nodes*/
	for (i_row = 0; i_row < nuMembers; i_row++) {
		cur = gRows[i_row];
		while (cur != NULL) {
			cur->index = mapNameToColIndex[cur->node_name];
			cur = cur->next;
		}
	}
	free(mapNameToColIndex);
}

/* splitByS helper, updates group members of target, according to group indicator (+-1)*/
static void update_group_members(double* vectorS, int* sourceGroupMembers,
		int* targetGroupMembers, int groupIndicator, int length) {
	int i, i_target;
	i_target = 0;

	for (i = 0; i < length; i++) {
		if (vectorS[i] == groupIndicator) {
			targetGroupMembers[i_target] = sourceGroupMembers[i];
			i_target++;
		}
	}
}

/* splitByS helper, updates the group fields*/
static void create_division_group(struct divisionGroup* g, int length,
		struct _spmat* mat, int* sumOfRows, int* groupMembers) {
	g->groupSize = length;
	g->groupSubmatrix = mat;
	g->sumOfRows = sumOfRows;
	g->groupMembers = groupMembers;
}

/*splitByS helper, free old div group*/
static void free_div_group(struct divisionGroup* g, int shouldUseFreell) {

	free(g->sumOfRows);
	free(g->groupMembers);

	if (shouldUseFreell == 1) {
		free_ll((struct _spmat*) g->groupSubmatrix);
	} else {
		/*free submatrix, don't use free_ll, since we don't want to free rows*/
		free_A(g->groupSubmatrix);
	}
	free(g);
}

/* splitByS helper, copies listS to listT*/
static void int_list_copy(int length, int *listT, int *listS) {
	int i;
	for (i = 0; i < length; i++) {
		listT[i] = listS[i];
	}
}

/* TODO: maybe delete splitByS helper, copies nodeS to nodeT*/
void spmat_node_copy(struct spmat_node* nodeT, struct spmat_node* nodeS) {
	nodeT->data = nodeS->data;
	nodeT->index = nodeS->index;
	nodeT->node_name = nodeS->node_name;
	nodeT->next = nodeS->next;
}

/* TODO: maybe delete splitByS helper, copies listS to listT*/
void spmat_node_list_copy(int length, struct spmat_node** listT,
		struct spmat_node** listS) {
	int i;
	for (i = 0; i < length; i++) {
		if (listS[i] != NULL) {
			listT[i] = (struct spmat_node*) malloc(sizeof(struct spmat_node));
			if (listT[i] == NULL)
				panic(ERROR_MALLOC_FAILED);
			spmat_node_copy(listT[i], listS[i]);
		}
	}
}

/* splitByS helper, deep copies spmatS to spmatS*/
void spmat_deep_copy(struct _spmat *spmatT, struct _spmat *spmatS) {
	struct spmat_node **gt_rows, **gs_rows;

	gs_rows = get_private(spmatS);
	gt_rows = get_private(spmatT);

	gt_rows = gs_rows;
	set_private(spmatS, NULL);

	/* copy*/
	spmatT->n = spmatS->n;
	spmatT->add_row = spmatS->add_row;
	spmatT->mult = spmatS->mult;
	spmatT->free = spmatS->free;
	spmatT->private = gt_rows;
}

/* splitByS helper, deep copies gs to gt*/
static void divisionGroup_deep_copy(int gtSize, struct divisionGroup* gt,
		struct divisionGroup* gs) {
	struct _spmat *gtMat, *gsMat;
	int *gtSumOfRows, *gtGroupMembers, *gsSumOfRows, *gsGroupMembers;

	/* allocate spmat*/
	gtMat = spmat_allocate_list_without_rows(gtSize);
	/* allocate sumOfRows*/
	gtSumOfRows = (int*) malloc(gtSize * sizeof(int));
	if (gtSumOfRows == NULL)
		panic(ERROR_MALLOC_FAILED);
	/* allocate groupMembers*/
	gtGroupMembers = (int*) malloc(gtSize * sizeof(int));
	if (gtGroupMembers == NULL)
		panic(ERROR_MALLOC_FAILED);

	gsMat = (struct _spmat*) gs->groupSubmatrix;
	gsSumOfRows = (int*) gs->sumOfRows;
	gsGroupMembers = (int*) gs->groupMembers;

	spmat_deep_copy(gtMat, gsMat);
	int_list_copy(gtSize, gtSumOfRows, gsSumOfRows);
	int_list_copy(gtSize, gtGroupMembers, gsGroupMembers);

	gt->groupSize = gs->groupSize;
	gt->groupMembers = gtGroupMembers;
	gt->groupSubmatrix = gtMat;
	gt->sumOfRows = gtSumOfRows;
}

/* splits g to groups, populates g1 and g2
 * if there's a group of size 0, g1 = g, g2 = NULL */
static struct divisionGroup* split_by_S(double* vectorS, struct divisionGroup* g,
		struct divisionGroup* g1) {
	int i, length, g1Size, g2Size, i1, i2;
	struct _spmat *g1Mat;
	struct _spmat *g2Mat;
	struct spmat_node **gRows, **g1Rows, **g2Rows;
	struct divisionGroup* g2;
	int *gSumOfRows, *g1SumOfRows, *g2SumOfRows;
	int *g1GroupMembers, *g2GroupMembers, *gGroupMembers;

	length = g->groupSize;
	gRows = get_private((struct _spmat*) g->groupSubmatrix);
	gSumOfRows = g->sumOfRows;
	gGroupMembers = g->groupMembers;

	/* if there's a group of size 0, g1 = g, g2 = NULL
	 * in this case, no need to free g*/
	g1Size = calc_size(vectorS, length);

	g2Size = length - g1Size;
	if (g1Size == length || g2Size == length) {
		divisionGroup_deep_copy(length, g1, g);
		g2 = NULL;
		/* free memory*/
		free_div_group(g, 1);
		return g2;
	}

	g2 = (struct divisionGroup*) malloc(sizeof(struct divisionGroup));
	if (g2 == NULL)
		panic(ERROR_MALLOC_FAILED);

	/* allocate spmats*/
	g1Mat = spmat_allocate_list(g1Size);
	g2Mat = spmat_allocate_list(g2Size);
	/* allocate sumOfRows*/
	g1SumOfRows = (int*) malloc(g1Size * sizeof(int));
	if (g1SumOfRows == NULL)
		panic(ERROR_MALLOC_FAILED);
	g2SumOfRows = (int*) malloc(g2Size * sizeof(int));
	if (g2SumOfRows == NULL)
		panic(ERROR_MALLOC_FAILED);
	/* allocate groupMembers*/
	g1GroupMembers = (int*) malloc(g1Size * sizeof(int));
	if (g1GroupMembers == NULL)
		panic(ERROR_MALLOC_FAILED);
	g2GroupMembers = (int*) malloc(g2Size * sizeof(int));
	if (g2GroupMembers == NULL)
		panic(ERROR_MALLOC_FAILED);
	/* allocate rows*/
	g1Rows = g1Mat->private;
	g2Rows = g2Mat->private;

	/* move rows from g to g1, g2*/
	i1 = 0;
	i2 = 0;
	for (i = 0; i < length; i++) {
		if (vectorS[i] == 1) {
			g1Rows[i1] = gRows[i];
			g1SumOfRows[i1] = gSumOfRows[i];
			i1++;
		} else {
			g2Rows[i2] = gRows[i];
			g2SumOfRows[i2] = gSumOfRows[i];
			i2++;
		}
	}

	/* edit rows and sum_of_rows, remove irrelevant nodes*/
	update_mat_rows(vectorS, g1Size, 1, g1Rows, g1SumOfRows);
	update_mat_rows(vectorS, g2Size, -1, g2Rows, g2SumOfRows);
	/* update group_members*/
	update_group_members(vectorS, gGroupMembers, g1GroupMembers, 1, length);
	update_group_members(vectorS, gGroupMembers, g2GroupMembers, -1, length);
	/* update index value=s of the nodes*/
	update_mat_rows_index(g1GroupMembers, g1Size, g1Rows);
	update_mat_rows_index(g2GroupMembers, g2Size, g2Rows);
	/* update spmats*/
	set_private(g1Mat, g1Rows);
	set_private(g2Mat, g2Rows);
	/* create divisionGroups*/
	create_division_group(g1, g1Size, g1Mat, g1SumOfRows,
			g1GroupMembers);
	create_division_group(g2, g2Size, g2Mat, g2SumOfRows,
			g2GroupMembers);
	/* free memory*/
	free_div_group(g, 0);
	return g2;
}

/*
 * creates a new division
 */
static struct division* new_division() {
	struct division* D = (struct division*) malloc(sizeof(struct division));
	if (D == NULL)
		panic(ERROR_MALLOC_FAILED);

	D->len = 0;
	D->divisions = NULL;
	return D;
}

/*  length is number of nodes in graph
 *  Returns a divisonGroup with all nodes in graph
 */
static struct divisionGroup* create_trivial_division(struct graph* inputGraph) {
	int i;
	int length;
	int* groupMembers;
	struct divisionGroup* group = (struct divisionGroup*) malloc(
			sizeof(struct divisionGroup));
	if (group == NULL)
		panic(ERROR_MALLOC_FAILED);

	length = inputGraph->numOfNodes;
	group->groupSize = length;
	group->groupSubmatrix = (inputGraph->A);

	group->sumOfRows = (int*) malloc(length * sizeof(int));
	groupMembers = (int*) malloc(length * sizeof(int));
	if (groupMembers == NULL)
		panic(ERROR_MALLOC_FAILED);

	for (i = 0; i < length; i++) {
		group->sumOfRows[i] = 0;
		groupMembers[i] = i;
	}
	group->groupMembers = groupMembers;

	return group;
}

/*Algorithm 3 as described in sp_project.pdf
 * return O, final division to groups
 * */
struct division* Algorithm3(struct graph* inputGraph) {
	struct divisionGroup* g;
	struct divisionGroup* g1;
	struct divisionGroup* g2;
	double* vectorS;

	struct division* P = new_division();
	struct division* O = new_division();
	add_groupDivision(P, create_trivial_division(inputGraph));

	/*until P is empty:*/
	while (P->len > 0) {

		g1 = (struct divisionGroup*) malloc(sizeof(struct divisionGroup));
		g = remove_first_group(P);

		vectorS = (double*) malloc(g->groupSize * sizeof(double)); /*vectorS is size of group g.*/
		if (vectorS == NULL)
			panic(ERROR_MALLOC_FAILED);

		Algorithm2(vectorS, g, inputGraph);
		modularity_maximization(inputGraph, vectorS, g);

		g2 = split_by_S(vectorS, g, g1);

		if (g2 == NULL) { /*meaning no split was done and group stayed the same*/
			add_groupDivision(O, g1);
		} else {
			update_division_post_split(g1, P, O);
			update_division_post_split(g2, P, O);
		}

		free(vectorS);

	}
	free_division_group(P);
	return O;

}

