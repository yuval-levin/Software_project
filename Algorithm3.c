#include<stdio.h>
#include <stdlib.h>
#include "modules.h"
#include "Algorithm2.h"
#include "ModularityMaximization.h"

//TODO: is include file.c ok? or should we do headers?
//TODO: add checks for all mallocs.
//adds in start
void add_groupDivision(struct division* D, struct divisionGroup* g) {
	struct node* add = (struct node*) malloc(sizeof(struct node));
	if (add == NULL)
		exit(1); //TODO: print error before exit.
	add->data.group = g;
	add->next = NULL;
	if (D->len == 0)
		D->divisions = add;
	else
		add->next = D->divisions;

	D->len = (D->len) + 1;
}

/* removes  first group from D and returns it*/
struct divisionGroup* removeFirstGroup(struct division* D) {
	struct node* group = D->divisions;
	struct node* nextGroup = group->next;

	group->next = NULL;
	D->divisions = nextGroup;
	D->len = (D->len) - 1;

	return group->data.group;
}

void updateDivisionPostSplit(struct divisionGroup* g, struct division* P,
		struct division* O) {
	if (g->groupSize == 1)
		add_groupDivision(O, g);
	else
		add_groupDivision(P, g);

}

/*splitByS helper, returns the size of g1*/
int calc_size (int* vectorS, int n) {
	int i, size;
	size = 0;

	for (i = 0; i < n; i++) {
		if (vectorS[i] == 1) {
			size++;
		}
	}
	return size;
}

/*splits g to groups, populating g1 and g2
 *if there's a group of size 0, g1 = g, g2 = NULL */
void splitByS(int* vectorS, struct divisionGroup* g, struct divisionGroup* g1, struct divisionGroup* g2) {
	int i, n, g1_size, g2_size, i1, i2;
	struct _spmat *g1_mat, *g2_mat;
	struct spmat_node **g_rows, **g1_rows, **g2_rows;
	int *g_sum_of_rows, *g1_sum_of_rows, *g2_sum_of_rows;

	n = g->groupSize;
	g_rows = get_private(g->groupSubmatrix);
	g_sum_of_rows = g->sumOfRows;

	/*if there's a group of size 0, g1 = g, g2 = NULL
	 *in this case, no need to free g*/
	g1_size = calc_size(vectorS, n);
	g2_size = n - g1_size;
	if (g1_size == n || g2_size == n) {
		g1 = g;
		g2 = NULL;
		return;
	}

	/*allocate spmats*/
	g1_mat = spmat_allocate_list(g1_size);
	g2_mat = spmat_allocate_list(g2_size);

	/*move rows from g to g1, g2*/
	i1 = 0;
	i2 = 0;
	for (i = 0; i < n; i++) {
		if (vectorS[i] == 1) {
			g1_rows[i1] = g_rows[i];
			g1_sum_of_rows[i1] == g_sum_of_rows[i];
			i1++;
		} else {
			g2_rows[i2] = g_rows[i];
			g2_sum_of_rows[i2] == g_sum_of_rows[i];
			i2++;
		}
	}

	/*edit rows and sum_of_rows, remove irrelevant nodes*/
	//TODO: edit rows and sum_of_rows

	/*update spmats*/
	set_private(g1_mat, g1_rows);
	set_private(g2_mat, g2_rows);

	/*create divisionGroups*/
	//TODO: write create_division_group
	create_division_group(g1_size, g1_mat, g1_sum_of_rows);
	create_division_group(g2_size, g2_mat, g2_sum_of_rows);

	/*free unnecessary memory*/
	//TODO: write free_div_group
	free_div_group(g);

}
struct division* new_division() {
	struct division* D = (struct division*) malloc(sizeof(struct division));
	if (D == NULL)
		exit(1); //TODO: print error before exit.
	D->len = 0;
	D->divisions = NULL;
}


/*PARAMS : n is number of nodes in graph
 * DESC : Returns a divisonGroup with all nodes in graph
 */
struct divisionGroup* createTrivialDivision(int n, struct graph* inputGraph) {
	int i;
	struct divisionGroup* group = (struct divisionGroup*)malloc(sizeof(struct divisionGroup));
	if (group == NULL)
		exit(1); //TODO: print error before exit.
	group->groupSize = n;
	group->groupSubmatrix = &(inputGraph->A);
	group->sumOfRows = (int*) malloc(n * sizeof(int));
	for (i = 0; i < n; i++) {
		group->sumOfRows[i] = 0;
	}
	return group;
}


struct division* Algorithm3(int numOfNodes, struct graph inputGraph) {
	struct divisionGroup* g;
	struct divisionGroup* g1;
	struct divisionGroup* g2;
	int* vectorS;
	struct division* P = new_division();
	struct division* O = new_division();
	add_groupDivision(P, createTrivialDivision(numOfNodes, &inputGraph));

	while (P->len > 0) {
		g = removeFirstGroup(P);
		vectorS = (int*) malloc(g->groupSize * sizeof(int)); //vectorS is size of group g.
		if (vectorS == NULL)
			exit(1); //TODO: print error before exit.
		Algorithm2(vectorS, g,&inputGraph);
		ModularityMaximization(vectorS, g);
		splitByS(vectorS, g1, g2);

		if (g2 == NULL)
			add_groupDivision(O, g1);
		else {
			updateDivisionPostSplit(g1, P, O);
			updateDivisionPostSplit(g2, P, O);
		}
		free(vectorS);
	}
	return O;

}

