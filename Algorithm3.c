#include <stdio.h>
#include <stdlib.h>
#include "modules.h"
#include "spmat.h"
#include "Algorithm2.h"
#include "ModularityMaximization.h"
#include <assert.h>

/*TODO: is include file.c ok? or should we do headers?*/
/*TODO: add checks for all mallocs.*/
/*adds in start*/
void add_groupDivision(struct division* D, struct divisionGroup* g) {
	struct node* add = (struct node*) malloc(sizeof(struct node));
	if (add == NULL)
		exit(1); /*TODO: print error before exit.*/
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

/* splitByS helper, returns the size of g1*/
int calc_size (double* vectorS, int n) {
	int i, size;
	size = 0;

	for (i = 0; i < n; i++) {
		if (vectorS[i] == 1) {
			size++;
		}
	}
	return size;
}

/* splitByS helper, updates the mat rows,
 * removes irrelevant nodes and updates sum of rows accordingly*/

void update_mat_rows(double* vectorS, int* g_group_members, int num_members, int group_indicator, struct spmat_node** gt_rows, int* gt_sum_of_rows) {
	int i_row, i_group_members;
	struct spmat_node *cur, *prev, *next;
	prev = NULL;

	for (i_row = 0; i_row < num_members; i_row++) {
		cur = gt_rows[i_row];
		i_group_members = 0;
		while (cur != NULL) {
			next = cur->next;
			/* search for the index in vectorS/groupMembers that fits cur*/
			while (i_group_members < num_members && g_group_members[i_group_members] != cur->index)
				i_group_members++;
			/* remove cur, prev remains the same*/
			if (vectorS[i_group_members] != group_indicator) {
				gt_sum_of_rows[i_row]--;
				if (prev == NULL) {
					gt_rows[i_row] = cur->next;
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

/* splitByS helper, updates group members of target, according to group indicator (+-1)*/
void update_group_members(double* vectorS, int* source_group_members, int* target_group_members, int group_indicator, int n) {
	int i, i_target;
	i_target = 0;

	for (i = 0; i < n; i++) {
		if (vectorS[i] == group_indicator) {
			target_group_members[i_target] = source_group_members[i];
			i_target++;
		}
	}
}

/* splitByS helper, updates the group fields*/
void create_division_group(struct divisionGroup* g, int size, struct _spmat* mat, int* sum_of_rows, int* group_members) {
	g->groupSize = size;
	g->groupSubmatrix = mat;
	g->sumOfRows = sum_of_rows;
	g->groupMembers = group_members;
}

/* splitByS helper, free old div group*/
void free_div_group(struct divisionGroup* g) {
	free(g->sumOfRows);
	free(g->groupMembers);

	/* free submatrix, don't use free_ll, since we don't want to free rows*/
	free(get_private((struct _spmat*)g->groupSubmatrix));
	free(g->groupSubmatrix);
	free(g);
}

/* splits g to groups, populates g1 and g2
 * if there's a group of size 0, g1 = g, g2 = NULL */
void splitByS(double* vectorS, struct divisionGroup* g, struct divisionGroup* g1, struct divisionGroup* g2) {
	int i;
	int n;
	int g1_size;
	int g2_size;
	int i1;
	int i2;
	struct _spmat *g1_mat;
	struct _spmat *g2_mat;
	struct spmat_node **g_rows, **g1_rows, **g2_rows;
	int *g_sum_of_rows, *g1_sum_of_rows, *g2_sum_of_rows;
	int *g1_group_members, *g2_group_members, *g_group_members;

	n = g->groupSize;
	g_rows = get_private((struct _spmat*)g->groupSubmatrix);
	g_sum_of_rows = g->sumOfRows;
	g_group_members = g->groupMembers;

	/* if there's a group of size 0, g1 = g, g2 = NULL
	 * in this case, no need to free g*/
	g1_size = calc_size(vectorS, n);
	g2_size = n - g1_size;
	if (g1_size == n || g2_size == n) {
		g1 = g;
		g2 = NULL;
		return;
	}

	/* allocate spmats*/
	g1_mat = spmat_allocate_list(g1_size);
	g2_mat = spmat_allocate_list(g2_size);
	/* allocate sumOfRows*/
	g1_sum_of_rows = (int*)malloc(g1_size * sizeof(int));
	assert(g1_sum_of_rows != NULL);							/* TODO: error module*/
	g2_sum_of_rows = (int*)malloc(g2_size * sizeof(int));
	assert(g2_sum_of_rows != NULL);							/* TODO: error module*/
	/* allocate groupMembers*/
	g1_group_members = (int*)malloc(g1_size * sizeof(int));
	assert(g1_group_members != NULL);						/* TODO: error module*/
	g2_group_members = (int*)malloc(g2_size * sizeof(int));
	assert(g2_group_members != NULL);						/* TODO: error module*/
	/* allocate rows*/
	g1_rows = (struct spmat_node**)malloc(g1_size * sizeof(struct spmat_node*));
	assert(g1_group_members != NULL);						/* TODO: error module*/
	g2_rows = (struct spmat_node**)malloc(g2_size * sizeof(struct spmat_node*));
	assert(g2_group_members != NULL);						/* TODO: error module*/

	/* move rows from g to g1, g2*/
	i1 = 0;
	i2 = 0;
	for (i = 0; i < n; i++) {
		if (vectorS[i] == 1) {
			g1_rows[i1] = g_rows[i];
			g1_sum_of_rows[i1] = g_sum_of_rows[i];
			i1++;
		} else {
			g2_rows[i2] = g_rows[i];
			g2_sum_of_rows[i2] = g_sum_of_rows[i];
			i2++;
		}
	}

	/* edit rows and sum_of_rows, remove irrelevant nodes*/
	update_mat_rows(vectorS, g_group_members, g1_size, 1, g1_rows, g1_sum_of_rows);
	update_mat_rows(vectorS, g_group_members, g2_size, -1, g2_rows, g2_sum_of_rows);

	/* update spmats*/
	set_private(g1_mat, g1_rows);
	set_private(g2_mat, g2_rows);

	/* update group_members*/
	update_group_members(vectorS, g_group_members, g1_group_members, 1, n);
	update_group_members(vectorS, g_group_members, g2_group_members, -1, n);

	/* create divisionGroups*/
	create_division_group(g1, g1_size, g1_mat, g1_sum_of_rows, g1_group_members);
	create_division_group(g2, g2_size, g2_mat, g2_sum_of_rows, g2_group_members);

	/* free memory*/
	free_div_group(g);
}

struct division* new_division() {
	struct division* D = (struct division*) malloc(sizeof(struct division));
	if (D == NULL)
		exit(1);  /*TODO: print error before exit.*/
	D->len = 0;
	D->divisions = NULL;
	return D;
}


/* PARAMS : n is number of nodes in graph
 * DESC : Returns a divisonGroup with all nodes in graph
 */
struct divisionGroup* createTrivialDivision(struct graph* inputGraph) {
	int i;
	int n;
	int* group_members;
	struct divisionGroup* group = (struct divisionGroup*)malloc(sizeof(struct divisionGroup));
	if (group == NULL)
		exit(1); /* TODO: error module*/
	n = inputGraph -> numOfNodes;
	group->groupSize = n;
	group->groupSubmatrix = (inputGraph->A);
	group->sumOfRows = (int*) malloc(n * sizeof(int));
	group_members = (int*)malloc(n * sizeof(int));
	assert(group_members != NULL);						/* TODO: error module*/
	for (i = 0; i < n; i++) {
		group->sumOfRows[i] = 0;
		group_members[i] = i;
	}
	group->groupMembers = group_members;

	return group;
}


struct division* Algorithm3(struct graph* inputGraph) {
	struct divisionGroup* g;
	struct divisionGroup* g1;
	struct divisionGroup* g2;
	double* vectorS;
	struct division* P = new_division();
	struct division* O = new_division();
	add_groupDivision(P, createTrivialDivision( inputGraph));
	/* TODO: calc sum of rows*/

	g1 = (struct divisionGroup*)malloc(sizeof(struct divisionGroup));
	g2 = (struct divisionGroup*)malloc(sizeof(struct divisionGroup));
	while (P->len > 0) {
		g = removeFirstGroup(P);
		vectorS = (double*) malloc(g->groupSize * sizeof(double)); /*vectorS is size of group g.*/
		if (vectorS == NULL)
			exit(1); /*TODO: print error before exit.*/
		Algorithm2(vectorS, g,inputGraph);
		printf("%s", "d\n");
		modularityMaximization(inputGraph,vectorS, g);
		printf("%s", "e\n");
		splitByS(vectorS, g, g1, g2);
		printf("%s", "f\n");

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

