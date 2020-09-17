#include <stdio.h>
#include <stdlib.h>
#include "modules.h"
#include "spmat.h"
#include "Algorithm2.h"
#include "ModularityMaximization.h"
#include "error_codes.h"

/*helper function to free a linked list*/
void freeList(struct spmat_node* cur) {
	struct spmat_node* next;
	while (cur != NULL) {
		next = cur->next;
		free(cur);
		cur = next;
	}
}

/*helper function to free the division groups*/
void freeGroup(struct divisionGroup* group) {
	int i, n;
	struct spmat_node** rows;

	rows = group->groupSubmatrix->private;
	n = group->groupSubmatrix->n;

	for (i = 0; i < n; i++) {
		freeList(rows[i]);
	}

	free(group->sumOfRows);
	free(group->groupMembers);
	free(group->groupSubmatrix);
	free(group); /*TODO: delete?*/
}

/*helper function to free the  divisions given to us by Algorithm 3
 * Which are inside of O */
void freeDivisionGroup(struct division* O) {
	int i, n;
	struct node* currDiv;
	struct node* nextDiv;
	n = O->len;

	currDiv = O->divisions;
	for (i = 0; i < n; i++) {
		nextDiv = currDiv->next;
		freeGroup(currDiv->data.group);
		free(currDiv);
		currDiv = nextDiv;
	}
	free(O);
}
/*adds in start*/
void add_groupDivision(struct division* D, struct divisionGroup* g) {
	struct node* add = (struct node*) malloc(sizeof(struct node));
	if (add == NULL)
		panic(ERROR_MALLOC_FAILED);

	add->data.group = g;
	add->next = NULL;
	if (D->len == 0)
		D->divisions = add;
	else {
		add->next = D->divisions;
		D->divisions = add; /*todo make sure no cycle*/
	}

	D->len = (D->len) + 1;
}

/* removes  first group from D and returns it*/
struct divisionGroup* removeFirstGroup(struct division* D) {
	struct node* group;
	struct node* nextGroup;

	group = D->divisions;
	nextGroup = group->next;
	group->next = NULL;
	D->divisions = nextGroup;
	D->len = (D->len) - 1;
	return group->data.group;
}

/*helper function for Algorithm 2.
 * After split,if group (here represented by g as input argument, but is  g1 or g2 in algo2) larger than 1 are put in P,
 * Any groups of size 1 are put in O */
void updateDivisionPostSplit(struct divisionGroup* g, struct division* P,
		struct division* O) {
	if (g->groupSize == 1)
		add_groupDivision(O, g);
	else
		add_groupDivision(P, g);

}

/* splitByS helper, returns the size of g1*/
int calc_size(double* vectorS, int n) {
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
void update_mat_rows(double* vectorS, int num_members, int group_indicator,
		struct spmat_node** gt_rows, int* gt_sum_of_rows) {
	int i_row;
	struct spmat_node *cur, *prev, *next;

	/* delete irrelevant nodes*/
	for (i_row = 0; i_row < num_members; i_row++) {
		prev = NULL;
		cur = gt_rows[i_row];
		while (cur != NULL) {
			next = cur->next;

			/* remove cur, prev remains the same*/
			if (vectorS[cur->index] != group_indicator) {
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

/* splitByS helper, updates nodes index according to the new mat */
void update_mat_rows_index(int* g_group_members, int num_members,
		struct spmat_node** g_rows) {
	int i_row, count, i;
	int* map_name_to_col_index;
	struct spmat_node *cur;

	/* create a map, will be used to update index field to fit the new column index in the new spmat*/
	count = 0;
	map_name_to_col_index = (int*) malloc(
			(g_group_members[num_members - 1] + 1) * sizeof(int));
	if (map_name_to_col_index == NULL)
		panic(ERROR_MALLOC_FAILED);

	for (i = 0; i < num_members; i++) {
		map_name_to_col_index[g_group_members[i]] = count;
		count++;
	}

	/* delete irrelevant nodes*/
	for (i_row = 0; i_row < num_members; i_row++) {
		cur = g_rows[i_row];
		while (cur != NULL) {
			cur->index = map_name_to_col_index[cur->node_name];
			cur = cur->next;
		}
	}
	free(map_name_to_col_index);
}

/* splitByS helper, updates group members of target, according to group indicator (+-1)*/
void update_group_members(double* vectorS, int* source_group_members,
		int* target_group_members, int group_indicator, int n) {
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
void create_division_group(struct divisionGroup* g, int size,
		struct _spmat* mat, int* sum_of_rows, int* group_members) {
	g->groupSize = size;
	g->groupSubmatrix = mat;
	g->sumOfRows = sum_of_rows;
	g->groupMembers = group_members;
}

/*splitByS helper, free old div group*/
void free_div_group(struct divisionGroup* g, int should_use_free_ll) {

	free(g->sumOfRows);
	free(g->groupMembers);

	if (should_use_free_ll == 1) {
		free_ll((struct _spmat*) g->groupSubmatrix);
	} else {
		/*free submatrix, don't use free_ll, since we don't want to free rows*/
		free(get_private((struct _spmat*) g->groupSubmatrix));
		free(g->groupSubmatrix);
	}
	free(g);
}

/* splitByS helper, copies list_s to list_t*/
void int_list_copy(int size, int *list_t, int *list_s) {
	int i;
	for (i = 0; i < size; i++) {
		list_t[i] = list_s[i];
	}
}

/* TODO: maybe delete splitByS helper, copies node_s to node_t*/
void spmat_node_copy(struct spmat_node* node_t, struct spmat_node* node_s) {
	node_t->data = node_s->data;
	node_t->index = node_s->index;
	node_t->node_name = node_s->node_name;
	node_t->next = node_s->next;
}

/* TODO: maybe delete splitByS helper, copies list_s to list_t*/
void spmat_node_list_copy(int size, struct spmat_node** list_t,
		struct spmat_node** list_s) {
	int i;
	for (i = 0; i < size; i++) {
		if (list_s[i] != NULL) {
			list_t[i] = (struct spmat_node*) malloc(sizeof(struct spmat_node));
			if (list_t[i] == NULL)
				panic(ERROR_MALLOC_FAILED);
			spmat_node_copy(list_t[i], list_s[i]);
		}
	}
}

/* splitByS helper, deep copies spmat_s to spmat_t*/
void spmat_deep_copy(struct _spmat *spmat_t, struct _spmat *spmat_s) {
	struct spmat_node **gt_rows, **gs_rows;

	gs_rows = get_private(spmat_s);
	gt_rows = get_private(spmat_t);

	gt_rows = gs_rows;
	set_private(spmat_s, NULL);

	/* copy*/
	spmat_t->n = spmat_s->n;
	spmat_t->add_row = spmat_s->add_row;
	spmat_t->mult = spmat_s->mult;
	spmat_t->free = spmat_s->free;
	spmat_t->private = gt_rows;
}

/* splitByS helper, deep copies gs to gt*/
void divisionGroup_deep_copy(int gt_size, struct divisionGroup* gt,
		struct divisionGroup* gs) {
	struct _spmat *gt_mat, *gs_mat;
	int *gt_sum_of_rows, *gt_group_members, *gs_sum_of_rows, *gs_group_members;

	/* allocate spmat*/
	gt_mat = spmat_allocate_list(gt_size);
	/* allocate sumOfRows*/
	gt_sum_of_rows = (int*) malloc(gt_size * sizeof(int));
	if (gt_sum_of_rows == NULL)
		panic(ERROR_MALLOC_FAILED);
	/* allocate groupMembers*/
	gt_group_members = (int*) malloc(gt_size * sizeof(int));
	if (gt_group_members == NULL)
		panic(ERROR_MALLOC_FAILED);

	gs_mat = (struct _spmat*) gs->groupSubmatrix;
	gs_sum_of_rows = (int*) gs->sumOfRows;
	gs_group_members = (int*) gs->groupMembers;

	spmat_deep_copy(gt_mat, gs_mat);
	int_list_copy(gt_size, gt_sum_of_rows, gs_sum_of_rows);
	int_list_copy(gt_size, gt_group_members, gs_group_members);

	gt->groupSize = gs->groupSize;
	gt->groupMembers = gt_group_members;
	gt->groupSubmatrix = gt_mat;
	gt->sumOfRows = gt_sum_of_rows;
}

/* splits g to groups, populates g1 and g2
 * if there's a group of size 0, g1 = g, g2 = NULL */
struct divisionGroup* splitByS(double* vectorS, struct divisionGroup* g,
		struct divisionGroup* g1) {
	int i, n, g1_size, g2_size, i1, i2;
	struct _spmat *g1_mat;
	struct _spmat *g2_mat;
	struct spmat_node **g_rows, **g1_rows, **g2_rows;
	struct divisionGroup* g2;
	int *g_sum_of_rows, *g1_sum_of_rows, *g2_sum_of_rows;
	int *g1_group_members, *g2_group_members, *g_group_members;

	n = g->groupSize;
	g_rows = get_private((struct _spmat*) g->groupSubmatrix);
	g_sum_of_rows = g->sumOfRows;
	g_group_members = g->groupMembers;

	/* if there's a group of size 0, g1 = g, g2 = NULL
	 * in this case, no need to free g*/
	g1_size = calc_size(vectorS, n);

	g2_size = n - g1_size;
	if (g1_size == n || g2_size == n) {
		divisionGroup_deep_copy(n, g1, g);
		g2 = NULL;
		/* free memory*/
		free_div_group(g, 1);
		return g2;
	}

	g2 = (struct divisionGroup*) malloc(sizeof(struct divisionGroup));
	if (g2 == NULL)
		panic(ERROR_MALLOC_FAILED);

	/* allocate spmats*/
	g1_mat = spmat_allocate_list(g1_size);
	g2_mat = spmat_allocate_list(g2_size);
	/* allocate sumOfRows*/
	g1_sum_of_rows = (int*) malloc(g1_size * sizeof(int));
	if (g1_sum_of_rows == NULL)
		panic(ERROR_MALLOC_FAILED);
	g2_sum_of_rows = (int*) malloc(g2_size * sizeof(int));
	if (g2_sum_of_rows == NULL)
		panic(ERROR_MALLOC_FAILED);
	/* allocate groupMembers*/
	g1_group_members = (int*) malloc(g1_size * sizeof(int));
	if (g1_group_members == NULL)
		panic(ERROR_MALLOC_FAILED);
	g2_group_members = (int*) malloc(g2_size * sizeof(int));
	if (g2_group_members == NULL)
		panic(ERROR_MALLOC_FAILED);
	/* allocate rows*/
	g1_rows = g1_mat->private;
	g2_rows = g2_mat->private;

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
	update_mat_rows(vectorS, g1_size, 1, g1_rows, g1_sum_of_rows);
	update_mat_rows(vectorS, g2_size, -1, g2_rows, g2_sum_of_rows);
	/* update group_members*/
	update_group_members(vectorS, g_group_members, g1_group_members, 1, n);
	update_group_members(vectorS, g_group_members, g2_group_members, -1, n);
	/* update index value=s of the nodes*/
	update_mat_rows_index(g1_group_members, g1_size, g1_rows);
	update_mat_rows_index(g2_group_members, g2_size, g2_rows);
	/* update spmats*/
	set_private(g1_mat, g1_rows);
	set_private(g2_mat, g2_rows);
	/* create divisionGroups*/
	create_division_group(g1, g1_size, g1_mat, g1_sum_of_rows,
			g1_group_members);
	create_division_group(g2, g2_size, g2_mat, g2_sum_of_rows,
			g2_group_members);
	/* free memory*/
	free_div_group(g, 0);

	return g2;
}

struct division* new_division() {
	struct division* D = (struct division*) malloc(sizeof(struct division));
	if (D == NULL)
		panic(ERROR_MALLOC_FAILED);

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
	struct divisionGroup* group = (struct divisionGroup*) malloc(
			sizeof(struct divisionGroup));
	if (group == NULL)
		panic(ERROR_MALLOC_FAILED);

	n = inputGraph->numOfNodes;
	group->groupSize = n;
	group->groupSubmatrix = (inputGraph->A);
	group->sumOfRows = (int*) malloc(n * sizeof(int));
	group_members = (int*) malloc(n * sizeof(int));
	if (group_members == NULL)
		panic(ERROR_MALLOC_FAILED);

	for (i = 0; i < n; i++) {
		group->sumOfRows[i] = 0;
		group_members[i] = i;
	}
	group->groupMembers = group_members;

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
	add_groupDivision(P, createTrivialDivision(inputGraph));

	while (P->len > 0) {

		g1 = (struct divisionGroup*) malloc(sizeof(struct divisionGroup));
		g = removeFirstGroup(P);

		vectorS = (double*) malloc(g->groupSize * sizeof(double)); /*vectorS is size of group g.*/
		if (vectorS == NULL)
			panic(ERROR_MALLOC_FAILED);

		Algorithm2(vectorS, g, inputGraph);
		printf("%s", "before mod max \n");
		modularityMaximization(inputGraph, vectorS, g);

		g2 = splitByS(vectorS, g, g1);

		if (g2 == NULL) { /*meaning no split was done and group stayed the same*/
			add_groupDivision(O, g1);
		} else {
			updateDivisionPostSplit(g1, P, O);
			updateDivisionPostSplit(g2, P, O);
		}

		free(vectorS);

	}
	freeDivisionGroup(P);
	return O;

}

