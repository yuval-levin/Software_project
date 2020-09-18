#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "modules.h"
#include "spmat.h"
#include "Algorithm3.h"
#include "error_codes.h"
#include <time.h> /*todo: remove time etc*/

void read_row(int i, int n, FILE* input, struct _spmat* A) {
	double* row;
	int k, j, cur;

	row = (double*) malloc(n * sizeof(double));
	if (row == NULL)
		panic(ERROR_MALLOC_FAILED);

	for (j = 0; j < n; j++) {
		k = fread(&cur, sizeof(int), 1, input);
		if (k != 1)
			panic(ERROR_READ_FAILED);
		row[j] = cur;
	}
	add_row_of_size_n(A, row, i, n);
	free(row);
}

void create_graph(FILE* input, struct graph* new_graph) {
	struct _spmat* A;
	long* vector_degrees;
	int i;
	int k;
	int n;
	int cur_deg;
	int deg_sum;
	/*allocating memory*/
	k = fread(&n, sizeof(int), 1, input);
	if (k != 1)
		panic(ERROR_READ_FAILED);
	A = spmat_allocate_list(n);
	vector_degrees = (long*) malloc(n * sizeof(long));
	if (vector_degrees == NULL)
		panic(ERROR_MALLOC_FAILED);

	/*reading input to struct*/
	deg_sum = 0;
	for (i = 0; i < n; i++) {
		k = fread(&cur_deg, sizeof(int), 1, input);
		if (k != 1)
			panic(ERROR_READ_FAILED);
		deg_sum += cur_deg;
		vector_degrees[i] = cur_deg;
		read_row(i, cur_deg, input, A);
	}

	/*initializing graph*/
	new_graph->A = A;
	new_graph->vectorDegrees = vector_degrees;
	new_graph->M = deg_sum;
	new_graph->numOfNodes = n;
}

void write_output_file(struct division* div, FILE* output) {
	int k, n, group_size;
	int* group_members;
	struct node* cur_node;
	struct divisionGroup* cur_group;

	n = div->len;
	k = fwrite(&n, sizeof(int), 1, output);
	if (k != 1)
		panic(ERROR_WRITE_FAILED);

	cur_node = div->divisions;

	/*write groups*/
	while (cur_node != NULL) {
		cur_group = cur_node->data.group;
		group_size = cur_group->groupSize;
		k = fwrite(&group_size, sizeof(int), 1, output);
		if (k != 1)
			panic(ERROR_WRITE_FAILED);

		group_members = cur_group->groupMembers;
		k = fwrite(group_members, sizeof(int), group_size, output);
		if (k != group_size)
			panic(ERROR_WRITE_FAILED);

		cur_node = cur_node->next;
	}
}

void print_result(struct division* div) {
	int n, group_size, i, j;
	int* group_members;
	struct node* cur_node;
	struct divisionGroup* cur_group;

	i = 1;

	n = div->len;
	printf("%s", "num of groups:  ");
	printf("%d \n", n);

	cur_node = div->divisions;

	/*write groups*/
	while (cur_node != NULL) {
		printf("%s", "group number ");
		printf("%d", i);
		cur_group = cur_node->data.group;
		group_size = cur_group->groupSize;
		printf("%s", " of size ");
		printf("%d \n", group_size);
		printf("%s", "group members: ");
		group_members = cur_group->groupMembers;

		for (j = 0; j < group_size; j++) {
			printf("%s", "  ");
			printf("%d", group_members[j]);
		}
		cur_node = cur_node->next;
		i++;
		printf("%s", "\n");
	}
}

int main(int args, char** argv) {
	struct graph* input_graph;
	FILE* input;
	FILE *output;
	struct division* final_division;
	clock_t start, end;

	start = clock();
	srand(time(NULL));

	if (args != 3)
		panic(ERROR_NUM_ARGS);
	input_graph = (struct graph*) malloc(sizeof(struct graph));
	if (input_graph == NULL)
		panic(ERROR_INPUT_NOT_FOUND);

	input = fopen(argv[1], "rb");
	if (input == NULL)
		panic(ERROR_OPEN_FAILED);
	create_graph(input, input_graph);
	fclose(input);
	/*setvbuf (stdout, NULL, _IONBF, 0);*/
	final_division = Algorithm3(input_graph);

	output = fopen(argv[2], "wb");
	if (output == NULL)
		panic(ERROR_OPEN_FAILED);
	write_output_file(final_division, output);
	fclose(output);
	print_result(final_division);

	/*TODO: free final_division and input_graph*/
	freeDivisionGroup(final_division); /*free O and inside*/
	free(input_graph->vectorDegrees);
	free(input_graph);
	end = clock();
	printf("Run took %f seconds\n", ((double) (end - start) / CLOCKS_PER_SEC));
	printf("%s", "done main\n");
	return 0; /*todo: check ok */
}
