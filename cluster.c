#include<stdio.h>
#include <stdlib.h>
#include "modules.h"

// @CR as discussed, should probably come in its own file.
int main(int args, char** argv) {
	graph* graph;

	input = fopen(argv[1], "r");
	assert(input!=NULL);			// TODO: error module
	create_graph(input, graph);
}

void read_row(int i, int n, FILE* input, struct _spmat* A) {
	int* row;
	int k, i, cur;
	for (i = 0; i < n; i++) {
		k = fread(&cur, sizeof(int), 1, input);
		assert(k==1); 				// TODO: error module
		row[i] = cur;
	}
	add_row_of_size_n(A, row, i, n, 1);
}

void create_graph(FILE* input, struct _graph* graph) {
	spmat* A;
	long* vector_degrees;
	int k, n, cur_deg, deg_sum;

	/*allocating memory*/
	k = fread(&n, sizeof(int), 1, input);
	assert(k==1); 					// TODO: error module
	graph = (graph*)malloc(sizeof(graph));
	assert(graph!=NULL)				// TODO: error module
	A = spmat_allocate_list(n);
	vector_degrees = (long*)malloc(n * sizeof(long*));
	assert(vector_degrees!=NULL);	// TODO: error module

	/*reading input to struct*/
	deg_sum = 0;
	for (i = 0; i < n; i++) {
		k = fread(&cur_deg, sizeof(int), 1, input);
		assert(k==1); 				// TODO: error module
		deg_sum += cur_deg;
		vector_degrees[i] = cur_deg;
		read_row(i, cur_deg, input, A);
	}

	/*initializing graph*/
	graph->A = A;
	graph->vectorDegrees = vector_degrees;
	graph->M = deg_sum;
}
