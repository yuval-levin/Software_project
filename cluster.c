#include <stdio.h>
#include <stdlib.h>
#include "modules.h"
#include "spmat.h"
#include "Algorithm3.h"
#include "error_codes.h"

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
	double* vectorDegrees;
	double* degreesDividedByM;
	int i;
	int k;
	int n;
	int curDeg;
	int degSum;
	/*allocating memory*/
	k = fread(&n, sizeof(int), 1, input);
	if (k != 1)
		panic(ERROR_READ_FAILED);
	A = spmat_allocate_list(n);
	vectorDegrees = (double*) malloc(n * sizeof(double));
	if (vectorDegrees == NULL)
		panic(ERROR_MALLOC_FAILED);

	degreesDividedByM = (double*) malloc(n * sizeof(double));
	if (degreesDividedByM == NULL)
			panic(ERROR_MALLOC_FAILED);
	/*reading input to struct*/
	degSum = 0;
	for (i = 0; i < n; i++) {
		k = fread(&curDeg, sizeof(int), 1, input);
		if (k != 1)
			panic(ERROR_READ_FAILED);
		degSum += curDeg;
		vectorDegrees[i] = curDeg;
		read_row(i, curDeg, input, A);
	}
	/*loop again, to calc degreesDividedByM. it Needs M, which is only calculated after prior loop*/
		for (i = 0; i < n; i++) {
			if(degSum < epsilon) panic(ERROR_DIVISION_BY_ZERO); /*when M of graph is zero */
			degreesDividedByM[i] = vectorDegrees[i] / degSum; }

		/*initializing graph*/
		new_graph->A = A;
		new_graph->vectorDegrees = vectorDegrees;
		new_graph->M = degSum;
		new_graph->degreesDividedByM = degreesDividedByM;
		new_graph->numOfNodes = n;
}

void write_output_file(struct division* div, FILE* output) {
	int k, length, groupSize;
	int* groupMembers;
	struct node* curNode;
	struct divisionGroup* curGroup;

	length = div->len;
	k = fwrite(&length, sizeof(int), 1, output);
	if (k != 1)
		panic(ERROR_WRITE_FAILED);

	curNode = div->divisions;

	/*write groups*/
	while (curNode != NULL) {
		curGroup = curNode->data.group;
		groupSize = curGroup->groupSize;
		k = fwrite(&groupSize, sizeof(int), 1, output);
		if (k != 1)
			panic(ERROR_WRITE_FAILED);

		groupMembers = curGroup->groupMembers;
		k = fwrite(groupMembers, sizeof(int), groupSize, output);
		if (k != groupSize)
			panic(ERROR_WRITE_FAILED);

		curNode = curNode->next;
	}
}


int main(int args, char** argv) {
	struct graph* inputGraph;
	FILE* input;
	FILE *output;
	struct division* finalDivision;

	if (args != 3)
		panic(ERROR_NUM_ARGS);
	inputGraph = (struct graph*) malloc(sizeof(struct graph));
	if (inputGraph == NULL)
		panic(ERROR_INPUT_NOT_FOUND);

	input = fopen(argv[1], "rb");
	if (input == NULL)
		panic(ERROR_OPEN_FAILED);
	create_graph(input, inputGraph);
	fclose(input);
	/*setvbuf (stdout, NULL, _IONBF, 0);*/
	finalDivision = Algorithm3(inputGraph);

	output = fopen(argv[2], "wb");
	if (output == NULL)
		panic(ERROR_OPEN_FAILED);

	write_output_file(finalDivision, output);
	/*print_result(final_division);*/
	fclose(output);

	free_division_group(finalDivision); /*free O and inside*/
	free(inputGraph->vectorDegrees);
	free(inputGraph->degreesDividedByM);
	free(inputGraph);

	printf("%s", "done main\n");
	return 0;
}
