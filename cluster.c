#include <stdio.h>
#include <stdlib.h>
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

void create_graph(FILE* input, struct graph* newGraph) {
	struct _spmat* A;
	double* vectorDegrees;
	double* degreesDividedByM;
	int i, k, n, curDeg, degSum;
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
		newGraph->A = A;
		newGraph->vectorDegrees = vectorDegrees;
		newGraph->M = degSum;
		newGraph->degreesDividedByM = degreesDividedByM;
		newGraph->numOfNodes = n;
}

void write_output_file(struct division* div, FILE* output) {
	int k, n, groupSize;
	int* groupMembers;
	struct node* curNode;
	struct divisionGroup* curGroup;

	n = div->len;
	k = fwrite(&n, sizeof(int), 1, output);
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

void print_result(struct division* div) {
	int n, groupSize, i, j;
	int* groupMembers;
	struct node* curNode;
	struct divisionGroup* curGroup;

	i = 1;

	n = div->len;
	printf("%s", "num of groups:  ");
	printf("%d \n", n);

	curNode = div->divisions;

	/*write groups*/
	while (curNode != NULL) {
		printf("%s", "group number ");
		printf("%d", i);
		curGroup = curNode->data.group;
		groupSize = curGroup->groupSize;
		printf("%s", " of size ");
		printf("%d \n", groupSize);
		printf("%s", "group members: ");
		groupMembers = curGroup->groupMembers;

		for (j = 0; j < groupSize; j++) {
			printf("%s", "  ");
			printf("%d", groupMembers[j]);
		}
		curNode = curNode->next;
		i++;
		printf("%s", "\n");
	}
}

int main(int args, char** argv) {
	struct graph* inputGraph;
	FILE* input;
	FILE *output;
	struct division* finalDivision;
	clock_t start, end;

	start = clock();
	srand(time(NULL));

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
	print_result(finalDivision);
	fclose(output);

	free_division_group(finalDivision); /*free O and inside*/
	free(inputGraph->vectorDegrees);
	free(inputGraph->degreesDividedByM);
	free(inputGraph);

	end = clock();
	printf("Run took %f seconds\n", ((double) (end - start) / CLOCKS_PER_SEC));
	printf("%s", "done main\n");
	return 0;
}
