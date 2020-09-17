/*
 * spmat.c
 * how to create a sparse matrix:
 * spmat* mat;
 * mat = spmat_allocate_list(n);
 * create_matrix(mat, input);
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "spmat.h"
#include "error_codes.h"

void add_row_ll(struct _spmat *A, const double *row, int i);
void free_ll(struct _spmat *A);

double sumTimesVectorHelper(struct spmat_node* cur_node, const double *v) {
	int index = 0;
	double sum = 0;
	struct spmat_node* cur;
	cur = cur_node;
	while (cur != NULL) {
		index = cur->index;
		sum += (cur->data) * (v[index]);
		cur = (struct spmat_node*) cur->next;
	}
	return sum;
}

/*
 * creating the matrix one row at a time
 */
void create_matrix(struct _spmat *A, FILE* input) {
	double *row;
	int i, k, n;

	n = A->n;
	fseek(input, 8, SEEK_SET); /*skipping dimensions*/
	row = (double*) malloc(n * sizeof(double));
	if (row == NULL)
		panic(ERROR_MALLOC_FAILED);

	for (i = 0; i < n; i++) {
		k = fread(row, sizeof(double), n, input);
		if (k != n)
			panic(ERROR_READ_FAILED);

		A->add_row(A, row, i);
	}

	free(row);
	rewind(input);
}

/*
 * mult implementation for linked list
 */
void mult_ll(const struct _spmat *A, const double *v, double *result) {
	int row_ind;
	double sum;
	struct spmat_node* cur_node;
	struct spmat_node** rows = (struct spmat_node**) A->private;

	/*calculate the result of each row*/
	for (row_ind = 0; row_ind < A->n; row_ind++) {
		cur_node = rows[row_ind];

		/*sum all the relevant multiplications*/
		sum = sumTimesVectorHelper(cur_node, v);
		result[row_ind] = sum;
	}
}

/*
 * helper for add_row_ll
 */
void add_row_of_size_n(struct _spmat *A, const double *row, int i, int n) {
	int first, j;
	struct spmat_node *cur, *prev;
	first = 0;

	if (n == 0) {
		((struct spmat_node**) A->private)[i] = NULL;
	}

	for (j = 0; j < n; j++) {
		cur = (struct spmat_node*) malloc(sizeof(struct spmat_node));
		if (cur == NULL)
			panic(ERROR_MALLOC_FAILED);

		cur->data = 1;
		cur->index = row[j];
		cur->node_name = row[j];
		cur->next = NULL;

		/*first struct spmat_node*/
		if (first == 0) {
			((struct spmat_node**) A->private)[i] = cur; /*inserted as the 1st node of the ith row*/
			first = 1; /*no updating prev.next*/
		} else {
			prev->next = (struct spmat_node*) cur; /*cur is not the 1st, update prev.next*/
		}
		prev = cur;
	}
}

/*
 * add_row implementation for linked list
 */
void add_row_ll(struct _spmat *A, const double *row, int i) {
	add_row_of_size_n(A, row, i, A->n);
}

struct spmat_node** get_private(struct _spmat* mat) {
	return mat->private;
}

void set_private(struct _spmat* mat, struct spmat_node** rows) {
	mat->private = rows;
}

/*
 * helper for free_ll, frees the nodes of the given row recursively
 */
void free_row_ll(struct spmat_node* row) {
	if (row != NULL) {
		free_row_ll((struct spmat_node*) row->next);
		free(row);
	}
}

/*
 * free implementation for linked list
 */
void free_ll(struct _spmat *A) {
	int i;
	struct spmat_node** rows;
	rows = A->private;

	if (rows == NULL)
		return;

	for (i = 0; i < A->n; i++) {
		free_row_ll(rows[i]);
	}
	free(A->private);
	free(A);
}

/*
 * allocating memory for rows.
 * populate the functions with linked list implementation.
 */
spmat* spmat_allocate_list(int n) {
	spmat* mat;
	struct spmat_node** rows;

	/*allocating memory*/
	mat = (spmat*) malloc(sizeof(spmat));
	if (mat == NULL)
		panic(ERROR_MALLOC_FAILED);

	rows = (struct spmat_node**) malloc(n * sizeof(struct spmat_node*));
	if (rows == NULL)
		panic(ERROR_MALLOC_FAILED);

	/*initializing*/
	mat->n = n;
	mat->add_row = &add_row_ll;
	mat->free = &free_ll;
	mat->mult = &mult_ll;
	mat->private = rows;

	return mat;
}
