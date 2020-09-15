#ifndef _SPMAT_H
#define _SPMAT_H

/*
 * spmat_node structure which implements the linked list.
 * index contains the column index.
 */
typedef struct spmat_node {
	int data;
	int index;
	int node_name;
	struct spmat_node *next;
} spmat_node;

typedef struct _spmat {
	/* Matrix size (n*n) */
	int		n;

	/* Adds row i the matrix. Called before any other call,
	 * exactly n times in order (i = 0 to n-1) */
	void	(*add_row)(struct _spmat *A, const double *row, int i);

	/* Frees all resources used by A */
	void	(*free)(struct _spmat *A);

	/* Multiplies matrix A by vector v, into result (result is pre-allocated) */
	void	(*mult)(const struct _spmat *A, const double *v, double *result);

	/* Private field for inner implementation.
	 * Should not be read or modified externally */
	void	*private;
} spmat;

/* Allocates a new linked-lists sparse matrix of size n */
spmat* spmat_allocate_list(int n);
void mult_ll(const struct _spmat *A, const double *v, double *result);
void add_row_of_size_n(struct _spmat *A, const double *row, int i, int n);
struct spmat_node** get_private(struct _spmat* mat);
void set_private(struct _spmat* mat, struct spmat_node** rows);
void free_ll(struct _spmat *A);

#endif
