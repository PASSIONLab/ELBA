#ifndef _MPI_UTILS_H_
#define _MPI_UTILS_H_

#include <stdio.h>
#include <mpi.h>

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1
#endif

#define CHECK_MPI(fn) { int err; err = (fn); if (err != MPI_SUCCESS) handle_mpi_error(err, __FILE__, __func__, __LINE__, 1, #fn); } /* calls MPI_Abort on an error */
#define CHECK_MPI_FAILED(fn) handle_mpi_error((fn), __FILE__, __func__, __LINE__, 0) /* returns true if an error */

static int handle_mpi_error(int err, const char *file, const char *func, int line, int abort, const char *code) {
	char msg[MPI_MAX_ERROR_STRING];
	int resultlen, errorClass, rank;
	if (err == MPI_SUCCESS) return 0;
	msg[0] = '\0';
	MPI_Error_string(err, msg, &resultlen);
	errorClass = -1;
	MPI_Error_class(err, &errorClass);
	rank = -1;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	fprintf(stderr, "Rank %d: MPI_Error! (%d) Class %d '%s': %s:%d (%s) %s\n", rank, err, errorClass, msg, file, line, func, code);
	if (abort)
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	return (err ? err : 1);
}

#endif