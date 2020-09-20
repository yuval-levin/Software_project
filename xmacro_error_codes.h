#ifndef ERROR_CODE
#error "xmacro_error_codes.h included without defining ERROR_CODE"
#endif

ERROR_CODE(ERROR_INPUT_NOT_FOUND, "ERROR: Input File not found")
ERROR_CODE(ERROR_MALLOC_FAILED, "ERROR: Memory allocation failed")
ERROR_CODE(ERROR_DIVISION_BY_ZERO, "ERROR: Invalid division by zero occured")
ERROR_CODE(ERROR_READ_FAILED, "ERROR: Reading from file failed")
ERROR_CODE(ERROR_WRITE_FAILED, "ERROR: Writing to file failed")
ERROR_CODE(ERROR_OPEN_FAILED, "ERROR: File opening failed")
ERROR_CODE(ERROR_NUM_ARGS, "ERROR: Cluster expected 3 Arguments, but received a different number")
ERROR_CODE(ERROR_LOOP_LIMIT_REACHED, "ERROR: Power iteration has entered infinite loop")
