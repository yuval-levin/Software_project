#ifndef ERROR_CODE
#error "xmacro_error_codes.h included without defining ERROR_CODE"
#endif

ERROR_CODE(ERROR_INPUT_NOT_FOUND, "Input File not found")
ERROR_CODE(ERROR_MALLOC_FAILED, "Memory allocation failed")
ERROR_CODE(ERROR_DIVISION_BY_ZERO, "Invalid division by zero occured")
ERROR_CODE(ERROR_READ_FAILED, "Reading from file failed")
ERROR_CODE(ERROR_WRITE_FAILED, "Writing to file failed")
ERROR_CODE(ERROR_OPEN_FAILED, "File opening failed")
ERROR_CODE(ERROR_NUM_ARGS, "Cluster expected 3 Arguments, but received a different number")
