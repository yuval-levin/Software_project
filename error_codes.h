#ifndef _ERROR_CODES_H_
#define _ERROR_CODES_H_

enum error_code {
    ERROR_SUCCESS = 0, /* so that everything that follows will get a non 0 value. */
	/*macro: search ERROR_CODE(x,y) and replace with x,*/
#define ERROR_CODE(x, y) x, 

#include "xmacro_error_codes.h"
#undef ERROR_CODE /*ERROR CODE no longer defined*/

    ERROR_INVALID_CODE /* so we'll have a non-valid error code. */
};

/* Prints the associated error message and exit the program. */
void panic(enum error_code error_code);

#endif

