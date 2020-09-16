#include <stdio.h>
#include <stdlib.h>
#include "error_codes.h"

static void print_and_exit(int code, const char * str) {
    fprintf(stderr, "%s\n",str);
    exit(code);
}

void panic(enum error_code error_code) {
    #define ERROR_CODE(_x_, _y_) do { \
        if ((_x_) == error_code) {    \
            print_and_exit(_x_, _y_); \
        }                             \
    } while(0);

    #include "xmacro_error_codes.h"
    #undef ERROR_CODE

    print_and_exit(ERROR_INVALID_CODE, "Unrecognized error code");
}
