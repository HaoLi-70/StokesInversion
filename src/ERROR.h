
#ifndef ERROR_h
#define ERROR_h

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <ctype.h>
#include <string.h>

/*--------------------------------------------------------------------------------*/

enum error_level {enum_warning, enum_error};

/*--------------------------------------------------------------------------------*/

#define MAX_MESSAGE_LENGTH 2000

/*--------------------------------------------------------------------------------*/

extern void Error(enum error_level error_lv, const char *routine_name,
                  const char *message_str);

/*--------------------------------------------------------------------------------*/

#endif /* ERROR_h */
