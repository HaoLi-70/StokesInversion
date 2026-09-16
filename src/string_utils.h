
#pragma once

/*--------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdbool.h>
#include <stddef.h>

/*--------------------------------------------------------------------------------*/

#define Key_Length 64
#define Max_Line_Length 512

/*--------------------------------------------------------------------------------*/

extern ptrdiff_t STR_COUNT_CHAR(const char *str, const char c);

extern ptrdiff_t STR_INDEX_CHAR(const char *str, const char c, size_t order);

extern size_t STR_TRIM_LEFT(char *str);

extern size_t STR_TRIM_RIGHT(char *str);

extern size_t STR_TRIM(char *str);

extern void STR_COPY(char *dest, size_t dest_size, const char *source_rank, 
    size_t srcsize, bool trim_flag);

extern int STR_SPLIT(char *dest, size_t dest_size, char *source_rank, 
    size_t source_size);

extern void STR_TOUPPER(char *str);

extern ptrdiff_t STR_ELEMENTS(const char *str);

extern ptrdiff_t STR_READ_LINE(char **line, size_t *size, FILE *fa);

/*--------------------------------------------------------------------------------*/
