
#include "string_utils.h"
#include <ctype.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Module:
        Bounds-checked string parsing, trimming, case conversion, token
        searches, and dynamically sized line input.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/


ptrdiff_t STR_COUNT_CHAR(const char *str, const char c){

    /*######################################################################
      Purpose:
        Count occurrences of a character in a string.
      Input parameters:
        str, the input string.
        c, the character.
      Return:
        Number of occurrences, or -1 for a NULL string.
    ######################################################################*/

    if(!str) return -1;
    
    ptrdiff_t count = 0;
    const char *p = str;

    while(*p){
      if(*p == c) count++;
      p++;
    }
    
    return count;
}

/*--------------------------------------------------------------------------------*/

ptrdiff_t STR_INDEX_CHAR(const char *str, const char c, size_t order){
    
    /*######################################################################
      Purpose:
        Find the requested occurrence of a character in a string.
      Input parameters:
        str, the input string.
        c, the character.
        order, the order of the character.
      Return:
        Zero-based index, or a negative status code when not found or 
        invalid.
    ######################################################################*/

    if(!str || *str=='\0') return -1;
    if(order == 0) return -2;
    
    size_t count = 0;
    const char *p = str;
    while(*p){
      if(*p == c && ++count == order)
        return p - str;
      p++;
    }

    return -3;
}

/*--------------------------------------------------------------------------------*/


static ptrdiff_t STR_INDEX_SPACE(const char *str, size_t order){
    
    /*######################################################################
      Purpose:
        Find the requested whitespace character in a string.
      Input parameters:
        str, the input string.
        order, the order of the space.
      Return:
        Zero-based index, or a negative status code when not found or 
        invalid.
    ######################################################################*/

    if(str == NULL || !*str) return -1;
    if(order == 0) return -2;

    size_t count = 0;
    const char *p = str;
    while(*p){
      if(isspace((unsigned char)*p) && ++count == order)
        return p - str;
      p++;
    }

    return -3;
}

/*--------------------------------------------------------------------------------*/

size_t STR_TRIM_LEFT(char *str){
    
    /*######################################################################
      Purpose:
        Remove leading whitespace from a string.
      Input parameters:
        str, the input string.
      Output parameters:
        str, the output string.
      Return:
        return the number of removed spaces.
    ######################################################################*/

    if(!str || !*str) return 0;
    
    char *p = str;
    while(*p && isspace((unsigned char)*p)){
      p++;
    }

    size_t count = (size_t)(p-str);
    if(count>0){
      memmove(str, p, strlen(p)+1);
    }

    return count;
}

/*--------------------------------------------------------------------------------*/

size_t STR_TRIM_RIGHT(char *str){
    
    /*######################################################################
      Purpose:
        Remove trailing whitespace from a string.
      Input parameters:
        str, the input string.
      Output parameters:
        str, the output string.
      Return:
        return the number of removed spaces.
    ######################################################################*/

    if(!str || *str=='\0') return 0;

    size_t length = strlen(str);
    size_t end = length;
    while(end > 0 && isspace((unsigned char)str[end-1])) end--;
    str[end] = '\0';

    return length-end;
}

/*--------------------------------------------------------------------------------*/

size_t STR_TRIM(char *str){

    /*######################################################################
      Purpose:
        Remove leading and trailing whitespace from a string.
      Input parameters:
        str, the input string.
      Output parameters:
        str, the output string.
      Return:
        return the number of removed spaces.
    ######################################################################*/

    if (!str || !*str) return 0;

    size_t length = strlen(str);
    size_t begin = 0;
    while(begin < length && isspace((unsigned char)str[begin])) begin++;
    size_t end = length;
    while(end > begin && isspace((unsigned char)str[end-1])) end--;

    size_t kept = end-begin;
    if(begin > 0 && kept > 0) memmove(str, str+begin, kept);
    str[kept] = '\0';

    return length-kept;
}

/*--------------------------------------------------------------------------------*/

void STR_COPY(char *dest, size_t dest_size, const char *source_rank, 
    size_t srcsize, bool trim_flag){
    
    /*######################################################################
      Purpose:
        Copy a bounded string into a destination buffer.
      Input parameters:
        dest_size, buffer size
        source_rank, the input string.
        srcsize, the source_rank size.
        trim_flag, if the flag > 0, remove the leading and trailing spaces.
      Output parameters:
        dest, the output string.
    ######################################################################*/

    if(!dest || dest_size == 0) return;
    if(!source_rank || srcsize == 0){
      dest[0] = '\0';
      return;
    }

    size_t size = dest_size>srcsize? srcsize : dest_size-1;

    memmove(dest, source_rank, size);

    dest[size] = '\0';
    
    if(trim_flag) STR_TRIM(dest);
    
    return;
}

/*--------------------------------------------------------------------------------*/

int STR_SPLIT(char *dest, size_t dest_size, char *source_rank, 
    size_t source_size){
    
    /*######################################################################
      Purpose:
        Extract the first whitespace-delimited element from a string.
      Input parameters:
        source_rank, the input string.
        dest_size, buffer size.
        source_size, capacity of the mutable source buffer.
      Output parameters:
        dest, the copied element.
        source_rank, the left elements.
      Return:
        1 when elements remain, 0 after the last element, or -1 on invalid 
        input.
    ######################################################################*/

    if(!dest || dest_size == 0 || !source_rank || source_size == 0) return -1;

    STR_TRIM(source_rank);
    size_t total_length = strlen(source_rank);

    ptrdiff_t space_idx = STR_INDEX_SPACE(source_rank, 1);
    if(space_idx < 0){  
      STR_COPY(dest, dest_size, source_rank, total_length, true);
      source_rank[0] = '\0';
      return 0;

    }else{
      size_t split = (size_t)space_idx;
      STR_COPY(dest, dest_size, source_rank, split, true);
      size_t begin = split;
      while(begin < total_length 
          && isspace((unsigned char)source_rank[begin])) begin++;
      size_t remaining = total_length-begin;
      if(remaining >= source_size) return -1;
      memmove(source_rank, source_rank+begin, remaining);
      source_rank[remaining] = '\0';
    }

    return 1;
}

/*--------------------------------------------------------------------------------*/

void STR_TOUPPER(char *str){
    
    /*######################################################################
      Purpose:
        Convert all characters in a string to uppercase.
      Input parameters:
        str, the input string.
      Output parameters:
        str, the output string.
    ######################################################################*/

    if(!str) return; 

    for(char *p=str; *p; p++){
      *p = (char)toupper((unsigned char)*p);
    }

    return;
}

/*--------------------------------------------------------------------------------*/

ptrdiff_t STR_ELEMENTS(const char *str){

    /*######################################################################
      Purpose:
        Count whitespace-delimited elements in a string.
      Input parameters:
        str, the input string.
      Return:
        the number of the elements.
    ######################################################################*/

    if(!str) return -1; 

    ptrdiff_t count = 0;
    bool in_element = false;
    for(const char *p=str; *p; p++){
      if(isspace((unsigned char)*p)){
        in_element = false;
      }else if(!in_element){
        count++;
        in_element = true;
      }
    }
    return count;
}

/*--------------------------------------------------------------------------------*/

ptrdiff_t STR_READ_LINE(char **line, size_t *size, FILE *fa){

    /*######################################################################
      Purpose:
        Read one logical line, growing the destination buffer as needed.
      Input parameters:
        fa, the file.
      Output parameters:
        line, the output string.
        size, buffer size.
      Return:
        Line length; -1 on EOF, -2 on allocation or length overflow,
          -3 for invalid input, or -4 on a file read error.
    ######################################################################*/

    if(!line || !size || !fa) return -3;

    if(*line != NULL && *size == 0) return -3;
    if(*line == NULL){
      *size = 512;
      *line = malloc(*size);
      if (!*line) return -2;
    }

    while(1){

      size_t len = 0;
      size_t nspaces;

      while(1){
        size_t available = *size-len;
        int chunk_capacity = available > (size_t)INT_MAX 
            ? INT_MAX : (int)available;
        if(!fgets(*line + len, chunk_capacity, fa)){
          if(len == 0) return ferror(fa) ? -4 : -1;
          break;
        }

        size_t chunk = strlen(*line+len);
        len += chunk;
        if (len > 0 && (*line)[len-1] == '\n') break;

        if(*size > SIZE_MAX/2) return -2;
        size_t new_size = (*size)*2;
        char *tmp = realloc(*line, new_size);
        if (!tmp) return -2;
        *line = tmp;
        *size = new_size;
      }

      if(len == 0) continue;

      nspaces = STR_TRIM_LEFT(*line);
      len -= nspaces;

      char first = (*line)[0];
      if(first=='#' || first=='!' || first=='*' || first=='\0'){
        continue;
      }

      char *comment = strpbrk(*line, "#!");
      if(comment){
        len = (size_t)(comment-*line);
        *comment = '\0';
      }
      nspaces = STR_TRIM_RIGHT(*line);
      len -= nspaces;

      if(len > (size_t)PTRDIFF_MAX) return -2;
      
      return (ptrdiff_t)len;
    }
}

/*--------------------------------------------------------------------------------*/
