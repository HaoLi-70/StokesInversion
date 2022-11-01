
#include "ERROR.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
     revision log:
        8 Sept. 2021.
     
     ######################################################################*/

/*--------------------------------------------------------------------------------*/

extern void Error(enum error_level error_lv, const char *routine_name,\
                  const char *message_str){
    
    /*######################################################################
      Purpose:
        error handler.
      Record of revisions:
        29 Nov. 2019.
      Input parameters:
        error_lv, type of the error.
        routine_name, the routine name. 
        message_str, the message to pring.
      Return:
        .
     ######################################################################*/

    switch(error_lv){
      
      case enum_error:
        fprintf(stderr, "\n-WARNING in routine %s\n %s\n", routine_name, \
          (message_str) ? message_str : " (Undocumented)\n");

           
        return;
      case enum_warning:
        fprintf(stderr, "\n-WARNING in routine %s\n %s\n", routine_name, \
          (message_str) ? message_str : " (Undocumented)\n");
      
        return;
      default:
              /*
           
               */
      return;
    }
    return;
}

/*--------------------------------------------------------------------------------*/

