
#pragma once

/*--------------------------------------------------------------------------------*/

#include <stdbool.h>
#include <mpi.h>

/*--------------------------------------------------------------------------------*/

#define MAX_MESSAGE_LENGTH 2000
#define MAX_BUFFER_SIZE 2500

extern char message_buffer[MAX_MESSAGE_LENGTH];

typedef enum ERR_LVL {ERR_LVL_WARNING, ERR_LVL_ERROR} ERR_LVL;

/*--------------------------------------------------------------------------------*/

extern _Noreturn void ABORTED(void);

extern int LOG_INIT(const char *filename);

extern void LOG_SET_CONSOLE(bool enabled);

extern int LOG_FINALIZE(void);

extern int LOG_WRITE(const char *msg, bool to_screen, bool verbose_flag);

extern void LOG_ERROR(ERR_LVL lvl, const char *routine, const char *msg);

extern bool FILE_EXIST(const char *filename);

/*--------------------------------------------------------------------------------*/
