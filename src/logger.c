
#include "logger.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Module:
        Process-local logging, checked file lifecycle management, and safe
        fatal-error termination for the MPI runtime.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

static FILE *log_file = NULL;
static bool log_write_error_reported = false;
static bool console_enabled = true;
static char buffer[MAX_BUFFER_SIZE];
char message_buffer[MAX_MESSAGE_LENGTH];

/*--------------------------------------------------------------------------------*/

_Noreturn void ABORTED(void){

    /*######################################################################
      Purpose:
        Terminate the MPI job after a fatal error.
      Input parameters:
        None.
    ######################################################################*/

    fflush(NULL);

    int initialized = 0;
    int finalized = 0;
    if(MPI_Initialized(&initialized) == MPI_SUCCESS && initialized
        && MPI_Finalized(&finalized) == MPI_SUCCESS && !finalized){
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    exit(EXIT_FAILURE);
}

/*--------------------------------------------------------------------------------*/

int LOG_INIT(const char *filename){

    /*######################################################################
      Purpose:
        Open the log file for appending.
      Input parameters:
        filename, path of the log file; NULL or an empty string disables it.
      Return:
        0 on success; -1 when the new file cannot be opened or the previous
          file cannot be closed.
    ######################################################################*/

    if(!filename || !filename[0]) return LOG_FINALIZE();

    FILE *new_file = fopen(filename, "a");
    if(!new_file) return -1;

    int status = LOG_FINALIZE();
    log_file = new_file;
    log_write_error_reported = false;

    return status == 0 ? 0 : -1;
}

/*--------------------------------------------------------------------------------*/

void LOG_SET_CONSOLE(bool enabled){

    /*######################################################################
      Purpose:
        Enable or suppress terminal output for the current process.
      Input parameters:
        enabled, true when LOG_WRITE may write to stderr.
      Note:
        File logging is unaffected. MPI programs use this to reserve terminal
        output for world rank 0 while retaining per-island and per-rank logs.
    ######################################################################*/

    console_enabled = enabled;
}

/*--------------------------------------------------------------------------------*/

int LOG_FINALIZE(void){

    /*######################################################################
      Purpose:
        Close the active log file.
      Input parameters:
        None.
      Return:
        0 on success; EOF if flushing or closing the file fails.
    ######################################################################*/

    if(!log_file) return 0;

    FILE *file = log_file;
    log_file = NULL;
    int flush_status = fflush(file);
    int close_status = fclose(file);

    return flush_status == 0 && close_status == 0 ? 0 : EOF;
}

/*--------------------------------------------------------------------------------*/

static int LOG_WRITE_LINE(FILE *file, const char *msg, size_t length){

    if(fwrite(msg, 1, length, file) != length || fputc('\n', file) == EOF
        || fflush(file) != 0) return -1;

    return 0;
}

/*--------------------------------------------------------------------------------*/

int LOG_WRITE(const char *msg, bool to_screen, bool verbose_flag){

    /*######################################################################
      Purpose:
        Write a message to the requested logging destinations.
      Input parameters:
        msg, the log message.
        to_screen, whether to also write the message to stderr.
        verbose_flag, whether this message is enabled.
      Return:
        0 on success or when disabled; -1 for an invalid message or a write
          failure.
    ######################################################################*/

    if(!verbose_flag) return 0;
    if(!msg) return -1;

    size_t length = strlen(msg);
    while(length > 0 && (msg[length-1] == '\n' || msg[length-1] == '\r')){
      length--;
    }

    int status = 0;
    if(to_screen && console_enabled
        && LOG_WRITE_LINE(stderr, msg, length) != 0) status = -1;
    if(log_file && LOG_WRITE_LINE(log_file, msg, length) != 0){
      status = -1;
      if(!log_write_error_reported){
        static const char warning[] =
            "-WARNING in routine LOG_WRITE: Failed to write the log file.";
        if(console_enabled){
          LOG_WRITE_LINE(stderr, warning, sizeof(warning)-1);
        }
        log_write_error_reported = true;
      }
    }

    return status;
}

/*--------------------------------------------------------------------------------*/

void LOG_ERROR(ERR_LVL lvl, const char *rname, const char *msg){

    /*######################################################################
      Purpose:
        Format and report a warning or fatal error.
      Input parameters:
        lvl, error severity.
        rname, function name.
        msg, the log message.
    ######################################################################*/

    const char *safe_rname = rname ? rname : "unknown";
    const char *safe_msg = msg ? msg : "no error message";
    int length = snprintf(buffer, sizeof(buffer), "%s in routine %s: %s",
        lvl == ERR_LVL_ERROR ? "-ERROR" : "-WARNING", safe_rname, safe_msg);
    if(length < 0){
      memcpy(buffer, "-ERROR: failed to format log message",
          sizeof("-ERROR: failed to format log message"));
    }else if((size_t)length >= sizeof(buffer)){
      static const char suffix[] = "... [truncated]";
      memcpy(buffer+sizeof(buffer)-sizeof(suffix), suffix, sizeof(suffix));
    }

    LOG_WRITE(buffer, true, true);

    if(lvl==ERR_LVL_ERROR) ABORTED();

    return;
}

bool FILE_EXIST(const char *filename){

    /*######################################################################
      Purpose:
        Check whether a file can be opened for reading.
      Input parameters:
        filename, path of the file.
      Return:
        true when the file is readable; otherwise false.
    ######################################################################*/

    if(!filename || !filename[0]) return false;

    FILE *file = fopen(filename, "r");
    
    if(file){
      fclose(file);
      return true;
    }
    return false;
}

/*--------------------------------------------------------------------------------*/
