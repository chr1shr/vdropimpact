#include <cstdlib>

#include "common.hh"

/** \brief Opens a file and checks the operation was successful.
 *
 * Opens a file, and checks the return value to ensure that the operation
 * was successful.
 * \param[in] filename the file to open.
 * \param[in] mode the cstdio fopen mode to use.
 * \return The file handle. */
FILE* safe_fopen(const char* filename,const char* mode) {
    FILE *temp=fopen(filename,mode);
    if(temp==NULL) {
        fprintf(stderr,"Error opening file \"%s\"\n",filename);
        exit(1);
    }
    return temp;
}

/** \brief Reads data from a file and checks that the operation was successful.
 *
 * Reads data from a file and checks the return value to ensure that the
 * operation was successful. If not successful, it prints an error message and
 * exits.
 * \param[in] ptr the memory to write to.
 * \param[in] size the size of each element to be read.
 * \param[in] count the number of elements to be read.
 * \param[in] fp the file handle to read from.
 * \param[in] p a description of what is being read, used in the error message
 *              if the operation isn't successful. */
void safe_fread(void *ptr,size_t size,size_t count,FILE *fp,const char* p) {
    if(fread(ptr,size,count,fp)!=count) {
        fprintf(stderr,"Can't read %s from file\n",p);
        exit(1);
    }
}

/** \brief Function for printing fatal error messages and exiting.
 *
 * Function for printing fatal error messages and exiting.
 * \param[in] p a pointer to the message to print.
 * \param[in] status the status code to return with. */
void fatal_error(const char *p,int code) {
    fprintf(stderr,"Error: %s\n",p);
    exit(code);
}
