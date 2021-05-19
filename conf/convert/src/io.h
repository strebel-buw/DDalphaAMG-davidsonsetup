#define RELOAD_SERIAL_COMP  0
#define RELOAD_SERIAL 1
#define RELOAD_MG4QCD 2
#define RELOAD_DDHMC 3

#define SAVE_SERIAL_COMP 0
#define SAVE_SERIAL 1
#define SAVE_MG4QCD 2
#define SAVE_DDHMC 3

/* Helps in defining globals */
#ifdef CONTROL
#define EXTERN
#else
#define EXTERN extern
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h> /* For "write" and "close" "off_t" */
#endif
#include <sys/types.h> /* For "off_t" */
#include <stdio.h>


void save_lattice( int flag, char *filename );
void reload_lattice( int flag, char *filename );
void dump_conf( void );