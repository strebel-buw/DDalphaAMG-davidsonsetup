#include "includes.h"

int nx, ny, nz, nt;
int sites_on_node;
su3_matrix *tlink;


int main ( int argc, char *argv[])
{

    printf ("%s(%i) -> %s(%i)\n", argv[2], atoi(argv[1]), argv[4], atoi(argv[3]));
    fflush (0);
    reload_lattice (atoi(argv[1]), argv[2]);
    save_lattice (atoi(argv[3]), argv[4]);
    
    //dump_conf();

    return 0;
}
