#ifndef _IO_CUBEP3M_H
#define _IO_CUBEP3M_H

#include <stdint.h>
#include "../particle.h"

//void gzip_file(char *filename);
void load_particles_cubep3m(char *prefix, struct particle **p, int64_t *num_p);

#endif /* _IO_CUBEP3M_H */
