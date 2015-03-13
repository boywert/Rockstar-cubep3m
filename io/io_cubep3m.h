#ifndef _IO_CUBEP3M_H
#define _IO_CUBEP3M_H

#include <stdint.h>
#include "../particle.h"

void load_particles_cubep3m_zip2015(char *filename, struct particle **p, int64_t *num_p);
void load_particles_cubep3m(char *filename, struct particle **p, int64_t *num_p);

#endif /* _IO_CUBEP3M_H */
