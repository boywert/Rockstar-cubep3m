#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <string.h>
#include "io_util.h"
#include "../check_syscalls.h"
#include "../particle.h"
#include "../config_vars.h"
#include "../universal_constants.h"
#include "../config.h"

struct cubep3m_header {
  int np_local;
  float a,t,tau;
  int nts;
  float dt_f_acc,dt_pp_acc,dt_c_acc;
  int cur_checkpoint,cur_projection,cur_halofind;
  float mass_p;
};


int string_replace_getblock(char *out, char *in, char *find, char *replace) {
  char *p,*q;
  char buf[1000] = {"\0"};
  int block;
  int i,j=0;
  int lenin = strlen(in);
  int lenfind = strlen(find);
  int lenreplace = strlen(replace);
  if((p=strstr(in,find)) != NULL) {
    int first_replace = p-&in[0];
    for(i=0;i<first_replace;i++) {
      out[i] = in[i];
      j++;
    }
    for(i=0;i<lenreplace;i++) {
      out[i+first_replace]=replace[i]; j++;
    }
    for(i=0;i<lenin-first_replace-lenfind;i++) {
      out[first_replace+lenreplace+i] = in[first_replace+lenfind+i];j++;
    }
    // get block by checking the number before .dat
    q = strstr(in,".dat");
    int last_block_pos = q-&in[0];
    int first_block_pos = first_replace+lenfind;
    int blocklen = last_block_pos-first_block_pos;
    strncpy(buf,&in[first_block_pos],blocklen);
    sscanf(buf,"%d",&block);
    return block;
  }
  else {
    sprintf(out,"%s",in);
    return -1;
  }
}

void cubep3m_read_PID(FILE *fp, int block, int np_local, struct particle **p, int64_t *num_p) {
  int i;
  // start after header
  for(i=0;i<np_local;i++) 
    fread(&((*p)[(*num_p)+i].id), sizeof(int64_t),1, fp);  

}
void cubep3m_read_xv(FILE *fp, int block, int np_local, float a, struct particle **p, int64_t *num_p) {
  float H0 = 100.;   //[h*km]/[sec*Mpc]
  float vunit_compute = BOX_SIZE * 1.5 * sqrt(Om) * H0 / (2.*(float)CUBEP3M_NP)/a;  //!km/s
  float  lunit_compute = BOX_SIZE/(2.*(float)CUBEP3M_NP); //Mpc/h
  int i,j,k;
  float offset[3] = {0.};
  float buf[6];

  for(k=0;k<CUBEP3M_NDIM;k++) 
    for(j=0;j<CUBEP3M_NDIM;j++) 
      for(i=0;i<CUBEP3M_NDIM;i++) {
	if(block == i+j*CUBEP3M_NDIM+k*CUBEP3M_NDIM*CUBEP3M_NDIM) {
	  offset[0] = (float)i*BOX_SIZE/(float)CUBEP3M_NDIM;
	  offset[1] = (float)j*BOX_SIZE/(float)CUBEP3M_NDIM;
	  offset[2] = (float)k*BOX_SIZE/(float)CUBEP3M_NDIM;
	}
      }
  // start after header
  for(i=0;i<np_local;i++) {
    fread(buf, sizeof(float),6, fp);
    for(j=0;j<3;j++) {
      buf[j] *= lunit_compute;
      buf[j] += offset[j];
    }
    for(j=3;j<6;j++)
      buf[j] *= vunit_compute;
    memcpy(&((*p)[(*num_p)+i].pos[0]),buf,sizeof(float)*6);    
  }
}


void load_particles_cubep3m(char *filename, struct particle **p, int64_t *num_p) {
  FILE *input;
  char xvfile[1024], PIDfile[1024];
  char buffer[1024] = {'\0'};
  int block;
  struct cubep3m_header header1,header2;

  block = string_replace_getblock(buffer,filename,"xvPID","xv");
  strcpy(xvfile,buffer);
  block = string_replace_getblock(buffer,filename,"xvPID","PID");
  strcpy(PIDfile,buffer);

  input = check_fopen(xvfile,"rb");
  fread(&header1, sizeof(struct cubep3m_header),1, input);
  *p = (struct particle *)check_realloc(*p, ((*num_p)+header1.np_local)*sizeof(struct particle), "Allocating particles.");

  cubep3m_read_xv(input, block, header1.np_local, header1.a, p, num_p);
  fclose(input);
  
  TOTAL_PARTICLES = (int64_t)CUBEP3M_NP*(int64_t)CUBEP3M_NP*(int64_t)CUBEP3M_NP;
  SCALE_NOW = header1.a;
  PARTICLE_MASS = Om*CRITICAL_DENSITY * pow(BOX_SIZE, 3) / TOTAL_PARTICLES;
  AVG_PARTICLE_SPACING = cbrt(PARTICLE_MASS / (Om*CRITICAL_DENSITY));


  if(CUBEP3M_PID == 1) {
    input = check_fopen(PIDfile,"rb");
    fread(&header2, sizeof(struct cubep3m_header),1, input);
    if(header1.np_local != header2.np_local) {
      printf("np_local not consistent.");
      exit(1);
    }
    cubep3m_read_PID(input, block, header1.np_local, p, num_p);
    fclose(input);
  }
  else if(CUBEP3M_PID == 0)
    IGNORE_PARTICLE_IDS = 1;
  else {
    fprintf(stderr, "[Error] Unrecognized CUBEP3M_PID = %" PRId64 "!\n", CUBEP3M_PID);
    exit(1);
  }
  *num_p += header1.np_local;
}
