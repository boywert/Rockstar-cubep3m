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

void rescale_xv(float *xv, int np_local, int block, float a) {
  float H0 = 100.;   //[h*km]/[sec*Mpc]
  float vunit_compute = BOX_SIZE * 1.5 * sqrt(Om) * H0 / (2.*(float)CUBEP3M_NP)/a;  //!km/s
  float  lunit_compute = BOX_SIZE/(2.*(float)CUBEP3M_NP); //Mpc/h
  int i,j,k;
  float offset1,offset2,offset3;
  offset1 = offset2 = offset3 = 0.;

  for(k=0;k<CUBEP3M_NDIM;k++) 
    for(j=0;j<CUBEP3M_NDIM;j++) 
      for(i=0;i<CUBEP3M_NDIM;i++) {
	if(block == i+j*CUBEP3M_NDIM+k*CUBEP3M_NDIM*CUBEP3M_NDIM) {
	  offset1 = (float)i*BOX_SIZE/(float)CUBEP3M_NDIM;
	  offset2 = (float)j*BOX_SIZE/(float)CUBEP3M_NDIM;
	  offset3 = (float)k*BOX_SIZE/(float)CUBEP3M_NDIM;
	}
	else {
	  printf("cannot find the block number in range\n");
	  exit(1);
	}
      }
  for(i=0;i<np_local;i++) {
    xv[6*i] *= lunit_compute;
    xv[6*i] += offset1;
    xv[i*6+1] *= lunit_compute;
    xv[i*6+1] += offset2; 
    xv[i*6+2] *= lunit_compute; 
    xv[1*6+2] += offset3;
    xv[i*6+3] *= vunit_compute;
    xv[i*6+4] *= vunit_compute;
    xv[i*6+5] *= vunit_compute;
  }
}


void load_particles_cubep3m(char *filename, struct particle **p, int64_t *num_p) {
  FILE *input;
  char xvfile[1024], PIDfile[1024];
  char buffer[1024] = {'\0'};
  float *xv;
  int64_t i, *PID;
  int block;
  struct cubep3m_header {
    int np_local;
    float a,t,tau;
    int nts;
    float dt_f_acc,dt_pp_acc,dt_c_acc;
    int cur_checkpoint,cur_projection,cur_halofind;
    float mass_p;
  } header1,header2;

  block = string_replace_getblock(buffer,filename,"xvPID","xv");
  strcpy(xvfile,buffer);
  block = string_replace_getblock(buffer,filename,"xvPID","PID");
  strcpy(PIDfile,buffer);

  input = check_fopen(xvfile,"rb");
  fread(&header1, sizeof(struct cubep3m_header),1, input);

  *p = (struct particle *)check_realloc(*p, ((*num_p)+header1.np_local)*sizeof(struct particle), "Allocating particles.");

  xv = malloc(sizeof(float)*header1.np_local*6);
  fread(xv, sizeof(float),6*header1.np_local, input);
  fclose(input);

  TOTAL_PARTICLES = (int64_t)CUBEP3M_NP*(int64_t)CUBEP3M_NP*(int64_t)CUBEP3M_NP;
  SCALE_NOW = header1.a;
  PARTICLE_MASS = Om*CRITICAL_DENSITY * pow(BOX_SIZE, 3) / TOTAL_PARTICLES;

  rescale_xv(xv, header1.np_local, block, header1.a);  
  AVG_PARTICLE_SPACING = cbrt(PARTICLE_MASS / (Om*CRITICAL_DENSITY));
  for(i=0;i<header1.np_local;i++) {
    memcpy(&((*p)[(*num_p)+i].pos[0]),&(xv[i*6]),sizeof(float)*6);
  }
  free(xv);
  if(CUBEP3M_PID == 1) {
    input = check_fopen(PIDfile,"rb");
    fread(&header2, sizeof(struct cubep3m_header),1, input);
    if(header1.np_local != header2.np_local) {
      printf("np_local not consistent.");
      exit(1);
    }
    PID = malloc(sizeof(int64_t)*header1.np_local);
    fread(PID, sizeof(int64_t),header1.np_local, input);
    fclose(input);

    for(i=0;i<header1.np_local;i++) {
      (*p)[(*num_p)+i].id = PID[i];
    }
  
    free(PID);
  }
  else if(CUBEP3M_PID == 0) {
    for(i=0;i<header1.np_local;i++) {
      (*p)[(*num_p)+i].id = (int64_t)block*CUBEP3M_NP*CUBEP3M_NP*CUBEP3M_NP/(CUBEP3M_NDIM*CUBEP3M_NDIM*CUBEP3M_NDIM) + i ;
    }
  }
  *num_p += header1.np_local;
}
