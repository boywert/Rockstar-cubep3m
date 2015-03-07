#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <string.h>
#include "../check_syscalls.h"
#include "stringparse.h"
#include "../particle.h"
#include "../config_vars.h"

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
  float RHO_CRIT_0 = 2.7755397e11;   //[h^2*Msun]/[Mpc^3]
  float Omega0 = 0.3;
  float vunit_compute = BOX_SIZE * 1.5 * sqrt(Omega0) * H0 / (2.*(float)CUBEP3M_NP)/a;  //!km/s
  float munit_compute = BOX_SIZE * BOX_SIZE * BOX_SIZE * Omega0 * RHO_CRIT_0 / ((float)CUBEP3M_NP*(float)CUBEP3M_NP*(float)CUBEP3M_NP); // ! Msun/h
  float  lunit_compute = BOX_SIZE/(2.*(float)CUBEP3M_NP); //Mpc/h
  int i,j,k;
  float offset1,offset2,offset3;
  PARTICLE_MASS = munit_compute;
  for(k=0;k<CUBEP3M_NDIM;k++) 
    for(j=0;j<CUBEP3M_NDIM;j++) 
      for(i=0;i<CUBEP3M_NDIM;i++) 
	if(block == i+j*CUBEP3M_NDIM+k*CUBEP3M_NDIM*CUBEP3M_NDIM) {
	  offset1 = (float)i*BOX_SIZE/(float)CUBEP3M_NDIM;
	  offset2 = (float)j*BOX_SIZE/(float)CUBEP3M_NDIM;
	  offset3 = (float)k*BOX_SIZE/(float)CUBEP3M_NDIM;
	}
  printf("offset1 = %f, offset2 = %f, offset3 = %f\n",offset1,offset2,offset3);
  for(i=0;i<np_local;i++) {
    xv[i] *= lunit_compute;
    xv[i] += offset1;
    xv[i+1] *= lunit_compute;
    xv[i+1] += offset2; 
    xv[i+2] *= lunit_compute; 
    xv[1+2] += offset3;
    xv[i+3] *= vunit_compute;
    xv[i+4] *= vunit_compute;
    xv[i+5] *= vunit_compute;
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
  printf("filename = %s\n",filename);
  block = string_replace_getblock(buffer,filename,"xvPID","xv");
  strcpy(xvfile,buffer);
  block = string_replace_getblock(buffer,filename,"xvPID","PID");
  strcpy(PIDfile,buffer);
  printf("block = %d, xv = %s, pid = %s\n",block,xvfile,PIDfile);
  input = check_fopen(xvfile,"rb");
  fread(&header1, sizeof(struct cubep3m_header),1, input);
  
  *p = (struct particle *)check_realloc(*p, ((*num_p)+header1.np_local)*sizeof(struct particle), "Allocating particles.");

  

  xv = malloc(sizeof(float)*header1.np_local*6);
  fread(xv, sizeof(float),6*header1.np_local, input);
  fclose(input);

  rescale_xv(xv, header1.np_local, block, header1.a);  
  for(i=0;i<100;i++)
    printf("%d, %f %f %f %f %f %f\n",i,xv[i],xv[i+1],xv[i+2],xv[i+3],xv[i+4],xv[i+5]);
  exit(1);
  for(i=0;i<header1.np_local;i++) {
    memcpy(&((*p)[(*num_p)+i].pos[0]),&(xv[i*6]),sizeof(float)*6);
  }
  free(xv);

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
  *num_p += header1.np_local;
}
