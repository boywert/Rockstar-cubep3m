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

void get_xvfile(char *buffer, char *string, int maxlen) {
  int64_t i=0, out=0, l=strlen(string);
  assert(snap < NUM_SNAPS);
  snprintf(buffer, maxlen, "");
  out=strlen(buffer);
  for (; (i<l)&&(out < (maxlen-1)); i++) {
    if (string[i] != '<') { buffer[out]=string[i]; buffer[out+1]=0; }
    else {
      if (!strncmp(string+i, "<xvPID>", 7)) {
	i+=6;
	snprintf(buffer+out, maxlen-out,"xv",);
      }
      else buffer[out] = string[i];
    }
    out = strlen(buffer);
  }
  buffer[out] = 0;
}
void get_PIDfile(char *buffer, char *string, int maxlen) {
  int64_t i=0, out=0, l=strlen(string);
  assert(snap < NUM_SNAPS);
  snprintf(buffer, maxlen, "");
  out=strlen(buffer);
  for (; (i<l)&&(out < (maxlen-1)); i++) {
    if (string[i] != '<') { buffer[out]=string[i]; buffer[out+1]=0; }
    else {
      if (!strncmp(string+i, "<xvPID>", 7)) {
	i+=6;
	snprintf(buffer+out, maxlen-out,"xv",);
      }
      else buffer[out] = string[i];
    }
    out = strlen(buffer);
  }
  buffer[out] = 0;
}
void rescale_xv(float *xv, int np_local) {
  float H0 = 100.;   //[h*km]/[sec*Mpc]
  float RHO_CRIT_0 = 2.7755397e11;   //[h^2*Msun]/[Mpc^3] 
}

void load_particles_cubep3m(char *filename, int block, struct particle **p, int64_t *num_p) {
  FILE *input;
  char xvfile[], PIDfile[1024];
  char buffer[1024];
  float *xv;
  int64_t i,n, *PID;
  struct cubep3m_header {
    int np_local;
    float a,t,tau;
    int nts;
    float dt_f_acc,dt_pp_acc,dt_c_acc;
    int cur_checkpoint,cur_projection,cur_halofind;
    float mass_p;
  } header1,header2;
  get_xvfile(buffer,filename,1024);
  sprintf(xvfile,"%s",buffer);
  get_PIDfile(buffer,filename,1024);
  sprintf(PIDfile,"%s",buffer);
  
  printf("xv = %s, pid = %s\n",xvfile,PIDfile);
  input = check_fopen(xvfile,"rb");
  fread(&header1, sizeof(struct cubep3m_header),1, input);
  
  *p = (struct particle *)check_realloc(*p, ((*num_p)+header1.np_local)*sizeof(struct particle), "Allocating particles.");

  xv = malloc(sizeof(float)*header1.np_local*6);
  fread(xv, sizeof(float),6*header1.np_local, input);
  fclose(input);
  
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
