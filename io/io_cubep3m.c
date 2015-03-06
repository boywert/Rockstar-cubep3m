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

void string_replace(char *out, char *in, char *find, char *replace) {
  char *p;
  int i,j;
  int lenin = strlen(in);
  int lenfind = strlen(find);
  int lenreplace = strlen(replace);
  // strcpy(out,"\0");
  j = 0;
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
      printf("i=%d,c=%c\n",i,in[first_replace+lenfind+i]);
      out[first_replace+lenreplace+i] = in[first_replace+lenfind+i];j++;
    }
    //strcat(out,"\0");
  }
  else
    sprintf(out,"%s",in);
  printf("expect = %d, real = %d, count = %d\n",lenin+lenreplace-lenfind,strlen(out),j);
}
void rescale_xv(float *xv, int np_local) {
  float H0 = 100.;   //[h*km]/[sec*Mpc]
  float RHO_CRIT_0 = 2.7755397e11;   //[h^2*Msun]/[Mpc^3] 
}

void load_particles_cubep3m(char *filename, int block, struct particle **p, int64_t *num_p) {
  FILE *input;
  char xvfile[1024], PIDfile[1024];
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
  printf("filename = %s\n",filename);
  string_replace(buffer,filename,"xvPID","xv");
  strcpy(xvfile,buffer);
  string_replace(buffer,filename,"xvPID","PID");
  strcpy(PIDfile,buffer);
  printf("xv = %s, pid = %s\n",xvfile,PIDfile);
  exit(0);
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

  /* input = check_fopen(PIDfile,"rb"); */
  /* fread(&header2, sizeof(struct cubep3m_header),1, input); */
  /* if(header1.np_local != header2.np_local) { */
  /*   printf("np_local not consistent."); */
  /*   exit(1); */
  /* } */
  /* PID = malloc(sizeof(int64_t)*header1.np_local); */
  /* fread(PID, sizeof(int64_t),header1.np_local, input); */
  /* fclose(input); */
  /* for(i=0;i<header1.np_local;i++) { */
  /*   (*p)[(*num_p)+i].id = PID[i]; */
  /* } */
  
  /* free(PID); */
  *num_p += header1.np_local;
}
