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

char *str_replace(char *string, char *sub, char *replace)
{
  printf("string = %s,sub=%s,replace=%s\n",string,sub,replace);
  if(!string || !sub || !replace) return NULL;
  char *pos = string; int found = 0;
  while((pos = strstr(pos, sub))){
    pos += strlen(sub);
    found++;
  }
  if(found == 0) return NULL;
  int size = ((strlen(string) - (strlen(sub) * found)) + (strlen(replace) * found)) + 1;
  char *result = (char*)malloc(size);
  pos = string;
  char *pos1;
  while((pos1 = strstr(pos, sub))){
    int len = (pos1 - pos);
    strncat(result, pos, len);
    strncat(result, replace, strlen(replace));
    pos = (pos1 + strlen(sub));
  }
  if(pos != (string + strlen(string)))
    strncat(result, pos, (string - pos));
  return result;
}
void rescale_xv(float *xv, int np_local) {
  float H0 = 100.;   //[h*km]/[sec*Mpc]
  float RHO_CRIT_0 = 2.7755397e11;   //[h^2*Msun]/[Mpc^3] 
}

void load_particles_cubep3m(char *filename, int block, struct particle **p, int64_t *num_p) {
  FILE *input;
  char *xvfile, *PIDfile;
  
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
  xvfile = str_replace(filename,"<cubep3m_magic>","xv");
  PIDfile = str_replace(filename,"<cubep3m_magic>","PID"); 

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
