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

#define CUBEP3M_BUFFER_SIZE 100000

int FORCE_BYTESWAP = 0;

struct cubep3m_header {
  int np_local;
  float a,t,tau;
  int nts;
  float dt_f_acc,dt_pp_acc,dt_c_acc;
  int cur_checkpoint,cur_projection,cur_halofind;
  float mass_p;
};
struct cubep3m_header_extend {
  float v_r2;
  float var2; // dummy
  float var3; // dummy
  float var4; // dummy
};

void byte_swap(void* input, int len) {
  unsigned char a;
  unsigned char* e = (unsigned char *)input;
  int half = len/2;
  int i;
  if(len%2 == 1) {
    fprintf(stderr,"Trying to swap odd-byte thing.\nExit\n");
    exit(1);
  }
#define SWAP(x,y) a=e[x]; e[x]=e[y]; e[y]=a;
  for(i=0;i<half;i++) {
    SWAP(i,len-i-1);
  }
#undef SWAP
}

static inline void swap_cubep3m_header(struct cubep3m_header* h) {
  byte_swap(&(h->np_local),4);
  byte_swap(&(h->a),4);
  byte_swap(&(h->t),4);
  byte_swap(&(h->tau),4);
  byte_swap(&(h->nts),4);
  byte_swap(&(h->dt_f_acc),4);
  byte_swap(&(h->dt_pp_acc),4);
  byte_swap(&(h->dt_c_acc),4);
  byte_swap(&(h->cur_checkpoint),4);
  byte_swap(&(h->cur_projection),4);
  byte_swap(&(h->cur_halofind),4);
  byte_swap(&(h->mass_p),4);
}
static inline void swap_cubep3m_header_extend(struct cubep3m_header_extend* h) {
  byte_swap(&(h->v_r2),4);
}

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
    out[j] = 0;
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

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

void cubep3m_read_PID(FILE *fp, int block, int np_local, struct particle **p, int64_t *num_p) {
  int64_t i,j;
  int64_t read_n;
  int64_t count = 0;
  int64_t buffer[CUBEP3M_BUFFER_SIZE];
  // start after header
  for (j=0;j<=np_local/CUBEP3M_BUFFER_SIZE;j++) {
    read_n = MIN(np_local-j*CUBEP3M_BUFFER_SIZE,CUBEP3M_BUFFER_SIZE);
    fread(buffer, sizeof(int64_t),read_n, fp);
    for(i=0;i<read_n;i++) {
      if(FORCE_BYTESWAP)
	byte_swap(&(buffer[i]),8);
      memcpy(&((*p)[(*num_p)+j*CUBEP3M_BUFFER_SIZE+i].id), &(buffer[i]),sizeof(int64_t));
      count++;
    }
  }
  if(count != np_local) {
    fprintf(stderr,"particle not consistent\nExit\n");
    exit(1);
  }
}

void cubep3m_read_zip2015(FILE *fp0, FILE *fp1, FILE *fp2, FILE *fp3, int block, int np_local, float a, struct particle **p, int64_t *num_p) {
  float H0 = 100.;   //[h*km]/[sec*Mpc]
  float vunit_compute = BOX_SIZE * 1.5 * sqrt(Om) * H0 / (2.*(float)CUBEP3M_NP)/a;  //!km/s
  float  lunit_compute = BOX_SIZE/(2.*(float)CUBEP3M_NP); //Mpc/h
  int i,j,k,n;
  float offset[3] = {0.};
  uint8_t rhoc_i1,xi1[3];
  int32_t rhoc_i4;
  int16_t vi2[3];
  int cur_p = 0;
  float mesh_scale = 2.*CUBEP3M_NP/CUBEP3M_NDIM/CUBEP3M_NCDIM;
  struct cubep3m_header_extend header1,header2;

  for(k=0;k<CUBEP3M_NDIM;k++)
    for(j=0;j<CUBEP3M_NDIM;j++)
      for(i=0;i<CUBEP3M_NDIM;i++) {
	if(block == i+j*CUBEP3M_NDIM+k*CUBEP3M_NDIM*CUBEP3M_NDIM) {
	  offset[0] = (float)i*BOX_SIZE/(float)CUBEP3M_NDIM;
	  offset[1] = (float)j*BOX_SIZE/(float)CUBEP3M_NDIM;
	  offset[2] = (float)k*BOX_SIZE/(float)CUBEP3M_NDIM;
	}
      }

  fread(&header1, sizeof(struct cubep3m_header_extend),1, fp0);
  fread(&header2, sizeof(struct cubep3m_header_extend),1, fp1);
  if(FORCE_BYTESWAP) {
    swap_cubep3m_header_extend(&header1);
    swap_cubep3m_header_extend(&header2);
  }
  if(header1.v_r2 != header2.v_r2) {
    printf("v_r2 not consistent.\n");
    exit(1);
  }
  for(k=0;k<CUBEP3M_NCDIM;k++)
    for(j=0;j<CUBEP3M_NCDIM;j++)
      for(i=0;i<CUBEP3M_NCDIM;i++) {
	fread(&rhoc_i1,sizeof(uint8_t),1,fp2);
	rhoc_i4 = (int32_t)rhoc_i1;
	if(rhoc_i1 == 255) {
	  fread(&rhoc_i4,sizeof(int32_t),1,fp3);
	  if(FORCE_BYTESWAP)
	    byte_swap(&rhoc_i4,4);
	}
	for(n=0;n<rhoc_i4;n++) {
	  fread(xi1,sizeof(uint8_t),3,fp0);
	  fread(vi2,sizeof(int16_t),3,fp1);
	  if(FORCE_BYTESWAP) {
	    byte_swap(&(vi2[0]),2);
	    byte_swap(&(vi2[1]),2);
	    byte_swap(&(vi2[2]),2);
	  }
	  (*p)[(*num_p)+cur_p].pos[0] = mesh_scale*(((float)xi1[0]+0.5)/256.+i)*lunit_compute + offset[0];
	  (*p)[(*num_p)+cur_p].pos[1] = mesh_scale*(((float)xi1[1]+0.5)/256.+j)*lunit_compute + offset[1];
	  (*p)[(*num_p)+cur_p].pos[2] = mesh_scale*(((float)xi1[2]+0.5)/256.+k)*lunit_compute + offset[2];
	  (*p)[(*num_p)+cur_p].pos[3] = (float)vi2[0]/header1.v_r2*vunit_compute;
	  (*p)[(*num_p)+cur_p].pos[4] = (float)vi2[1]/header1.v_r2*vunit_compute;
	  (*p)[(*num_p)+cur_p].pos[5] = (float)vi2[2]/header1.v_r2*vunit_compute;
	  cur_p++;
	}
      }
  if(cur_p != np_local) {
    printf("Counted particles are not consistent.\n");
    exit(1);
  }
  /* FILE* testfp; */
  /* testfp = check_fopen("test.xv","wb+"); */
  /* for(i=0;i<np_local;i++) */
  /*   fwrite(&((*p)[(*num_p)+i].pos[0]),sizeof(float),6,testfp); */
  /* fclose(testfp); */
}

void cubep3m_read_xv(FILE *fp, int block, int np_local, float a, struct particle **p, int64_t *num_p) {
  float H0 = 100.;   //[h*km]/[sec*Mpc]
  float vunit_compute = BOX_SIZE * 1.5 * sqrt(Om) * H0 / (2.*(float)CUBEP3M_NP)/a;  //!km/s
  float  lunit_compute = BOX_SIZE/(2.*(float)CUBEP3M_NP); //Mpc/h
  int i,j,k;
  int read_n;
  float offset[3] = {0.};
  float buffer[6*CUBEP3M_BUFFER_SIZE];

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
  for (k=0;k<=np_local/CUBEP3M_BUFFER_SIZE;k++) {
    read_n = MIN(np_local-k*CUBEP3M_BUFFER_SIZE,CUBEP3M_BUFFER_SIZE);
    fread(buffer, sizeof(float),read_n*6, fp);
    for(i=0;i<read_n;i++) {
      for(j=0;j<3;j++) {
	if(FORCE_BYTESWAP)
	  byte_swap(&(buffer[6*i+j]),4);
	buffer[6*i+j] *= lunit_compute;
	buffer[6*i+j] += offset[j];
      }
      for(j=3;j<6;j++) {
	if(FORCE_BYTESWAP)
	  byte_swap(&(buffer[6*i+j]),4);
	buffer[6*i+j] *= vunit_compute;
      }
      memcpy(&((*p)[(*num_p)+k*CUBEP3M_BUFFER_SIZE+i].pos[0]),&(buffer[6*i]),sizeof(float)*6);
    }
  }
}
void load_particles_cubep3m_zip2015(char *filename, struct particle **p, int64_t *num_p) {
  FILE *zip_fp[4],*input;
  char zip0_file[1024],zip1_file[1024],zip2_file[1024],zip3_file[1024], PIDfile[1024];
  char buffer[1024] = {'\0'};
  int i,block;
  struct cubep3m_header header1,header2;
  block = string_replace_getblock(buffer,filename,"xvPID","zip0_");
  strcpy(zip0_file,buffer);
  block = string_replace_getblock(buffer,filename,"xvPID","zip1_");
  strcpy(zip1_file,buffer);
  block = string_replace_getblock(buffer,filename,"xvPID","zip2_");
  strcpy(zip2_file,buffer);
  block = string_replace_getblock(buffer,filename,"xvPID","zip3_");
  strcpy(zip3_file,buffer);
  block = string_replace_getblock(buffer,filename,"xvPID","PID");
  strcpy(PIDfile,buffer);
  printf("file = %s\n",PIDfile);
  zip_fp[0] = check_fopen(zip0_file,"rb");
  zip_fp[1] = check_fopen(zip1_file,"rb");
  zip_fp[2] = check_fopen(zip2_file,"rb");
  zip_fp[3] = check_fopen(zip3_file,"rb");

  fread(&header1, sizeof(struct cubep3m_header),1, zip_fp[0]);
  fread(&header2, sizeof(struct cubep3m_header),1, zip_fp[1]);
  if(header1.mass_p != 8.)
    FORCE_BYTESWAP = 1;
  if(FORCE_BYTESWAP) {
    swap_cubep3m_header(&header1);
    swap_cubep3m_header(&header2);
    if(header1.mass_p != 8.) {
      printf("mass_p = %g line:%d\nExit\n",header1.mass_p,__LINE__ );
      exit(1);
    }
  }
  if(header1.np_local != header2.np_local) {
    printf("np_local not consistent.");
    exit(1);
  }
  *p = (struct particle *)check_realloc(*p, ((*num_p)+header1.np_local)*sizeof(struct particle), "Allocating particles.");

  cubep3m_read_zip2015(zip_fp[0],zip_fp[1],zip_fp[2],zip_fp[3], block, header1.np_local, header1.a, p, num_p);
  for(i=0;i<4;i++)
    fclose(zip_fp[i]);

  TOTAL_PARTICLES = (int64_t)CUBEP3M_NP*(int64_t)CUBEP3M_NP*(int64_t)CUBEP3M_NP;
  SCALE_NOW = header1.a;
  PARTICLE_MASS = Om*CRITICAL_DENSITY * pow(BOX_SIZE, 3) / TOTAL_PARTICLES;
  AVG_PARTICLE_SPACING = cbrt(PARTICLE_MASS / (Om*CRITICAL_DENSITY));

  if(CUBEP3M_PID == 1) {
    input = check_fopen(PIDfile,"rb");
    fread(&header2, sizeof(struct cubep3m_header),1, input);
    if(FORCE_BYTESWAP)
      swap_cubep3m_header(&header2);
    if(header1.np_local != header2.np_local) {
      printf("np_local not consistent.\n");
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
  if(header1.mass_p < 7.999 || header1.mass_p > 8.001) {
    FORCE_BYTESWAP = 1;
  }
  if(FORCE_BYTESWAP) {
    swap_cubep3m_header(&header1);
    if(header1.mass_p < 7.999 || header1.mass_p > 8.001)) {
      printf("mass_p = %g\nExit\n",header1.mass_p);
      exit(1);
    }
  }
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
    if(FORCE_BYTESWAP)
      swap_cubep3m_header(&header2);
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
