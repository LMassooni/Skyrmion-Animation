#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#define FRANDOM (rand()/(RAND_MAX+1.0))
#define PI M_PI
#define L 70
#define N (L*L)
#define D 1.
#define J 1.
#define B 0
#define H 0.4
#define p 500000

int right[N], left[N], up[N], down[N];
const int prog_bar_length = 30;

void update(double progress,double total) {
	int num = progress*prog_bar_length/total;
	printf("\r[");
	for (int i=0;i<num;i++){
		printf("#");
	}
	for(int i=0;i<prog_bar_length- num;i++){
		printf(" ");
	}
printf("]%lf%% Done",(progress/(total-1))*100);
fflush(stdout);
}

//calcula os primeiros vizinhos da melhor forma que eu consegui pensar (condições periodicas de contorno)
void vizinhos(){
    for (int i = 0; i<N; i++){
        right[i] = i+1;
        left[i] = i-1;
        up[i] = i-L;
        down[i] = i+L;
        if(i<L){  //primeira linha
            up[i] = i+L*L-L;
        }
        if(i>=L*(L-1)){ //ultima linha
        down[i] = i-(L*L-L);
        }
        if(i%L==0){ //primeira coluna
        left[i] = i+L-1;
        }
        if(i%L==L-1){//ultima coluna
        right[i] = i-L+1;
        }
}
}


void initGnuplot(){
  printf("set term wxt size %d,%d\n", 1000, 1000);
  printf("set size -1\n");
  printf("set autoscale fix\n");
  printf("unset xtics\n");
  printf("unset ytics\n");
  printf("unset key\n");
  printf("unset colorbox\n");

  printf("set palette defined ( -1 \"dark-spring-green\", 1 \"dark-goldenrod\" )\n");
  
}

void printGnuplot(double s[N][2]){
  int ix, iy;
  double sz,sx,sy;
  printf("plot '-' using 1:2:3 with image, '-' using 1:2:4:5 w vectors filled lw 2 lc \"white\"\n");
  for (int a = 0; a < N ;a++) {
      ix = a%L;
      iy = a/L;
      sz = cos(s[a][1]);
printf("%d %d %lf\n", ix, iy, sz);
  }
  printf("e\n");
  for(int j=0;j<N;j++){
  ix = j%L;
  iy = j/L;
  sz = cos(s[j][1]);
  sx = sin(s[j][1])*cos(s[j][0]);
  sy = sin(s[j][1])*sin(s[j][0]);
  printf("%d %d %lf %lf %lf\n", ix, iy, sz, sx,sy);
  }
 printf("e\n");
 }

void metropolis_otimizado(double s[N][2],double temp, double *ene){
double dH,dHx,dHy,dHz;
double newphi,newtheta;
int sitio;
for(int a  =0; a<N; a++){
sitio = FRANDOM*N;
newphi = s[sitio][0]+(1.-2.*FRANDOM)*1.;
newtheta = acos(1.-2.*FRANDOM);
 dHx = -J*(sin(newtheta)*cos(newphi) - sin(s[sitio][1])*cos(s[sitio][0])) * (
    		 sin(s[up[sitio]][1]) * cos(s[up[sitio]][0]) +
     		 sin(s[down[sitio]][1]) * cos(s[down[sitio]][0]) +
    		 sin(s[left[sitio]][1]) * cos(s[left[sitio]][0]) +
    		 sin(s[right[sitio]][1]) * cos(s[right[sitio]][0])
 )
 -D*((sin(newtheta)*sin(newphi)-sin(s[sitio][1])*sin(s[sitio][0]))*
 (cos(s[right[sitio]][1])-cos(s[left[sitio]][1]))-
 (cos(newtheta)-cos(s[sitio][1]))*
 (sin(s[right[sitio]][1])*sin(s[right[sitio]][0])
 -sin(s[left[sitio]][1])*sin(s[left[sitio]][0])));
 
dHy = -J*(sin(newtheta)*sin(newphi) - sin(s[sitio][1])*sin(s[sitio][0])) * (
		    sin(s[up[sitio]][1]) * sin(s[up[sitio]][0]) +
		    sin(s[down[sitio]][1]) * sin(s[down[sitio]][0]) +
		    sin(s[left[sitio]][1]) * sin(s[left[sitio]][0]) +
		    sin(s[right[sitio]][1]) * sin(s[right[sitio]][0])
		)
  -D*((cos(newtheta)-cos(s[sitio][1]))*
 (sin(s[down[sitio]][1])*cos(s[down[sitio]][0])
 -sin(s[up[sitio]][1])*cos(s[up[sitio]][0]))-
  (sin(newtheta)*cos(newphi)-sin(s[sitio][1])*cos(s[sitio][0]))*
 (cos(s[down[sitio]][1])-cos(s[up[sitio]][1])));

 dHz = -J*(cos(newtheta)-cos(s[sitio][1]))*(
  			cos(s[up[sitio]][1])+cos(s[down[sitio]][1])+
 			cos(s[left[sitio]][1])+cos(s[right[sitio]][1]))
 -H*(cos(newtheta)-cos(s[sitio][1]))
 +B*(pow(cos(newtheta),2)-pow(cos(s[sitio][1]),2));
 
 dH = dHx+dHy+dHz;
      if (dH <=0 || FRANDOM < exp (-dH / temp)){
          s[sitio][0] = newphi;
          s[sitio][1] = newtheta;
          *ene += dH;

   }
   }
    }

    




void energia_inicial(double s[N][2],double *ene){
  int sitio;
  *ene=0;
  for (sitio = 0; sitio < N; sitio++){
  *ene += -J/2.*(sin(s[sitio][1])
  			*(sin(s[up[sitio]][1])*cos(s[sitio][0]-s[up[sitio]][0])
  			+sin(s[down[sitio]][1])*cos(s[sitio][0]-s[down[sitio]][0])
  			+sin(s[left[sitio]][1])*cos(s[sitio][0]-s[left[sitio]][0])
  			+sin(s[right[sitio]][1])*cos(s[sitio][0]-s[right[sitio]][0]))+
 	 	   cos(s[sitio][1])
 	 	   *(cos(s[up[sitio]][1])
 	 	   +cos(s[down[sitio]][1])
 	 	   +cos(s[left[sitio]][1])
 	 	   +cos(s[right[sitio]][1])))
 -D*(sin(s[sitio][1])*sin(s[sitio][0])*cos(s[right[sitio]][1])
 -cos(s[sitio][1])*sin(s[right[sitio]][1])*sin(s[right[sitio]][0])
 +cos(s[sitio][1])*sin(s[down[sitio]][1])*cos(s[down[sitio]][0])
 -sin(s[sitio][1])*cos(s[sitio][0])*cos(s[down[sitio]][1]))
  -H*cos(s[sitio][1])+B*cos(s[sitio][1])*cos(s[sitio][1]);
}
}

	void overrelax(double s[N][2]){
	double heffx, heffy, heffz,prod;
	double nsx,nsz,nsy;
	int i, sit;
	for(i=0;i<N;i++){
	sit = FRANDOM*N;
	heffx = J*(sin(s[up[sit]][1])*cos(s[up[sit]][0])
		+sin(s[down[sit]][1])*cos(s[down[sit]][0])
		+sin(s[left[sit]][1])*cos(s[left[sit]][0])
		+sin(s[right[sit]][1])*cos(s[right[sit]][0]))
	-D*(-cos(s[up[sit]][1])+cos(s[down[sit]][1]));
		
	heffy = J*(sin(s[up[sit]][1])*sin(s[up[sit]][0])
		+sin(s[down[sit]][1])*sin(s[down[sit]][0])
		+sin(s[left[sit]][1])*sin(s[left[sit]][0])
		+sin(s[right[sit]][1])*sin(s[right[sit]][0]))
	+D*(-cos(s[left[sit]][1])
		+cos(s[right[sit]][1]));
		
	heffz = J*(cos(s[up[sit]][1])
		+cos(s[left[sit]][1])
		+cos(s[down[sit]][1])
		+cos(s[right[sit]][1]))
	+D*(-cos(s[up[sit]][0])*sin(s[up[sit]][1])
	+cos(s[down[sit]][0])*sin(s[down[sit]][1])
	+sin(s[left[sit]][0])*sin(s[left[sit]][1])
	-sin(s[right[sit]][0])*sin(s[right[sit]][1]))+H-B;
		
	prod = 2.*(cos(s[sit][0])*sin(s[sit][1])*heffx
	+sin(s[sit][0])*sin(s[sit][1])*heffy
	+cos(s[sit][1])*heffz)
	/(heffx*heffx+heffy*heffy+heffz*heffz);

	nsz = -cos(s[sit][1])+prod*heffz;
	nsx = -cos(s[sit][0])*sin(s[sit][1])+prod*heffx;
	nsy = -sin(s[sit][0])*sin(s[sit][1])+prod*heffy;

	s[sit][0] = atan2(nsy,nsx);
	s[sit][1] = acos(nsz);	
}
}

void grandezas(double s[N][2], double *magx, double *magy, double *magz){
  int sitio;
  *magx=0;
  *magy=0;
  *magz=0;
  for (sitio = 0; sitio < N; sitio++){
    *magx += sin(s[sitio][1])*sin(s[sitio][0]);
    *magy += sin(s[sitio][1])*cos(s[sitio][0]);
    *magz += cos(s[sitio][1]);
    
}
*magx=*magx/N;
*magy=*magy/N;	
*magz=*magz/N;
}

int main(){
double s[N][2];
    double T = 0.1*J;
    int i;
    double ene;
    srand(time(NULL));
    vizinhos();
initGnuplot();
//for(int k=0;k<n;k++){
//update(k,n);

for (i =0; i<N;i++){
    s[i][0]=FRANDOM*2.*PI;
    s[i][1]=FRANDOM*PI;
}

energia_inicial(s, &ene);

for(int j = 0; j<p;j++){
            metropolis_otimizado(s, T, &ene);
           for(int x = 0; x<5; x++)
            	overrelax(s);
            if(j%100==0)
            	printGnuplot(s);

}

    return 0;
}