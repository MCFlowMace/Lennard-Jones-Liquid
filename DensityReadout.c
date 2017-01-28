#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>


int main(int argc, char* argv[]) { 


	char* posIn;
	int divisionsX,divisionsY,divisionsZ;
	float box[3];
	int indicator = 0;
	
	
	posIn = argv[1];
	sscanf(argv[2],"%d",&divisionsX); 
	sscanf(argv[3],"%d",&divisionsY);
	sscanf(argv[4],"%d",&divisionsZ);
	sscanf(argv[5],"%f",&box[0]);
	sscanf(argv[6],"%f",&box[1]);			//Box Dimension automatisch gewinnen ?!
	sscanf(argv[7],"%f",&box[2]);
	
	int npart,i;
	float posX,posY,posZ;
	int frames = 0;
	int error = 0;
	
	
	float* densityX = (float*) calloc(divisionsX,sizeof(float));
	float* densityY = (float*) calloc(divisionsY,sizeof(float));
	float* densityZ = (float*) calloc(divisionsZ,sizeof(float));
	
	
	FILE* xin = fopen(posIn,"r");
	
	fscanf(xin,"%d\n",&npart); //read number of particles
	fscanf(xin, "%*[^\n]\n", NULL);
	
	error = fscanf(xin,"atom %f %f %f \n",&posX,&posY,&posZ);
	
	while( error!= EOF) { //Solange es noch etwas zu lesen gibt? 
		for(i=0; i<npart;i++) {
			error = fscanf(xin,"atom %f %f %f \n",&posX,&posY,&posZ);
			densityX[(int)(posX / (box[0] / divisionsX))] += 1; 
			densityY[(int)(posY / (box[1] / divisionsY))] += 1; 
			densityZ[(int)(posZ / (box[2] / divisionsZ))] += 1; 
		}
		frames++; //Hiermit dann mitzählen wie viele Frames abgelaufen - Am Ende die Dichte damit verrechnen und ausgeben
	}
	
	if (divisionsX != 0) {
		FILE* densX = fopen("./densitiesX.xyz","w");
		
		if(!densX) {
		printf("File not found!");
		exit(-1);
		}
	
		fprintf(densX,"%d\n",npart);
		fprintf(densX,"Cell densities in X direction: Cell length %f \n", box[0] / divisionsX );
		
		for(i=0; i < divisionsX; i++) {
		fprintf(densX, "%d %8.8f \n",i, densityX[i] / (float)(frames));
		}
	}
	
	
	if (divisionsY != 0) {
		FILE* densY = fopen("./densitiesY.xyz","w");
		
		if(!densY) {
		printf("File not found!");
		exit(-1);
		}
	
		fprintf(densY,"%d\n",npart);
		fprintf(densY,"Cell densities in Y direction: Cell length %f \n", box[1] / divisionsY );
		
		for(i=0; i < divisionsY; i++) {
		fprintf(densY, "%d %8.8f \n",i, densityY[i] / (float)(frames));
		}
	}
	
	if (divisionsZ != 0) {
		FILE* densZ = fopen("./densitiesZ.xyz","w");
		
		if(!densZ) {
		printf("File not found!");
		exit(-1);
		}
	
		fprintf(densZ,"%d\n",npart);
		fprintf(densZ,"Cell densities in X direction: Cell length %f \n", box[2] / divisionsZ );
		
		for(i=0; i < divisionsZ; i++) {
		fprintf(densZ, "%d %8.8f \n",i, densityZ[i] / (float)(frames));
		}
	}
	
}