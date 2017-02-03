#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>


int main(int argc, char* argv[]) { 


	char* posIn;
	int divisionsX,divisionsY,divisionsZ,mode;
	float box[3];
	int indicator = 0;
	
	
	posIn = argv[1];
	sscanf(argv[2],"%d",&divisionsX); 
	sscanf(argv[3],"%d",&divisionsY);
	sscanf(argv[4],"%d",&divisionsZ);
	sscanf(argv[5],"%f",&box[0]);
	sscanf(argv[6],"%f",&box[1]);			//Box Dimension automatisch gewinnen ?!
	sscanf(argv[7],"%f",&box[2]);
	sscanf(argv[8],"%d",&mode);
	
	
	int npart,i,j,z;
	float posX,posY,posZ;
	int frames = 0;
	int error = 0;
	
	
	switch(mode) {
	
		case 1:  {
	
			float* densityX = (float*) calloc(divisionsX,sizeof(float));
			float* densityY = (float*) calloc(divisionsY,sizeof(float));
			float* densityZ = (float*) calloc(divisionsZ,sizeof(float));
	
			FILE* xin = fopen(posIn,"r");
	
			fscanf(xin,"%d\n",&npart); //read number of particles
			fscanf(xin, "%*[^\n]\n", NULL);
	
			error = fscanf(xin,"atom %f %f %f \n",&posX,&posY,&posZ);
	
			while( error != EOF) { //Solange es noch etwas zu lesen gibt? 
				for(i=0; i<npart;i++) {
					posX += box[0]/2;
					posY += box[1]/2;
					posZ += box[2]/2;
					densityX[(int)(posX / (box[0] / divisionsX))] += 1 / (box[1] * box[2] * (box[0] / divisionsX ) ); 
					densityY[(int)(posY / (box[1] / divisionsY))] += 1 / (box[0] * box[2] * (box[1] / divisionsX ) ); 
					densityZ[(int)(posZ / (box[2] / divisionsZ))] += 1 / (box[0] * box[1] * (box[2] / divisionsX ) ); 
					error = fscanf(xin,"atom %f %f %f \n",&posX,&posY,&posZ);
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
		fprintf(densZ,"Cell densities in Z direction: Cell length %f \n", box[2] / divisionsZ );
		
		for(i=0; i < divisionsZ; i++) {
		fprintf(densZ, "%d %8.8f \n",i, densityZ[i] / (float)(frames));
		}
	}
		}
		break;
		
		case 2: {

		
				float cellDensities[divisionsX][divisionsY][divisionsZ];
				float average,var,varVar,averageTemp;
		
		
				for (i=0;i< divisionsX ;i++) {
					for (j = 0; j < divisionsY; j++) {
						for (z = 0; z < divisionsZ; z++) {
							cellDensities[i][j][z] = 0;			//Initialisiere cellDensities Array
						}
					}
				}
		
		
				FILE* xin = fopen(posIn,"r");
	
				fscanf(xin,"%d\n",&npart); //read number of particles
				fscanf(xin, "%*[^\n]\n", NULL);
		
				error = fscanf(xin,"atom %f %f %f \n",&posX,&posY,&posZ);
			
		
				while( error != EOF) { //Solange es noch etwas zu lesen gibt? 
		
					for(i=0; i<npart;i++) {
						posX += box[0]/2;
						posY += box[1]/2;
						posZ += box[2]/2;
						cellDensities[(int)(posX / (box[0] / divisionsX))][(int)(posY / (box[1] / divisionsY))][(int)(posZ / (box[2] / divisionsZ))] +=  1 / ( (box[0] / divisionsX) * (box[1] / divisionsY) * (box[2] / divisionsZ) ) ;
						error = fscanf(xin,"atom %f %f %f \n",&posX,&posY,&posZ);
					}
					frames++; //Hiermit dann mitzählen wie viele Frames abgelaufen - Am Ende die Dichte damit verrechnen und ausgeben
		
		
		
					//Berechnung roh_bar / m^2 / m^4 
			
			
					//average entspricht roh_bar
					for (i=0;i< divisionsX ;i++) {
						for (j = 0; j < divisionsY; j++) {
							for (z = 0; z < divisionsZ; z++) {
								averageTemp += cellDensities[i][j][z];			
							}
						}
					}
			
					averageTemp /= divisionsX * divisionsY * divisionsZ;
					average += averageTemp;
		
		
					//var entspricht m^2 
					for (i=0;i< divisionsX ;i++) {
						for (j = 0; j < divisionsY; j++) {
							for (z = 0; z < divisionsZ; z++) {
								var += (averageTemp - cellDensities[i][j][z]) * (averageTemp - cellDensities[i][j][z]) ;			
							}
						}
					}
		
		
					//varVar entspricht m^4
					for (i=0;i< divisionsX ;i++) {
						for (j = 0; j < divisionsY; j++) {
							for (z = 0; z < divisionsZ; z++) {
								varVar += (averageTemp - cellDensities[i][j][z]) * (averageTemp - cellDensities[i][j][z]) * (averageTemp - cellDensities[i][j][z])* (averageTemp - cellDensities[i][j][z]) ;			
							}
						}
					}
			
					for (i=0;i< divisionsX ;i++) {
						for (j = 0; j < divisionsY; j++) {
							for (z = 0; z < divisionsZ; z++) {
								cellDensities[i][j][z] = 0;			//Reset Cell Densities
							}
						}
					}
		
				}
		
				average /= frames;
				var /= (frames*divisionsX*divisionsY*divisionsZ);
				varVar /= (frames*divisionsX*divisionsY*divisionsZ);
		
		
		
				FILE* binder = fopen("./BinderParameter.xyz","w");
		
				if(!binder) {
					printf("File not found!");
					exit(-1);
				}
		
				fprintf(binder,"%d\n",npart);
				fprintf(binder,"Moments of density: Average density %f | Variance of average density | Variance of the variance \n",(float)(npart)/(box[0]*box[1]*box[2]) ) ;
				
				fprintf(binder, "%8.8f %8.3f %8.3f \n", average, var , varVar );
				fprintf(binder, "Binder parameter: %f",varVar/(var * var));
				
			}
		
		break;
	
		case 3: {
			
			FILE* xin = fopen(posIn,"r");
			
			fscanf(xin,"%d\n",&npart); //read number of particles
			fscanf(xin, "%*[^\n]\n", NULL);
			
			
			float* densityZ = (float*) calloc(divisionsZ,sizeof(float));
			
			FILE* file = fopen("./HistogramDensities.xyz","w");
		
			if(!file) {
				printf("File not found!");
				exit(-1);
			}
				
			fprintf(file,"%d\n",npart);
		
			error = fscanf(xin,"atom %f %f %f \n",&posX,&posY,&posZ);
			
			while( error != EOF) { //Solange es noch etwas zu lesen gibt? 
		
				for(i=0; i<npart;i++) {
					posZ += box[2]/2;
					densityZ[(int)(posZ / (box[2] / divisionsZ))] += 1 / (box[0] * box[1] * (box[2] / divisionsX ) ); 
					error = fscanf(xin,"atom %f %f %f \n",&posX,&posY,&posZ);
				}
				
				frames++; //Hiermit dann mitzählen wie viele Frames abgelaufen - Am Ende die Dichte damit verrechnen und ausgeben
				
				
				
				for (i=0;i<divisionsZ;i++) {
					fprintf(file,"%f \n", densityZ[i] );
					densityZ[i]=0;
				}
				
				
				
			}
			
			
		
			
	
	
		}
		break;
	}
	
	}
	
