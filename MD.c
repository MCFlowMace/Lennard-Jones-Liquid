#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>


void initStart(const int npart, const float temp, float** positions_ptr, float** ppositions_ptr, float** velocities_ptr, float** f_ptr, const float dt, const float *box) {
	
	int i,j,k,counter;
	
	float* sumv = (float*)calloc(3,sizeof(float));
	float sumv2 = 0;
	int dimensions = (int)ceil(cbrtf(npart));
	float dx = box[0]/dimensions;
	float dy = box[1]/dimensions;
	float dz = box[2]/dimensions;
	
	*positions_ptr = (float*) calloc(3*npart,sizeof(float));
	*velocities_ptr = (float*) calloc(3*npart,sizeof(float));
	*ppositions_ptr = (float*) calloc(3*npart,sizeof(float));
	*f_ptr = (float*) calloc(3*npart,sizeof(float));
	
	float* positions = *positions_ptr;
	float* velocities = *velocities_ptr;
	float* ppositions = *ppositions_ptr;
	
	fprintf(stderr,"dimensions: %d dx: %f\n",dimensions,dx);
	
	//srand(time(NULL));
	srand(42);
	
	for (i = 0; i < npart; i++) {
		for(j = 0; j < 3; j++){
			
		//positions[i*3 + j] = (rand()/(float)RAND_MAX)*box;
		velocities[i*3 + j] = (rand()/(float)RAND_MAX - 0.5);
		sumv[j] += velocities[i*3+j];
		sumv2 += velocities[i*3+j]*velocities[i*3+j];
		
		}
		
	}
	
	counter=0;
	for (i=0; i<dimensions; i++)
		for(j=0; j<dimensions; j++)
			for(k=0; k<dimensions; k++) {
				
				if(counter == npart) {
					j=dimensions;
					i=dimensions;
					break;
				}
				
				positions[(dimensions*dimensions*i+dimensions*j+k)*3] = k*dx +dx/2 -box[0]/2;
				positions[(dimensions*dimensions*i+dimensions*j+k)*3+1] = j*dy +dy/2 -box[1]/2;
				positions[(dimensions*dimensions*i+dimensions*j+k)*3+2] = i*dz +dz/2 -box[2]/2;
				counter++;
				//fprintf(stderr,"counter: %d i= %d j=%d k=%d\n",counter,i,j,k);
			}
	
	sumv[0] /= npart;
	sumv[1] /= npart;
	sumv[2] /= npart;
	sumv2 /= npart;
	
	float vCMS2 = sumv[0]*sumv[0] + sumv[1]*sumv[1] + sumv[2]*sumv[2];
	
	fprintf(stdout,"Velocity Center of mass: %f frame 0\n",vCMS2);
	
	float fs = sqrt(3*temp/sumv2);
	
	for (i = 0; i < npart; i++) {
		for(j = 0; j < 3; j++) {
		velocities[i*3 + j] -= sumv[j];
		velocities[i*3 + j] *= fs;
		ppositions[i*3 + j] = positions[i*3 + j] - velocities[i*3 + j]*dt;
		}
	}
	
	free(sumv);
	
}

void initContinued(int *npart, float** positions_ptr, float** velocities_ptr, float** f_ptr, FILE* xin, FILE* vin) {
	
	int i;
	
	srand(420);
	
	fscanf(xin,"%d\n",npart); //read number of particles
	
	int npartCheck;
	
	fscanf(vin,"%d\n",&npartCheck);
	
	
	if(npartCheck != *npart) {
		fprintf(stdout, "Error: Input files incompatible!\n");
		exit(-1);
	}
	
	fprintf(stdout,"Particles: %d\n",*npart);
	
	*positions_ptr = (float*) calloc(3*(*npart),sizeof(float));
	*velocities_ptr = (float*) calloc(3*(*npart),sizeof(float));
	*f_ptr = (float*) calloc(3*(*npart),sizeof(float));
	
	float* positions = *positions_ptr;
	float* velocities = *velocities_ptr;
	
	//skip second line
	fscanf(xin, "%*[^\n]\n", NULL);
	fscanf(vin, "%*[^\n]\n", NULL);
	
	
	for(i=0; i<*npart; i++) {
		fscanf(xin, "atom %f %f %f\n",&positions[3*i],&positions[3*i+1],&positions[3*i+2]);
		fscanf(vin, "atom %f %f %f\n",&velocities[3*i], &velocities[3*i+1], &velocities[3*i+2]);
		//fprintf(stderr,"i: %d\n",i);
	}
	
	
}


float gauss(float sigma, float mean) {
	
	float r =2.0;
	float v1,v2;
	
	while( r >= 1) {
		v1=2*(rand()/(float)RAND_MAX)-1;
		v2=2*(rand()/(float)RAND_MAX)-1;
		r=v1*v1+v2*v2;
		
	}
	float l= v1*sqrt(-2*log(r)/r);
	
	return mean+sigma*l;
	
}

#ifdef CUDA


void errorHandler  (cudaError_t error, int line){
  if(error != cudaSuccess)
  {
    // print the CUDA error message and exit
    fprintf(stderr,"CUDA error: %s in line number %d\n", cudaGetErrorString(error),line);
    exit(-1);
  }
}

__global__ void particleForce(float* f, const float* positions, const int npart, const float boxX, const float boxY, const float boxZ, const float cutoff2, const int n, const float ecut, float* energies) {
	
	
	int threadId = blockIdx.x*blockDim.x + threadIdx.x;
	
	if(threadId >= npart)
		return;
	
	int j;
	
	float rx;
	float ry;
	float rz;
	float r2;
	float r2i;
	float r6i;
	float ff; //lennard-jones potential
	
	f[3*threadId] = 0;
	f[3*threadId+1] =0;
	f[3*threadId+2] =0;
	
	energies[threadId] =0;
			
	for(j=0; j<npart; j++) {
			
						if(j==threadId) continue;
						
						rx = positions[3*threadId] - positions[3*j];
						ry = positions[3*threadId+1] - positions[3*j+1];
						rz = positions[3*threadId+2] - positions[3*j+2];
						rx -= boxX* round(rx/boxX);
						ry -= boxY* round(ry/boxY);
						rz -= boxZ* round(rz/boxZ);
						
						r2 = rx*rx + ry*ry + rz*rz;
						
						if(r2 < 0.1) printf("Not good! i %d j %d distance %f frame %d\n",threadId,j,r2,n);
						if(r2 < cutoff2) {
							r2i = 1/r2;
							r6i = r2i*r2i*r2i;
							ff = 48*r2i*r6i*(r6i - 0.5);
							f[3*threadId] += ff*rx;
							f[3*threadId+1] += ff*ry;
							f[3*threadId+2] += ff*rz;
							
							if(j>threadId)
								energies[threadId] += 4*r6i*(r6i-1) - ecut;
						}
						
	}
}

void calcForce(float* en, float* f, const float* positions, const int npart, const float *box, const float cutoff2, const int n, const float ecut) {
	
	
	float *devPtr_forces;
	float *devPtr_positions;
	float *devPtr_energies;
	int i;
	
	float *energies = (float*)calloc(npart,sizeof(float));

	errorHandler(cudaMalloc((void**)&devPtr_forces, 3*npart* sizeof(float)),__LINE__);
	errorHandler(cudaMalloc((void**)&devPtr_positions, 3*npart* sizeof(float)),__LINE__);
	errorHandler(cudaMalloc((void**)&devPtr_energies, npart*sizeof(float)),__LINE__);
	
	//errorHandler(cudaMemcpy(devPtr_forces, f, 3*npart * sizeof(float), cudaMemcpyHostToDevice),__LINE__);
	errorHandler(cudaMemcpy(devPtr_positions, positions, 3*npart * sizeof(float), cudaMemcpyHostToDevice),__LINE__);
	
	int threadsPerBlock = 512;
	int blocks = npart/threadsPerBlock + 1;
	
	//fprintf(stderr,"starting GPU calc\n");
	
	particleForce<<<blocks, threadsPerBlock>>>(devPtr_forces, devPtr_positions, npart, box[0], box[1], box[2], cutoff2, n,ecut, devPtr_energies);
	errorHandler( cudaPeekAtLastError(),__LINE__);
	
	errorHandler(cudaMemcpy(f, devPtr_forces, 3*npart * sizeof(float), cudaMemcpyDeviceToHost),__LINE__);
	errorHandler(cudaMemcpy(energies,devPtr_energies, npart*sizeof(float), cudaMemcpyDeviceToHost),__LINE__);
	
	errorHandler(cudaFree(devPtr_forces),__LINE__);
	errorHandler(cudaFree(devPtr_positions),__LINE__);
	errorHandler(cudaFree(devPtr_energies),__LINE__);
	
	*en =0;
	for(i=0; i<npart; i++) {
		*en += energies[i];
	}
	
	free(energies);

	
}


/*__global__ void integrateParticle1(const int npart, const float* f, float* positions, float* velocities, const float dt) {
	
	
	
	int threadId = blockIdx.x*blockDim.x + threadIdx.x;
	
	if(threadId >= npart)
		return;
		
	int j;
		
	for( j=0; j<3; j++) {
		positions[3*threadId+j] += dt* velocities[3*threadId+j] +dt*dt* f[3*threadId+j]/2;
		velocities[3*threadId+j] += dt*f[3*threadId+j]/2;
	}
	
}

__global__ void integrateParticle2(const int npart, const float* f, float* velocities, const float dt) {
	
	
	
	int threadId = blockIdx.x*blockDim.x + threadIdx.x;
	
	if(threadId >= npart)
		return;
		
	int j;
		
	for( j=0; j<3; j++) {
		velocities[3*threadId+j] += dt*f[3*threadId+j]/2;
	}
	
}

void integrateVelVerlet(const int part, const float* f, const int npart, float* positions, float* velocities, const float dt, const float temp_target, const float nu) {
	
	int i,j;
	
	if(part ==1) {
		
		
		float *devPtr_forces;
		float *devPtr_positions;
		float *devPtr_velocities;
	
		errorHandler(cudaMalloc((void**)&devPtr_forces, 3*npart* sizeof(float)),__LINE__);
		errorHandler(cudaMalloc((void**)&devPtr_positions, 3*npart* sizeof(float)),__LINE__);
		errorHandler(cudaMalloc((void**)&devPtr_velocities, 3*npart*sizeof(float)),__LINE__);
		
		errorHandler(cudaMemcpy(devPtr_forces, f, 3*npart * sizeof(float), cudaMemcpyHostToDevice),__LINE__);
		errorHandler(cudaMemcpy(devPtr_positions, positions, 3*npart * sizeof(float), cudaMemcpyHostToDevice),__LINE__);
		errorHandler(cudaMemcpy(devPtr_velocities, velocities, 3*npart *sizeof(float), cudaMemcpyHostToDevice),__LINE__);
		
		int threadsPerBlock = 512;
		int blocks = npart/threadsPerBlock + 1;
		
		//fprintf(stderr,"starting GPU calc\n");
		
		integrateParticle1<<<blocks, threadsPerBlock>>>(npart, devPtr_forces, devPtr_positions, devPtr_velocities, dt);
		errorHandler( cudaPeekAtLastError(),__LINE__);
		
		errorHandler(cudaMemcpy(positions, devPtr_positions, 3*npart * sizeof(float), cudaMemcpyDeviceToHost),__LINE__);
		errorHandler(cudaMemcpy(velocities, devPtr_velocities, 3*npart * sizeof(float), cudaMemcpyDeviceToHost),__LINE__);
		
		errorHandler(cudaFree(devPtr_forces),__LINE__);
		errorHandler(cudaFree(devPtr_positions),__LINE__);
		errorHandler(cudaFree(devPtr_velocities),__LINE__);
		
				
	}
	
	if(part ==2) {
		
		float *devPtr_forces;
		float *devPtr_velocities;
	
		errorHandler(cudaMalloc((void**)&devPtr_forces, 3*npart* sizeof(float)),__LINE__);
		errorHandler(cudaMalloc((void**)&devPtr_velocities, 3*npart*sizeof(float)),__LINE__);
		
		errorHandler(cudaMemcpy(devPtr_forces, f, 3*npart * sizeof(float), cudaMemcpyHostToDevice),__LINE__);
		errorHandler(cudaMemcpy(devPtr_velocities, velocities, 3*npart *sizeof(float), cudaMemcpyHostToDevice),__LINE__);
		
		int threadsPerBlock = 512;
		int blocks = npart/threadsPerBlock + 1;
		
		//fprintf(stderr,"starting GPU calc\n");
		
		integrateParticle2<<<blocks, threadsPerBlock>>>(npart, devPtr_forces, devPtr_velocities, dt);
		errorHandler( cudaPeekAtLastError(),__LINE__);
		
		errorHandler(cudaMemcpy(velocities, devPtr_velocities, 3*npart * sizeof(float), cudaMemcpyDeviceToHost),__LINE__);
		
		errorHandler(cudaFree(devPtr_forces),__LINE__);
		errorHandler(cudaFree(devPtr_velocities),__LINE__);
		
		
		//temp_current /= 3*npart;
		float sigma = sqrt(temp_target);
		
		for(i=0; i<npart; i++) {
			if(rand()/(float)RAND_MAX < nu*dt) {
				for(j=0; j<3; j++)
					velocities[3*i+j] = gauss(sigma,0);
			
			}
		}
	
	}
}*/

#else
void calcForce(float* en,float* f, const float* positions, const int npart, const float *box, const float cutoff2, const int n, const float ecut) {
	
	int i,j;
	
	*en = 0;
	
	for(i=0; i <3*npart; i++) {
		f[i] = 0;
	}
	
	float rx;
	float ry;
	float rz;
	float r2;
	float r2i;
	float r6i;
	float ff; //lennard-jones potential
	
	for(i=0; i < npart-1; i++) {
		for(j=i+1; j<npart; j++) {
				
				rx = positions[3*i] - positions[3*j];
				ry = positions[3*i+1] - positions[3*j+1];
				rz = positions[3*i+2] - positions[3*j+2];
				rx -= box[0]* round(rx/box[0]);
				ry -= box[1]* round(ry/box[1]);
				rz -= box[2]* round(rz/box[2]);
				
				r2 = rx*rx + ry*ry + rz*rz;
				
				if(r2 < 0.1) fprintf(stderr,"Not good! i %d j %d distance %f frame %d\n",i,j,r2,n);
				if(r2 < cutoff2) {
					r2i = 1/r2;
					r6i = r2i*r2i*r2i;
					ff = 48*r2i*r6i*(r6i - 0.5);
					f[3*i] += ff*rx;
					f[3*j] -= ff*rx;
					f[3*i+1] += ff*ry;
					f[3*j+1] -= ff*ry;
					f[3*i+2] += ff*rz;
					f[3*j+2] -= ff*rz;
					
					*en += 4*r6i*(r6i-1) - ecut;
				}
			
		}
	}
}


/*void integrateVelVerlet(const int part, const float* f, const int npart, float* positions, float* velocities, const float dt, const float temp_target, const float nu) {
	
	int i,j;
	
	if(part ==1) {
		
		for(i=0; i<npart; i++)
			for(j=0; j<3; j++) {
				positions[3*i+j] = positions[3*i+j] + dt* velocities[3*i+j] +dt*dt* f[3*i+j]/2;
				velocities[3*i+j] += dt*f[3*i+j]/2;
			}
				
	}
	
	if(part ==2) {
		
		//float temp_current =0;
		for(i=0; i<npart; i++) {
			for(j=0; j<3; j++){
				velocities[3*i+j] = velocities[3*i+j] +dt* f[3*i+j]/2;
				//temp_current += velocities[3*i+j] *velocities[3*i+j];
			}
		}
		
		//temp_current /= 3*npart;
		float sigma = sqrt(temp_target);
		
		for(i=0; i<npart; i++) {
			if(rand()/(float)RAND_MAX < nu*dt) {
				for(j=0; j<3; j++)
					velocities[3*i+j] = gauss(sigma,0);
			
			}
		}
	
	}
}*/
#endif

void integrateVelVerlet(const float en, const float* box, const int part, const float* f, const int npart, float* positions, float* velocities, const float dt, const float temp_target, const float nu, const int frame) {
	
	int i,j;
	
	if(part ==1) {
		
		for(i=0; i<npart; i++)
			for(j=0; j<3; j++) {
				
				positions[3*i+j] += dt* velocities[3*i+j] +dt*dt* f[3*i+j]/2; //update position
				
				//apply periodic boundary conditions
				if(positions[3*i+j] > box[j]/2 || positions[3*i+j] < -box[j]/2) {
					//float pposition = positions[3*i+j];
					//positions[3*i+j] -= floor(positions[3*i+j]/box[j])*box[j];
					positions[3*i+j] -= round(positions[3*i+j]/box[j])*box[j];
					//fprintf(stdout,"Out of box positive particle %d in frame %d. Calculated position %f, corrected position %f\n",i,frame,pposition, positions[3*i+j]);
				}

				
				velocities[3*i+j] += dt*f[3*i+j]/2;
			}	
				
	}
	
	if(part ==2) {
		
		float temp_current =0;
		for(i=0; i<npart; i++) {
			for(j=0; j<3; j++){
				velocities[3*i+j] = velocities[3*i+j] +dt* f[3*i+j]/2;
				temp_current += velocities[3*i+j] *velocities[3*i+j];
			}
		}
		
		temp_current /= 3*npart;
		float sigma = sqrt(temp_target);
		
		float* sumv = (float*)calloc(3,sizeof(float));
		for(i=0; i<npart; i++) {
			if(rand()/(float)RAND_MAX < nu*dt) {
				for(j=0; j<3; j++) {
					velocities[3*i+j] = gauss(sigma,0);
					sumv[j] += velocities[3*i+j];
				}
			}
		}
		
		sumv[0] /=npart;
		sumv[1] /=npart;
		sumv[2] /=npart;
		
		float sumvSquared = sumv[0]*sumv[0] + sumv[1]*sumv[1] + sumv[2]*sumv[2];	
		float etot = (en + 0.5*temp_current)/npart;
		
		fprintf(stdout,"Vel CMS: %f Temp: %f Energy: %f frame: %d\n",sumvSquared, temp_current, etot, frame);
	
		
		//fprintf(stdout,"Velocity Center of mass: %f frame: %d\n",sumv2, frame);
	
	}
}


void integrateVerlet(const float *box, const float* f, const float en, const int npart, float* positions, float* ppositions, float* velocities, const float dt, const int frame) {
	
	int i,j;
	
	float* sumv = (float*) calloc(3,sizeof(float));
	float sumv2 = 0;
	float xx = 0;
	
	for(i=0; i<npart; i++) {
		for(j=0; j<3; j++) {
			xx = 2*positions[3*i+j] - ppositions[3*i+j] + dt*dt*f[3*i+j];
			velocities[3*i+j] = (xx - ppositions[3*i+j])/(2*dt);
			sumv[j] += velocities[3*i+j];
			sumv2 += velocities[3*i+j]*velocities[3*i+j];
			ppositions[3*i+j] = positions[3*i+j];
			positions[3*i+j] = xx;
			
	
		}
	}
	
		sumv[0] /=npart;
		sumv[1] /=npart;
		sumv[2] /=npart;
		
		float sumvSquared = sumv[0]*sumv[0] + sumv[1]*sumv[1] + sumv[2]*sumv[2];
		float temp = sumv2/(3*npart);
		float etot = (en + 0.5*sumv2)/npart;
		
		fprintf(stdout,"Vel CMS: %f Temp: %f Energy: %f frame: %d\n",sumvSquared, temp, etot, frame);
	
	
	free(sumv);
}

//void sample(FILE* xres, FILE* vres, FILE* fres, int npart, float* positions, float* velocities, float* f) {
void sample(FILE* xres, FILE* vres, int npart, float* positions, float* velocities) {
	
	int i;
	
	for(i=0; i<npart; i++) {
		//fprintf(xres, " %d %8.8f %8.8f %8.8f \n", i,positions[3*i],positions[3*i+1],positions[3*i+2]);
		fprintf(xres, " atom %8.8f %8.8f %8.8f \n", positions[3*i],positions[3*i+1],positions[3*i+2]);
	}
	
	for(i=0; i<npart; i++) {
		//fprintf(vres, " %d %8.8f %8.8f %8.8f \n", i,velocities[3*i],velocities[3*i+1],velocities[3*i+2]);
		fprintf(vres, " atom %8.8f %8.8f %8.8f \n",velocities[3*i],velocities[3*i+1],velocities[3*i+2]);
		
	}
}

/*void sample(FILE* xres, FILE* vres, int npart, float* positions, float* velocities, float box, int nslabs, int* density) {
	
	int i;
	
	for(i=0; i<npart; i++) {
		//fprintf(xres, " %d %8.8f %8.8f %8.8f \n", i,positions[3*i],positions[3*i+1],positions[3*i+2]);
		fprintf(xres, " atom %8.8f %8.8f %8.8f \n", positions[3*i],positions[3*i+1],positions[3*i+2]);
	}
	
	for(i=0; i<npart; i++) {
		//fprintf(vres, " %d %8.8f %8.8f %8.8f \n", i,velocities[3*i],velocities[3*i+1],velocities[3*i+2]);
		fprintf(vres, " atom %8.8f %8.8f %8.8f \n",velocities[3*i],velocities[3*i+1],velocities[3*i+2]);
	}
	
		//Dichte Messung - liefert Teilchenzahl pro Schicht in jedem Durchgang als XYZ Datei
		//Eindimensional & in Z-Richtung
	for (i=0; i<npart; i++) {
		density[(int)(positions[3*i+2] / (box / nslabs))] += 1;  //Verteile Atome auf Slabs
	}
	
	
	/* for(i=0; i<npart; i++) {
		//fprintf(fres, " %d %8.8f %8.8f %8.8f \n", i,f[3*i],f[3*i+1],f[3*i+2]);
		fprintf(fres, " atom %8.8f %8.8f %8.8f \n",f[3*i],f[3*i+1],f[3*i+2]);
	}
	
}*/

int main(int argc, char* argv[])
{
	int npart,nslabs;
	float temp,dt,tmax,cutoff, cutoff2,nu;
	float box[3];
	//cutoff for lennard jones
	//char* xpath;
	//char* vpath;
	
	char* inXpath;
	char* inVpath;
	//char* fpath;
	int thermostat =0;
	int continued = 0;
	//int densityMeasurement = 0;
		
	if(argc==9) {
		
			sscanf(argv[1],"%d",&npart); //units: https://en.wikipedia.org/wiki/Lennard-Jones_potential#Dimensionless_.28reduced.29_units
			sscanf(argv[2],"%f",&temp);
			sscanf(argv[3],"%f",&dt);
			sscanf(argv[4],"%f",&tmax);
			sscanf(argv[5],"%f",&box[0]);
			sscanf(argv[6],"%f",&box[1]);
			sscanf(argv[7],"%f",&box[2]);
			sscanf(argv[8],"%f",&cutoff);
			//xpath = argv[9];
			//vpath = argv[10];
			//endXpath = argv[10];
			//endVpath = argv[11];
			printf("Input: Particles: %d Temperature: %f timestep: %f Maxtime: %f BoxX: %f BoxY: %f BoxZ: %f cutoff: %f\n",npart,temp,dt,tmax,box[0],box[1],box[2],cutoff);
		} else if(argc==10) {
			
			sscanf(argv[1],"%d",&npart); //units: https://en.wikipedia.org/wiki/Lennard-Jones_potential#Dimensionless_.28reduced.29_units
			sscanf(argv[2],"%f",&temp);
			sscanf(argv[3],"%f",&dt);
			sscanf(argv[4],"%f",&tmax);
			sscanf(argv[5],"%f",&box[0]);
			sscanf(argv[6],"%f",&box[1]);
			sscanf(argv[7],"%f",&box[2]);
			sscanf(argv[8],"%f",&cutoff);
			//xpath = argv[9];
			//vpath = argv[10];
			//endXpath = argv[11];
			//endVpath = argv[12];
			sscanf(argv[9],"%f",&nu);
			thermostat =1;
			printf("Input: Particles: %d Temperature: %f timestep: %f Maxtime: %f BoxX: %f BoxY: %f BoxZ: %f cutoff: %f collision-probability: %f run with thermostat\n",npart,temp,dt,tmax,box[0],box[1],box[2],cutoff,nu);
		} else if(argc==11) {
			
			//sscanf(argv[1],"%d",&npart); //units: https://en.wikipedia.org/wiki/Lennard-Jones_potential#Dimensionless_.28reduced.29_units
			sscanf(argv[1],"%f",&temp);
			sscanf(argv[2],"%f",&dt);
			sscanf(argv[3],"%f",&tmax);
			sscanf(argv[4],"%f",&box[0]);
			sscanf(argv[5],"%f",&box[1]);
			sscanf(argv[6],"%f",&box[2]);
			sscanf(argv[7],"%f",&cutoff);
			//xpath = argv[8];
			//vpath = argv[9];
			sscanf(argv[8],"%f",&nu);
			//endXpath = argv[11];
			//endVpath = argv[12];
			inXpath = argv[9];
			inVpath = argv[10];
			
			thermostat =1;
			continued =1;
			printf("Input: Temperature: %f timestep: %f Maxtime: %f BoxX: %f BoxY: %f BoxZ: %f cutoff: %f collision-probability: %f continued run with thermostat reading path: %s reading path: %s\n",temp,dt,tmax,box[0],box[1],box[2],cutoff, nu, inXpath, inVpath);
		/*} else if(argc==12) {
			
			//sscanf(argv[1],"%d",&npart); //units: https://en.wikipedia.org/wiki/Lennard-Jones_potential#Dimensionless_.28reduced.29_units
			sscanf(argv[1],"%f",&temp);
			sscanf(argv[2],"%f",&dt);
			sscanf(argv[3],"%f",&tmax);
			sscanf(argv[4],"%f",&box[0]);
			sscanf(argv[5],"%f",&box[1]);
			sscanf(argv[6],"%f",&box[2]);
			sscanf(argv[7],"%f",&cutoff);
			//xpath = argv[8];
			//vpath = argv[9];
			sscanf(argv[8],"%f",&nu);
			//endXpath = argv[11];
			//endVpath = argv[12];
			inXpath = argv[9];
			inVpath = argv[10];
			sscanf(argv[11],"%d",&nslabs);
			
			thermostat =1;
			continued =1;
			densityMeasurement=1;
			printf("Input: Temperature: %f timestep: %f Maxtime: %f BoxX: %f BoxY: %f BoxZ: %f cutoff: %f collision-probability: %f continued run with thermostat reading path: %s reading path: %s Number of slabs %d \n",temp,dt,tmax,box[0],box[1],box[2],cutoff, nu, inXpath, inVpath,nslabs);*/
		} else {
			
			printf("Syntax: ./MD (*not in continued run* <number of particles>) <temperature> <dt> <tmax> <boxX> <boxY> <boxZ> <cutoffDistance> (*optional* <collision probability> <Readpath pos> <Readpath vel>)\n");
			exit(-1);
		}
		
	cutoff2 = cutoff*cutoff;
	FILE* xres = fopen("./pos.xyz","w");
	FILE* vres = fopen("./vel.xyz","w");
//	FILE* fres = fopen(fpath,"w");
	
	//FILE* dres = fopen("./density.xyz","w");
	//if(!dres) {printf("File not found!");}
	
	
	if(!xres || !vres) {
		printf("File not found!");
		exit(-1);
	}
	
	//int* density = (int*) calloc(nslabs,sizeof(int));
 	
	float en;
	
	float *positions =NULL;
	float *ppositions = NULL;
	float *velocities = NULL;
	float *f = NULL;
	
	if(continued) { 		
		
		FILE* xin = fopen(inXpath,"r");
		FILE* vin = fopen(inVpath,"r");
		
		if(!xin || !vin) {
			printf("File not found!");
			exit(-1);
		}
		
		initContinued(&npart, &positions, &velocities, &f,xin,vin);
		fprintf(stderr,"pos %f %f %f\n", positions[0],positions[1],positions[2]);
		fclose(xin);
		fclose(vin);

	} else {
		initStart(npart, temp, &positions, &ppositions, &velocities, &f,dt,box);
	}
	
	fprintf(xres,"%d\n",npart);
	fprintf(xres,"generated by my MD simulation\n");
	fprintf(vres,"%d\n",npart);
	fprintf(vres,"generated by my MD simulation\n");
	
	float t = 0;
	
	float cutoffi = 1/cutoff2;
	float cutoff6i = cutoffi*cutoffi*cutoffi;
	float ecut = 4*cutoff6i*(cutoff6i-1);

	
	if(thermostat) calcForce(&en,f, positions, npart, box, cutoff2,(int)(t/dt), ecut);
	
	sample(xres, vres, npart, positions, velocities);
	while(t < tmax) {
		
		if(!thermostat) {
			calcForce(&en, f, positions, npart, box, cutoff2,(int)(t/dt), ecut);
			integrateVerlet(box, f, en, npart, positions, ppositions, velocities, dt,(int)(t/dt));
			t += dt;
			sample(xres, vres, npart, positions, velocities);
			fprintf(stderr,"Stage %% %f\r",(double)t/tmax*100.0);
			
		} else {
			
			integrateVelVerlet(en,box, 1,f, npart, positions, velocities, dt, temp,nu,(int)(t/dt));
			calcForce(&en, f, positions, npart, box, cutoff2,(int)(t/dt), ecut);
			integrateVelVerlet(en,box, 2,f, npart, positions, velocities, dt, temp,nu,(int)(t/dt));
			t += dt;
			if(t/dt==10000) sample(xres, vres, npart, positions, velocities);
			//if(!densityMeasurement) sample(xres, vres, npart, positions, velocities);
				//else sample(xres, vres, npart, positions, velocities,box[2],nslabs,density);
			fprintf(stderr,"Stage %% %f\r",(double)t/tmax*100.0);
		}
	}
	
	/*if(densityMeasurement) {
		fprintf(dres, "Number of particles in cells of length %f / %d \n", box[2],nslabs);
		
		for(int i=0; i<nslabs;i++){							
			fprintf(dres, "Slab %d % d \n", i, density[i] / ( (int) (tmax/dt) ) );						//Gebe Anzahl Atome pro Slab in einer Zeile nebeneinander aus 
		}
		fclose(dres);
	}*/
	
	fclose(xres);
	fclose(vres);
	//fclose(fres);
	free(f);
	
	FILE* xend = fopen("./posEnd.xyz","w");
	FILE* vend = fopen("./velEnd.xyz","w");
	
	if(!xend || !vend) {
		printf("File not found!");
		exit(-1);
	}
	
	fprintf(xend,"%d\n",npart);
	fprintf(xend,"generated by my MD simulation\n");
	fprintf(vend,"%d\n",npart);
	fprintf(vend,"generated by my MD simulation\n");
	
	sample(xend, vend, npart, positions, velocities);
	
	fclose(xend);
	fclose(vend);
	
	free(positions);
	free(ppositions);
	free(velocities);
	
}
