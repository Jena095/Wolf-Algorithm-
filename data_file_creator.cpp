# include <iostream>
# include <cstdlib>
# include <cstdio>
# include <cmath>

# define sigma		16.0
# define rho		45.92
# define beta		4.0
# define step_len   1.25e-3


class trajectory{

public:

	double x,y,z;

public:

	void incrementor();
	double distance_calculator(trajectory*);  
};


void trajectory :: incrementor(){

	double kx1, kx2, kx3, kx4, ky1, ky2, ky3, ky4, kz1, kz2, kz3, kz4;


	//FILE* fp;
	//fp = fopen("trsajectory_data","w+");

	kx1 = sigma*(this->y - this->x);
	ky1 = this->x*(rho - this->z) - this->y;
	kz1 = this->x * this->y - beta* this->z;


	kx2 = sigma*((this->y + ky1*step_len*0.5)   - (this->x + kx1*step_len*0.5));
	ky2 = (this->x + kx1*step_len*0.5) * (rho - (this->z + kz1*step_len*0.5)) - (this->y + ky1*step_len*0.5);
	kz2 = (this->x + kx1*step_len*0.5) * (this->y +ky1*step_len*0.5) - beta*(this->z + kz1*step_len*0.5);


	kx3 = sigma*((this->y + ky2*step_len*0.5)   - (this->x + kx2*step_len*0.5));  
	ky3 = (this->x + kx2*step_len*0.5) * (rho - (this->z + kz2*step_len*0.5)) - (this->y + ky2*step_len*0.5);
	kz3 = (this->x + kx2*step_len*0.5) * (this->y +ky2*step_len*0.5) - beta*(this->z + kz2*step_len*0.5);


	kx4 = sigma*((this->y + ky3*step_len)   - (this->x + kx3*step_len));  
	ky4 = (this->x + kx3*step_len) * (rho - (this->z + kz3*step_len)) - (this->y + ky3*step_len);
	kz4 = (this->x + kx3*step_len) * (this->y +ky3*step_len) - beta*(this->z + kz3*step_len);


	this->x  = this->x + step_len * (kx1 + 2*(kx2 + kx3) + kx4 ) * 0.166666667;
	this->y  = this->y + step_len * (ky1 + 2*(ky2 + ky3) + ky4 ) * 0.166666667;
	this->z  = this->z + step_len * (kz1 + 2*(kz2 + kz3) + kz4 ) * 0.166666667;

	return;
}


double trajectory :: distance_calculator(trajectory* T){
	return (sqrt((this->x - T->x)*(this->x - T->x) + (this->y - T->y)*(this->y - T->y) + (this->z - T->z)*(this->z - T->z)));
}


int main(){

	trajectory* one = new trajectory();

	double d_0,ddd,d_1,lambda,a,lambda_avg;

	int flag = 0;

	/*FILE* fp;
	fp = fopen("loren_data", "w+");
*/

	one->x = -5.60;
	one->y = -9.0;
	one->z = 16.0;
	a = 1.0;

	char filename_x[1000];

	
	for(int j = 0; j < 20; j++)
	{


		FILE* fpx;

		sprintf(filename_x,"loren_data_x_%d",j);

		fpx = fopen(filename_x, "w+");

		for(int i = 0; i < 1000; i++){
			one->incrementor();
		}


		for(long i = 1; i < 100000; i++)
		{
				one->incrementor();
				fprintf(fpx,"%lf\n",one->x);

		}

		fclose(fpx);
	}

	
	delete(one);

	return (0);

}







