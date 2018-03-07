#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <omp.h>
#include <string>
#include <cstdlib>
#include <math.h>


/*
// t1--> time separation in components of vector  v0={x0,x1,x2,...xd}
// t2--> time separtion  in vectors  v0   v1   v2  v3
// t3-->  time for transient

// distance calculation is replaced by ((x-xi)^2 +( y-yi)^2 + (z- zi)^2 + ...)^(0.5)


i have changed certain condition in the requirement like
if adj_index + t3 is grater than row
if the min len is greater than a given cutoff limit
and i am planning to add other requirements soon

date--- 04/11/2017

date----03/01/2018

the parameter t1 and t2 are equated as suggested by Mr Soumitro bannerjee 
and this is how it was described in the wolf algorithm. I further implement the 
additional condition of angle criteria in this code.

int t2 = 0;
int t3 = 0;
*/


//int t1 = 0;

//# define  t1  60
//# define  t2  60
# define  t3  		40
# define  tol 		1e-1
# define  scale_min 2.0
# define  scale_max 75.0
# define  angle_max 0.3
//# define dimension 7

using namespace std;

void data_input(vector<double> &vec, char filename[1000])
{

	ifstream in_stream;

	ofstream out_stream;

	double tmp, average;

	//Check input and output file.

	in_stream.open(filename);

	if (in_stream.fail())
	{
		cout << "Input file opening failed";

		exit(1);
	}

	while (in_stream >> tmp)

	{
		vec.push_back(tmp);
	}
}

void correalation_cal(vector<double> &vec)

{
	double average, tmp;

	for (int i = 0; i < vec.size(); i++)
	{
		average = average + tmp;
	}

	for (int i = 0; i < vec.size(); i++)
	{
		vec[i] = vec[i] - average;
	}

	double numerator, denominator;

	int size = vec.size();

	FILE *fp;

	fp = fopen("correlation.dat", "w+");

	for (int i = 0; i < 1000; i++)
	{
		numerator = 0.0;

		denominator = 0.0;

		for (int j = 0; j + i < size; j++)
		{
			numerator = numerator + vec[j] * vec[j + i];

			denominator = denominator + vec[j] * vec[j];
		}

		fprintf(fp, "%lf\t\t%d\n", (numerator / denominator), i);
	}
}

void vector_creator(double **arr, int &dimension, long int &row, vector<double> &vec,int t1)
{

	/*for(long int i = 0; i < vec.size() - t1*dimension; i = i + t1)
	{
		cout << i << "\t\t" << vec[i] << endl;

	}

	exit(1);*/

	#pragma omp parallel for schedule(dynamic) 
	for (long int i = 0; i < row; i++)
	{
		arr[i] = new double[dimension];

		for (int j = 0; j < dimension; j++)
		{
			//cout << "works fine" << "\t\trow:" << row << "\t\ti:" << i << "\t\t" << "\t\tj:" << j << "\t\t" <<i*t2 + t1*j << endl;

			arr[i][j] = vec[(i * t1) + ( t1 * j)];

			//cout << arr[i][j] << "\t\t";
		}
	}

	//exit(1);

	return;
}
/*
double min_distance_calc(double **arr, int &dimension, long int &row, int &index, int &adj_trajec_index)
{
	double len, tmp_var;

	len = 1000000000.0;
	tmp_var = 0;

	for (long int i = 0; i < row; i++)
	{
		if (i != index)
		{
			for (int j = 0; j < dimension; j++)
			{
				/*if(tmp_var < fabs(arr[index][j] - arr[i][j]))
				{
					tmp_var = fabs(arr[index][j] - arr[i][j]) ;
					//cout << i<<"\n";
				}
				
				tmp_var = tmp_var + (arr[index][j] - arr[i][j]) * (arr[index][j] - arr[i][j]);
			}

			tmp_var = sqrt(tmp_var);

			if (len > tmp_var && len > 0.01)
			{
				len = tmp_var;

				adj_trajec_index = i;
			}
		}
	}

	return (len);
}

double distance_calc(double **arr, int &dimension, long int &row, int &index, int &adj_trajec_index)
{
	double tmp_var = 0;

	for (int j = 0; j < dimension; j++)
	{
		tmp_var = tmp_var + (arr[index][j] - arr[adj_trajec_index][j]) * (arr[index][j] - arr[adj_trajec_index][j]);
	}

	tmp_var = sqrt(tmp_var);
	return (tmp_var);
}

void lyapunov_expo(double **arr, int &dimension, long int &row, vector<double> &lyp)
{
	double len_min, len_max, lyp_exp;

	int base_index = 0, adj_trajec_index;

	long int prop_time = row;

	double time_step = 0;

	lyp_exp = 0.0;

	FILE *fp;


	fp = fopen("lyp_data", "a+");

	for (long int i = 0; i < prop_time; i = i + t3)
	{

		base_index = i;

		len_min = min_distance_calc(arr, dimension, row, base_index, adj_trajec_index);

		cout << prop_time << "\t\t" << i << endl;

		if (adj_trajec_index + t3 < row && base_index + t3 < row)
		{

			base_index = base_index + t3;
			adj_trajec_index = adj_trajec_index + t3;
			len_max = distance_calc(arr, dimension, row, base_index, adj_trajec_index);
		}

		if (len_max >= len_min)
		{

			lyp_exp = lyp_exp + (log(len_max / len_min)); // * (1.0/(t3*time_step*t2*1.0e-2)) ;
			time_step = time_step + 1.0;
		}
	}

	lyp_exp = lyp_exp * (1.0 / (t3 * time_step * t2 * 1e-3));
	lyp.push_back(lyp_exp);
	cout << "lyp_exp_avg: " << lyp_exp << "\t\t"
		 << "time_step:  " << time_step << endl; //<<"prop_time:  "<<prop_time <<endl;

	//	cout<<"\nlyp_exp:"<<lyp_exp<<endl;

	fprintf(fp, "%lf\n", lyp_exp);
	fclose(fp);
}
*/
//the function min distance calc only calculates the minimum distance between the two vecotrs
double min_distance_calc(double **arr, int &dimension, long int &row, int &index, int &adj_trajec_index)
{
	double len, tmp_var;

	len = 1000000000.0;
	tmp_var = 0;

	for (long int i = 0; i < row; i++)
	{
		if (i != index)
		{
			for (int j = 0; j < dimension; j++)
			{
				/*if(tmp_var < fabs(arr[index][j] - arr[i][j]))
				{
					tmp_var = fabs(arr[index][j] - arr[i][j]) ;
					//cout << i<<"\n";
				}
				*/
				tmp_var = tmp_var + (arr[index][j] - arr[i][j]) * (arr[index][j] - arr[i][j]);
			}

			tmp_var = sqrt(tmp_var);

			if (len > tmp_var && len > scale_min)
			{
				len = tmp_var;

				adj_trajec_index = i;
			}
		}
	}

	return (len);
}

double min_dis_calc(double **arr, int &dimension, long int &row, int &base_index, int &adj_trajec_index,vector<int> & discarded_vec_point)
{
	double len, tmp_var;

	int flag = 0;

	len = 1000000000.0;
	tmp_var = 0.0;

	for (long int i = 0; i < row; i++)
	{
		flag = 0;
		for(long int t = 0; t < discarded_vec_point.size(); t++)
		{
			if( i == discarded_vec_point[t] )
			{
				flag = 1;
				break;
			}
		}

		if ( flag == 0)
		{
			for (int j = 0; j < dimension; j++)
			{
				/*if(tmp_var < fabs(arr[index][j] - arr[i][j]))
				{
					tmp_var = fabs(arr[index][j] - arr[i][j]) ;
					//cout << i<<"\n";
				}
				*/
				tmp_var = tmp_var + (arr[base_index][j] - arr[i][j]) * (arr[base_index][j] - arr[i][j]);
			}

			tmp_var = sqrt(tmp_var);

			if (len > tmp_var && len > scale_min)
			{
				len = tmp_var;

				adj_trajec_index = i;
			}
		}
	}

	return (len);
}

bool angle_chk_func(double **arr, double *angle_arr,int & dimension, int &base_index,int &adj_trajec_index)
{
	double mod_angle_arr = 0.0, mod_new_angle_arr = 0.0, dot_product = 0.0;
	for(int i = 0; i < dimension; i++)
	{
		mod_angle_arr = mod_angle_arr + angle_arr[i] * angle_arr[i];

		mod_new_angle_arr = mod_new_angle_arr + (arr[adj_trajec_index][i] - arr[base_index][i]) * (arr[adj_trajec_index][i] - arr[base_index][i]) ;

		dot_product = dot_product + (angle_arr[i]) * (arr[adj_trajec_index][i] - arr[base_index][i]) ;
	}

	mod_new_angle_arr = sqrt(mod_new_angle_arr);
	mod_angle_arr     = sqrt(mod_angle_arr); 

	//cout << "the angle found is :  " << acos( dot_product / (mod_angle_arr * mod_new_angle_arr )) << endl;

	if ( acos( dot_product / (mod_angle_arr * mod_new_angle_arr )) <= angle_max )
	{
		return(true);
	}
	else
	{
		return(false);
	}
}

double distance_calc(double **arr, int &dimension, long int &row, int &index, int &adj_trajec_index)
{
	double tmp_var = 0;

	for (int j = 0; j < dimension; j++)
	{
		tmp_var = tmp_var + (arr[index][j] - arr[adj_trajec_index][j]) * (arr[index][j] - arr[adj_trajec_index][j]);
	}

	tmp_var = sqrt(tmp_var);
	return (tmp_var);
}

void lyapunov_expo(double **arr, int &dimension, long int &row, vector<double> &lyp,int t1)
{
	double len_min, len_max, lyp_exp = 0.0, time_step = 0.0;

	int base_index = 0, adj_trajec_index;
	long int counter;

	vector<int> discarded_vec_point;
	double angle_arr[dimension];

	FILE *fp;

	fp = fopen("lyp_data", "a+");

	// the loop controls propagation of vetors in time

	for (long int i = 0; i < row; i = i + t3)
	{
		counter = 0;
		//cout << counter <<"\t\trow:"<< i << endl;
		base_index = i;
		discarded_vec_point.push_back(base_index);

	// for the first time do not consider the angle criteria.
		if(i == 0)
		{
			//base_index = i;
			len_min = min_distance_calc(arr, dimension, row, base_index, adj_trajec_index);
		}
	// once the minimum distance has been calculated for the first time we move on to consider the minimum angle criteria
		else
		{
	// the code below checks both the minimum distance and minimum angle criteria.
			while(true)
			{
				len_min = min_dis_calc(arr,dimension, row, base_index, adj_trajec_index,discarded_vec_point);
				//cout << len_min << endl;

				if(angle_chk_func(arr,angle_arr,dimension,base_index,adj_trajec_index))
				{
					discarded_vec_point.clear();
					break;
				}
				else
				{
					counter += 1;
					//cout << "i am in" << "\t" << counter << "\t" << row <<endl;
					if(counter >= row)
					{
						//cout << "breakup" << endl;
						discarded_vec_point.clear();
						break;
					}
					//cout << counter << endl ;
					discarded_vec_point.push_back(adj_trajec_index);
				}
	
			}
		}

		//cout << "i am out" << endl;

	//cout << prop_time << "\t\t" << i << endl;

	// the if statement below checks whether further propagation is possible or not

		if (adj_trajec_index + t3 < row && base_index + t3 < row)
		{
	// the code below propagates the vecotr in time for t3 units
			
			base_index = base_index + t3;
			adj_trajec_index = adj_trajec_index + t3;
			for(int D = 0 ; D < dimension; D++)
			{
				angle_arr[D] = arr[adj_trajec_index][D] - arr[base_index][D];
			}

	// the code below calculates the distance between the two vectors after they have propagated in time
			len_max = distance_calc(arr, dimension, row, base_index, adj_trajec_index);
		}

	// the code below checks whether the vectors separated in time " caution -------- why shall i choose that the len min is greater than len max condition ? 
	//and what will happen if i carry out the calculation without using this condition"
		if (len_max >= len_min && len_max < scale_max)
		{

			lyp_exp = lyp_exp + (log(len_max / len_min)); // * (1.0/(t3*time_step*t2*1.0e-2)) ;
			time_step = time_step + 1.0;
		}
	}


	// when all the lyapunov exponent has been calculated we take the average to compute the avg lypunov exponent.

	lyp_exp = lyp_exp * (1.0 / (t3 * time_step * t1 * 1.25e-3));

	// store the avg lyapunov exponent in a vector so that we can do statistics on the data set. We use many data file and store the lypunov exponent from each of the data file

	lyp.push_back(lyp_exp);

	cout << "lyp_exp_avg: " << lyp_exp << "\t\t"
		 << "time_step:  " << time_step << endl; //<<"prop_time:  "<<prop_time <<endl;

	//	cout<<"\nlyp_exp:"<<lyp_exp<<endl;

	// printing the data obtained in a file

	fprintf(fp, "%lf\n", lyp_exp);


	fclose(fp);
}

int main()
{
	/*vector<double> vec;
	vector<double> lyp;*/
	int dimension = 3;

	/*long int row = 0;char filename[1000];*/
	


	FILE *fp;

	fp = fopen("t3_40_shift_1.25.dat", "a+");

	#pragma omp parallel for schedule(dynamic)
			for (int t1 = 25; t1 <36; t1 = t1 + 1)
			{
				//cout << "begin" << endl;
				//int t2 = t1;
				vector<double> vec;
				vector<double> lyp;
				long int row = 0;
				char filename[1000];

		//#pragma omp parallel for schedule(dynamic) 
				for (int j = 0; j < 20; j++)
				{

					sprintf(filename, "loren_data_x_%d", j);

			//data input convertes the data stored in the file to a 1 dimensional vector.
					data_input(vec, filename);

					long int size = vec.size() - 1;

					cout << "vec.size:" << size << endl;
					//	correalation_cal(vec);

			// calculating the number of rows(v0,v1,v2...) of dimension (2d+1).
					row = floor((size - (t1*dimension)));
					row = (row / (t1 - 1));

			//creating the 2d array to store all the vectors created.
					double **arr = new double *[row];

			//using the vector creator fuction to input elements from 1D data vector to 2D array
					vector_creator(arr, dimension, row, vec,t1);

			//using the lyapunov expo function i calculate the lyapunov exponent using the wolf algorithm
					lyapunov_expo(arr, dimension, row, lyp, t1);

					for (long int i = 0; i < row; i++)
					{
						delete[] arr[i];
					}

					delete[] arr;

					vec.clear();
				}

				//cout << "\n\nlyp.size():\t" << lyp.size() << endl;

				double lyp_expo_avg = 0.0;
				double lyp_expo_avg_sqre = 0.0;
				double list_size = lyp.size();

				if(list_size > 0){

					for (long int i = 0; i < list_size; i++)
					{

						lyp_expo_avg = lyp_expo_avg + lyp[i];

						lyp_expo_avg_sqre = lyp_expo_avg_sqre + lyp[i] * lyp[i];
					}

					lyp_expo_avg = lyp_expo_avg / (list_size);
				
					lyp_expo_avg_sqre = lyp_expo_avg_sqre / (list_size);
				
					cout << "\navg__lyp__expo:" << lyp_expo_avg << "\ts_dev: " << fabs(sqrt(lyp_expo_avg_sqre - (lyp_expo_avg * lyp_expo_avg)) / lyp_expo_avg) << endl;
					if(!isnan(lyp_expo_avg))
					{
						//#pragma omp critical
						fprintf(fp, "%d\t\t%lf\t\t\t%lf\n", t1, lyp_expo_avg, fabs(sqrt(lyp_expo_avg_sqre - (lyp_expo_avg * lyp_expo_avg)) / lyp_expo_avg));
					}
				}

				lyp.clear();
			}

	return 0;
}
