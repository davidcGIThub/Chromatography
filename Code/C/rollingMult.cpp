#include <iostream>

class RollingMult
{
    public:
        double* func(double array[], double factor[])
	{
		int len = sizeof(array)/sizeof(*array);
		int window = sizeof(factor)/sizeof(*factor);
		double arr[len-window+1];
		for(int i = 0; i <= (len-window); i++)
		{
			int sum = 0;
			for(int j = i; j < (i+window); j++)
			{	
				sum = sum + array[j] * factor[j-i];
			}	
			arr[i] = sum;
		}	
		return arr;
        }
};


/*extern "C" 
{
	RollingMult* RollingMult_new()
	{ 
		return new RollingMult(); 
	}
    
	void RollingMult_func(RollingMult* rollingMult)
	{ 
		rollingMult->func(double*, double); 
	}
}
*/

int main()
{
	RollingMult roll = RollingMult();
	double tempvec[] = {1.0,2.0,3,4,5,6,7,8,9,10};
	double kern[] = {1,2,3};
	double* a = roll.func(tempvec,kern);
	for(int i = 0; i < sizeof(a)/sizeof(*a); i++)
	{
		std::cout << a[i] << std::endl;
	} 
}


