#include <iostream>
#include <math.h>

#define PHI (-1.0+sqrt(5.0))/2.0		//define golden ratio

typedef float (*PFN)(float x);

float gss( PFN f, float x0, float x1, float crit ) {	//golden section search[x0,x1]
	float x2 = x0 + (x1-x0)* PHI;
	float x3;
	while( fabs(x0 - x1) > crit ) {				//区间小于误差限
		x3 = x0 + (x2-x0)* PHI;					//x3取距x0较远的黄金分割点
		if( (*f)(x3) < (*f)(x2) )				
		{		x2 = x3;  x1 = x2; }			//以[x0,x3,x2]为新区间
		else
		{		x0 = x1;  x1 = x3; }			//x0,x1大小关系颠倒,以[x3,x2,x1]为新区间
	}
	return (x0+x1)/2.0;
}


float fibs ( PFN f, float x0, float x1, int n )	{	//Fibnoacci search
	float fibn[100], x2, x3;
	int k;
	fibn[0] = 0;
	fibn[1] = 1;
	for( k = 2; k <= n; k++ )
		fibn[k] = fibn[k - 1] + fibn[k - 2];		//creat fibonacci sequence
	x2 = x0 + (x1 - x0)* fibn[n-1]/ fibn[n];
	for( k = 1; k < n-2; k++ ) {
		x3 = x0 + (x2-x1)* fibn[n-k-1]/fibn[n-k];
		if( (*f)(x3) < (*f)(x2) ) {
			x1 = x2;
			x2 = x3;
		} else {
			x0 = x1;
			x1 = x3;
		}
	}
	return x3;
}