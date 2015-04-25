#include <iostream>
#include <math.h>

#define PHI (-1.0+sqrt(5.0))/2.0		//define golden ratio

typedef float (*PFN)(float x);

float gss( PFN f, float x0, float x1, float crit ) {	//golden section search[x0,x1]
	float x2 = x0 + (x1-x0)* PHI;
	float x3;
	while( fabs(x0 - x1) > crit ) {				//����С�������
		x3 = x0 + (x2-x0)* PHI;					//x3ȡ��x0��Զ�Ļƽ�ָ��
		if( (*f)(x3) < (*f)(x2) )				
		{		x2 = x3;  x1 = x2; }			//��[x0,x3,x2]Ϊ������
		else
		{		x0 = x1;  x1 = x3; }			//x0,x1��С��ϵ�ߵ�,��[x3,x2,x1]Ϊ������
	}
	return (x0+x1)/2.0;
}


