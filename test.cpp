#include<iostream>
#include"NMSimplex.h"
using namespace std;

float myfunc( float* x ){	//自定义函数
	float fx = ( x[0] * x[0] -13 )* (x[0] * x[0] - 13 );
	return fx;
}
int main() {
	float fmin;
	float start[1] = { 13 };
	NMSimplex nms( myfunc, start, 1 );
	nms.simplex( &fmin );
	cout<< start[0] <<"\n" << fmin <<endl;
	return 0;
}
