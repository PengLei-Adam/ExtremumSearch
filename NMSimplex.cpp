#include <math.h>
#include "NMSimplex.h"

//初始化函数，输入初始变量数组，计算初始值存入各点，计算所有点函数值并存入
void NMSimplex::init( ) {
	float pn, qn;
	int i,j;

	pn = SCALE* ( sqrt(float(n+1)) -1+n )/( n*sqrt(2.) );
	qn = SCALE* ( sqrt(float(n+1)) -1 )/( n*sqrt(2.) );
	for( j = 0; j < n; j++ ) v[0].x[j] = xmin[j];
	for( i = 1; i < n+1; i++ ) {
		for( j = 0; j < n; j++ ) {
			if( i -1 == j )		v[i].x[j] = pn + xmin[j];
			else		v[i].x[j] = qn + xmin[j];
		}
	}
	for( i = 0; i < n+1; i++ )
		calcVertex( &v[i] );
}
//计算顶点函数值
void NMSimplex::calcVertex( PT * pt ){
	pt ->f = (*fn)( &(pt->x[1]));
}
//计算S,H,G顶点索引
void NMSimplex::findVertexIndexes() {
	int i, j;
	vs = vg = 0;
	for( i = 0; i < n+1; i++ ) {
		if( v[i].f < v[vs].f ) vs = i;
		if( v[i].f > v[vg].f ) vg = i;
	}
	vh = 0;
	for( j = 0; j < n+1; j ++ )
		if( v[i].f > v[vh].f && v[i].f < v[vg].f ) 
			vh = i;
}
//计算中点
void NMSimplex::calcCenterPt() {
	int i, j;
	float cent;
	for( j = 0; j < n; j++ ){
		cent = 0.0;
		for( i = 0; i < n+1; i++ )
			if( i != vg )  cent += v[i].x[j];
		vm.x[j] = cent/n;
	}
}
//计算反射点
void NMSimplex::calcReflectPt() {
	for( int j = 0; j < n; j++ )
		vr.x[j] = (1 + ALPHA)* vm.x[j] - ALPHA* v[vg].x[j];
	calcVertex( &vr);
}
//计算扩展点
void NMSimplex::calcExtendPt() {
	for( int j = 0; j < n; j++ )
		ve.x[j] = (1 + GAMMA)* vr.x[j] - GAMMA* vm.x[j];
	calcVertex( &ve);
}
//计算压缩点
void NMSimplex::calcConstractPt() {
	for( int j = 0; j < n; j++ )
		vc.x[j] = BETA* v[vg].x[j] + (1 - BETA)* vm.x[j];
	calcVertex( &vc );
}
 
//经过反射、扩展和压缩，没有找到比G点小的点，要进行收缩计算
void NMSimplex::shrink() {
	int i, j;
	for( i = 0; i < n+1; i++ ){
		if( i != vs ){
			for( j = 0; j < n; j++ )
				v[i].x[j] = v[vs].x[j] +( v[i].x[j] - v[vs].x[j] )/2.0;
			calcVertex(&v[i]);
		}
	}
}
//通过均方差判别迭代终止
bool NMSimplex::checkout() {
	int i;
	float fsum, favg, s;
	fsum = 0.0;
	for( i = 0; i < n+1; i++ )
		fsum = fsum + v[i].f;
	favg = fsum/( n+1.);
	s = 0.0;
	for( i = 0; i < n+1; i++ )
		s += ( v[i].f - favg )*( v[i].f - favg )/n;
	s = sqrt(s);
	return  (s<VERR ? true:false);
}
//主算法
void NMSimplex::simplex( float *fmin ){		//参数为结果地址
	int i, j, m;
	for( m = 0; m < MAXIT; m++ ) {
		findVertexIndexes();		//查找G,H,S点索引
		calcCenterPt();				//计算M点
		calcReflectPt();			//计算反射点R
		if( vr.f > v[vs].f && vr.f <= v[vh].f )
			v[vg] = vr;				//用R点替换G点（=重载？）
		if( vr.f <= v[vs].f ) {
			calcExtendPt();			//计算扩展点E
			v[vg] = ve.f <= v[vs].f ? ve : vr;	//E、R中较小点替换G点
		}
		if( vr.f > v[vh].f ) {
			calcConstractPt();		//计算压缩点C
			if( vc.f < v[vg].f ) v[vg] = vc;	//用C点替换G点
			else shrink();
		}
		if( checkout() ) break;
	}
	vs = 0;
	for( i = 0; i < n+1; i++ )
		if( v[i].f < v[vs].f ) vs = i;	//寻找最小值点
	for( j = 0; j < n; j++ )
		xmin[j] = v[vs].x[j];			//所有变量更新为最小值变量
	*fmin = (*fn)(xmin);
}