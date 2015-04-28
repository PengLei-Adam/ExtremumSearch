#include <math.h>
#include "NMSimplex.h"

//��ʼ�������������ʼ�������飬�����ʼֵ������㣬�������е㺯��ֵ������
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
//���㶥�㺯��ֵ
void NMSimplex::calcVertex( PT * pt ){
	pt ->f = (*fn)( &(pt->x[1]));
}
//����S,H,G��������
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
//�����е�
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
//���㷴���
void NMSimplex::calcReflectPt() {
	for( int j = 0; j < n; j++ )
		vr.x[j] = (1 + ALPHA)* vm.x[j] - ALPHA* v[vg].x[j];
	calcVertex( &vr);
}
//������չ��
void NMSimplex::calcExtendPt() {
	for( int j = 0; j < n; j++ )
		ve.x[j] = (1 + GAMMA)* vr.x[j] - GAMMA* vm.x[j];
	calcVertex( &ve);
}
//����ѹ����
void NMSimplex::calcConstractPt() {
	for( int j = 0; j < n; j++ )
		vc.x[j] = BETA* v[vg].x[j] + (1 - BETA)* vm.x[j];
	calcVertex( &vc );
}
 
//�������䡢��չ��ѹ����û���ҵ���G��С�ĵ㣬Ҫ������������
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
//ͨ���������б������ֹ
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
//���㷨
void NMSimplex::simplex( float *fmin ){		//����Ϊ�����ַ
	int i, j, m;
	for( m = 0; m < MAXIT; m++ ) {
		findVertexIndexes();		//����G,H,S������
		calcCenterPt();				//����M��
		calcReflectPt();			//���㷴���R
		if( vr.f > v[vs].f && vr.f <= v[vh].f )
			v[vg] = vr;				//��R���滻G�㣨=���أ���
		if( vr.f <= v[vs].f ) {
			calcExtendPt();			//������չ��E
			v[vg] = ve.f <= v[vs].f ? ve : vr;	//E��R�н�С���滻G��
		}
		if( vr.f > v[vh].f ) {
			calcConstractPt();		//����ѹ����C
			if( vc.f < v[vg].f ) v[vg] = vc;	//��C���滻G��
			else shrink();
		}
		if( checkout() ) break;
	}
	vs = 0;
	for( i = 0; i < n+1; i++ )
		if( v[i].f < v[vs].f ) vs = i;	//Ѱ����Сֵ��
	for( j = 0; j < n; j++ )
		xmin[j] = v[vs].x[j];			//���б�������Ϊ��Сֵ����
	*fmin = (*fn)(xmin);
}