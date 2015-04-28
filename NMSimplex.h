/**Nelder-Mead Method to Search Minimun for Multivariate Function**/
#ifndef _NMSIMPLEX_H_
#define _NMSIMPLEX_H_

#define		ALPHA		1.0			//������߶α���
#define		GAMMA		1.0			//��չ���߶α���
#define		BETA		0.5			//ѹ�����߶α���
#define		SCALE		0.5			//��ʼ�㵽���������ζ���Ĳ���
#define		MATSIZ		50			//������������
#define		MAXIT		600			//����������
#define		VERR		1.0e-8		//�����������

struct PT {			//�����ζ������ݽṹ
	float f;		//���㺯��ֵ
	float x[MATSIZ];//��������
};
typedef float (*PFN)(float*);		//����ָ�����Ͷ��壨��ĩ*��ʾ����Ϊ���飩

class NMSimplex {
private:
	PT v[MATSIZ+1], vm, vr, ve, vc;	//�����ζ��㣬�е㣬����㣬��չ���������
	int vs;							//��Сֵ��������
	int vg;							//���ֵ��������
	int vh;							//��2���ֵ��������
	PFN fn;							//����ָ��
	int n;							//��������
 	float *xmin;					//��Сֵ�����������

	void init();					//��ʼ��
	void calcVertex( PT * );				//���㵥���ζ��㺯��ֵ
public:
	void findVertexIndexes();		//����S,H,G���������
//���������M,R,E,C
	void calcCenterPt();			//�е�
	void calcReflectPt();			//�����
	void calcExtendPt();			//��չ��
	void calcConstractPt();			//ѹ����
	void shrink();					//��������
	bool checkout();					//�б������ֹ
	void simplex(float *fm);						//���㷨

//������������
	NMSimplex() {}
	NMSimplex( PFN f, float *start, int num):	//startΪ����������ʼ����ֵ
		fn(f), n(num), xmin(start){
		init();	}
	~NMSimplex() {}

};

#endif