/**Nelder-Mead Method to Search Minimun for Multivariate Function**/
#ifndef _NMSIMPLEX_H_
#define _NMSIMPLEX_H_

#define		ALPHA		1.0			//反射点线段比例
#define		GAMMA		1.0			//扩展点线段比例
#define		BETA		0.5			//压缩点线段比例
#define		SCALE		0.5			//初始点到各个单纯形顶点的步长
#define		MATSIZ		50			//函数变量个数
#define		MAXIT		600			//最大迭代次数
#define		VERR		1.0e-8		//计算允许误差

struct PT {			//单纯形定点数据结构
	float f;		//顶点函数值
	float x[MATSIZ];//顶点坐标
};
typedef float (*PFN)(float*);		//函数指针类型定义（最末*表示参数为数组）

class NMSimplex {
private:
	PT v[MATSIZ+1], vm, vr, ve, vc;	//单纯形顶点，中点，反射点，扩展点和收缩点
	int vs;							//最小值顶点索引
	int vg;							//最大值顶点索引
	int vh;							//第2最大值顶点索引
	PFN fn;							//函数指针
	int n;							//变量个数
 	float *xmin;					//最小值点变量数组结果

	void init();					//初始化
	void calcVertex( PT * );				//计算单纯形顶点函数值
public:
	void findVertexIndexes();		//查找S,H,G顶点的索引
//计算特殊点M,R,E,C
	void calcCenterPt();			//中点
	void calcReflectPt();			//反射点
	void calcExtendPt();			//扩展点
	void calcConstractPt();			//压缩点
	void shrink();					//收缩计算
	bool checkout();					//判别迭代终止
	void simplex(float *fm);						//主算法

//构造析构函数
	NMSimplex() {}
	NMSimplex( PFN f, float *start, int num):	//start为函数变量初始估计值
		fn(f), n(num), xmin(start){
		init();	}
	~NMSimplex() {}

};

#endif