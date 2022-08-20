#include "stdafx.h"
#include <iostream>
#include<cmath>
#include "struct.h"
using namespace std;

/*************************************
BLHToXYZ
目的：大地坐标到笛卡尔坐标的转换算法

参数：
a   椭球长半轴
f   椭球扁率
blh 大地坐标系坐标
xyz ecef坐标系坐标
**************************************/
void BLHToXYZ(double a, double f, BLH *blh, XYZ *xyz)
{
	double N, e2, B, L, H;
	B = blh->b;
	L = blh->l;
	H = blh->h;

	e2 = 2 * f - f*f;
	N = a / (sqrt(1 - e2*sin(B)*sin(B)));
	xyz->x = (N + H)*cos(B)*cos(L);
	xyz->y = (N + H)*cos(B)*sin(L);
	xyz->z = (N*(1 - e2) + H)*sin(B);
}

/*************************************
XYZToBLH
目的：笛卡尔坐标到大地坐标的转换算法

参数：
a   椭球长半轴
f   椭球扁率
B   大地经度（单位：弧度）
blh 大地坐标系坐标
xyz ecef坐标系坐标

返回值：1=正常 0=死循环
**************************************/
int XYZToBLH(double a, double f, BLH *blh, XYZ *xyz)
{
	int j;
	double delta_Z[2], N, e2, B, L, H, X, Y, Z;
	X = xyz->x;
	Y = xyz->y;
	Z = xyz->z;

	if ((X * X + Y * Y + Z * Z) < 1e-8)
	{
		blh->b = 0.0;
		blh->l = 0.0;
		blh->h = -a;
		return 0;
	}
	e2 = 2 * f - f*f;
	delta_Z[0] = e2*Z;

	blh->l = atan2(Y, X);
	for (j = 0; j < 10; j++)
	{
		blh->b = atan2(Z + delta_Z[0], sqrt(X*X + Y*Y));
		N = a / (sqrt(1 - e2*sin(blh->b)*sin(blh->b)));
		blh->h = sqrt(X*X + Y*Y + (Z + delta_Z[0])*(Z + delta_Z[0])) - N;
		delta_Z[1] = N*e2*(Z + delta_Z[0]) / (sqrt(X*X + Y*Y + (Z + delta_Z[0])*(Z + delta_Z[0])));

		if (abs(delta_Z[1] - delta_Z[0]) < 1e-10) break;
		else delta_Z[0] = delta_Z[1];
	}

	if (abs(delta_Z[1] - delta_Z[0]) >= 1e-10)
	{
		cout << "XYZ2BLH进入死循环！" << endl;
		return 0;
	}
	else
	{
		blh->b = atan2(Z + delta_Z[1], sqrt(X*X + Y*Y));
		N = a / (sqrt(1 - e2*sin(blh->b)*sin(blh->b)));
		blh->h = sqrt(X*X + Y*Y + (Z + delta_Z[1])*(Z + delta_Z[1])) - N;
		return 1;
	}
}

/*************************************
XYZToNEU
目的：空间直角坐标到站心坐标的转换算法

参数：
delta        基线
StaPos_ref   站心xyz坐标
n            站心坐标n
e            站心坐标e
u            站心坐标u
**************************************/
void XYZToNEU(XYZ delta, XYZ StaPos_ref, double *n, double *e, double *u)
{
	BLH blh_ref;
	double b, l, h;
	XYZToBLH(CGCS2000_a, CGCS2000_f, &blh_ref, &StaPos_ref);
	b = blh_ref.b;
	l = blh_ref.l;
	h = blh_ref.h;
	*n = -sin(b)*cos(l)*delta.x - sin(b)*sin(l)*delta.y + cos(b)*delta.z;
	*e = -sin(l)*delta.x + cos(l)*delta.y;
	*u = cos(b)*cos(l)*delta.x + cos(b)*sin(l)*delta.y + sin(b)*delta.z;
}