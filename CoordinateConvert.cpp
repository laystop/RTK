#include "stdafx.h"
#include <iostream>
#include<cmath>
#include "struct.h"
using namespace std;

/*************************************
BLHToXYZ
Ŀ�ģ�������굽�ѿ��������ת���㷨

������
a   ���򳤰���
f   �������
blh �������ϵ����
xyz ecef����ϵ����
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
Ŀ�ģ��ѿ������굽��������ת���㷨

������
a   ���򳤰���
f   �������
B   ��ؾ��ȣ���λ�����ȣ�
blh �������ϵ����
xyz ecef����ϵ����

����ֵ��1=���� 0=��ѭ��
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
		cout << "XYZ2BLH������ѭ����" << endl;
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
Ŀ�ģ��ռ�ֱ�����굽վ�������ת���㷨

������
delta        ����
StaPos_ref   վ��xyz����
n            վ������n
e            վ������e
u            վ������u
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