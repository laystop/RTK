#include "stdafx.h"
#include "struct.h"
#include<iostream>
#include<cmath>
using namespace std;

/*************************************
ComputeAzimuthElement
Ŀ�ģ������λ�Ǻ͸߶Ƚ�

������
a        ���򳤰���
f        �������
RecPos   ��վxyz����
SatPos   ����xyz����
azimuth  ��λ��
element  �߶Ƚ�
**************************************/
void ComputeAzimuthElement(double a, double f, XYZ *RecPos, XYZ *SatPos, double *azimuth, double *element)
{
	double N, E, U;
	BLH Rec;
	XYZToBLH(a, f, &Rec, RecPos);
	N = -sin(Rec.b)*cos(Rec.l)*(SatPos->x - RecPos->x) - sin(Rec.b)*sin(Rec.l)*(SatPos->y - RecPos->y) + cos(Rec.b)*(SatPos->z - RecPos->z);
	E = -sin(Rec.l)*(SatPos->x - RecPos->x) + cos(Rec.l)*(SatPos->y - RecPos->y);
	U = cos(Rec.b)*cos(Rec.l)*(SatPos->x - RecPos->x) + cos(Rec.b)*sin(Rec.l)*(SatPos->y - RecPos->y) + sin(Rec.b)*(SatPos->z - RecPos->z);
	*azimuth = atan2(E, N);
	*element = atan2(U, sqrt(N*N + E*E));
}

/**************************************************************
IonoFree
Ŀ�ģ�˫Ƶ�û�����������

������
RecPos       ���ջ�λ��
RawData      ԭʼ����
Calculation  ������

����ֵ��0=������ 1=�ɹ�
***************************************************************/
int IonoFree(XYZ *RecPos, OBS *Obs, CALCULATION *Calculation)
{
	double k_GPS, k_BDS;
	int i, prn;
	BLH Rec;

	XYZToBLH(CGCS2000_a, CGCS2000_f, &Rec, RecPos);
	if (Rec.h<1e-8 || Rec.h>1e6)
	{
		return 0;
	}
	else
	{
		k_GPS = GPS_L1*GPS_L1 / GPS_L2 / GPS_L2;
		k_BDS = BDS_B1*BDS_B1 / BDS_B3 / BDS_B3;

		for (i = 0; i < Obs->GPS_SatNum; i++)
		{
			prn = Obs->GPS_Sat[i].Prn;
			if (Obs->GPS_Sat[i].Psr[0]>1 && Obs->GPS_Sat[i].Psr[1]>1) //˫Ƶ�۲�ֵ����
			{
				if (abs(Obs->GPS_Sat[i].Psr[0] - Obs->GPS_Sat[i].Psr[1])>100) //����α������
					Calculation->GPS_POSnVEL[prn - 1].iono_flag = 1;
				else
				{
					Calculation->GPS_POSnVEL[prn - 1].Psr_InonCorrect = (Obs->GPS_Sat[i].Psr[1] - k_GPS*Obs->GPS_Sat[i].Psr[0]) / (1 - k_GPS);
					Calculation->GPS_POSnVEL[prn - 1].iono_flag = 2;
				}
			}
			else
				Calculation->GPS_POSnVEL[prn - 1].iono_flag = 3;
		}
		for (i = 0; i < Obs->BDS_SatNum; i++)
		{
			prn = Obs->BDS_Sat[i].Prn;
			if (Obs->BDS_Sat[i].Psr[0]>1 && Obs->BDS_Sat[i].Psr[1] > 1) //˫Ƶ�۲�ֵ����
			{
				if (abs(Obs->BDS_Sat[i].Psr[0] - Obs->BDS_Sat[i].Psr[1])>100)
					Calculation->BDS_POSnVEL[prn - 1].iono_flag = 1;
				else
				{
					Calculation->BDS_POSnVEL[prn - 1].Psr_InonCorrect = (Obs->BDS_Sat[i].Psr[1] - k_BDS*Obs->BDS_Sat[i].Psr[0]) / (1 - k_BDS);
					Calculation->BDS_POSnVEL[prn - 1].iono_flag = 2;
				}
			}
			else
				Calculation->BDS_POSnVEL[prn - 1].iono_flag = 3;
		}
		return 1;
	}
}

/**************************************************************
Hopfield
Ŀ�ģ����������

������
RecPos   ��վλ��
element  �߶Ƚ�

����ֵ�����������ֵ 0=����ʧ��
***************************************************************/
double Hopfield(XYZ *RecPos, double element)
{
	double H0, T0, p0, RH0, hd, hw, Kd, Kw, e, T, p, RH;
	BLH Rec;
	if ((RecPos->x*RecPos->x + RecPos->y*RecPos->y + RecPos->z*RecPos->z) < 1e-8)
	{
		Rec.b = 0.0;
		Rec.l = 0.0;
		Rec.h = -CGCS2000_a;
	}
	else
		XYZToBLH(CGCS2000_a, CGCS2000_f, &Rec, RecPos);
	H0 = 0.0;
	T0 = 15 + 273.16;
	p0 = 1013.25;
	RH0 = 0.5;
	if (Rec.h<1e-8 || Rec.h>1.7e4)
	{
		return 0;
	}
	else
	{
		T = T0 - 0.0065*(Rec.h - H0);
		p = p0*pow((1 - 0.0000226*(Rec.h - H0)), 5.225);
		RH = RH0*exp(-0.0006396*(Rec.h - H0));
		hd = 40136 + 148.72*(T0 - 273.16);
		hw = 11000.0;
		Kd = 155.2*1e-7*p / T*(hd - Rec.h);
		e = RH*exp(-37.2465 + 0.213166*T - 0.000256908*T*T);
		Kw = 155.2*1e-7 * 4810 / T / T*e*(hw - Rec.h);
		element = element*180.0 / PI;
		return Kd / sin(sqrt(element*element + 6.25) / 180.0*PI) + Kw / sin(sqrt(element*element + 2.25) / 180.0*PI);
	}
}

/**************************************************************
EarthRotate
Ŀ�ģ�������ת����

������
RecPos   ��վλ��
tao      �źŴ���ʱ��
***************************************************************/
void EarthRotate(XYZ *Pos, double &tao)
{
	Pos->x = cos(tao*EarthRotation)*Pos->x + sin(tao*EarthRotation)*Pos->y;
	Pos->y = -sin(tao*EarthRotation)*Pos->x + cos(tao*EarthRotation)*Pos->y;
	Pos->z = Pos->z;
}