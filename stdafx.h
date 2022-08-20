// stdafx.h : ��׼ϵͳ�����ļ��İ����ļ���
// ���Ǿ���ʹ�õ��������ĵ�
// �ض�����Ŀ�İ����ļ�
//

#pragma once

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>
#include<map>


// TODO:  �ڴ˴����ó�����Ҫ������ͷ�ļ�
#include "struct.h"
#include "lambda.h"
#include<string>
using namespace std;

//ʱ��ϵͳת��
int CommonTimeToMJDTime(int &year, int &month, int &day, int &hour, int &minute, double &second, int &days, double &fracday);
int MJDTimeToCommonTime(int &days, double &fracday, int &year, int &month, int &day, int &hour, int &minute, double &second);
int GPSTimeToMJDTime(int &week, double &secondofweek, int &days, double &fracday);
int MJDTimeToGPSTime(int &days, double &fracday, int &week, double &secondofweek);
int CommonTimeToGPSTime(int &year, int &month, int &day, int &hour, int &minute, double &second, int &week, double &secondofweek);

//����ϵת��
void BLHToXYZ(double a, double f, BLH *blh, XYZ *xyz);
int XYZToBLH(double a, double f, BLH *blh, XYZ *xyz);
void XYZToNEU(XYZ delta, XYZ StaPos_ref, double *n, double *e, double *u);

//��������
int MatrixAdd(int line_a, int row_a, int line_b, int row_b, double a[], double b[], double c[]);
int MatrixSub(int line_a, int row_a, int line_b, int row_b, double a[], double b[], double c[]);
int MatrixMul(int line_a, int row_a, int line_b, int row_b, double a[], double b[], double c[]);
int MatrixInv(int n, double a[], double b[]);
int MatrixTrans(int line_a, int row_a, double a[], double b[]);

//���ļ�
int ReadHeader_RENIX(ifstream &fin, OBS *Obs);
int ReadObs_RENIX(ifstream &fin, OBS *Obs);
int ReadEph_RENIX(ifstream &fin, map<int, map<GPSTIME, EPHEM>> &map_eph);
string Remove_space(const string& src);
void InitObs(OBS *Obs);

//spp������
void ComputeAzimuthElement(double a, double f, XYZ *RecPos, XYZ *SatPos, double *azimuth, double *element);
int IonoFree(XYZ *RecPos, OBS *Obs, CALCULATION *Calculation);
double Hopfield(XYZ *RecPos, double element);
void EarthRotate(XYZ *Pos, double &tao);

//spp
int GPS_SatPositionVelocity(GPSTIME t, int prn, CALCULATION *Calculation, map<int, map<GPSTIME, EPHEM>> map_eph, double &tao);
int GPS_SatClock(GPSTIME t_sv, int prn, CALCULATION *Calculation, map<int, map<GPSTIME, EPHEM>> map_eph);
int BDS_SatPositionVelocity(GPSTIME t, int prn, CALCULATION *Calculation, map<int, map<GPSTIME, EPHEM>> map_eph, double &tao);
int BDS_SatClock(GPSTIME t_sv, int prn, CALCULATION *Calculation, map<int, map<GPSTIME, EPHEM>> map_eph);
int GPS_SignalTransmit(OBS Obs, CALCULATION *Calculation, map<int, map<GPSTIME, EPHEM>> map_eph);
int BDS_SignalTransmit(OBS Obs, CALCULATION *Calculation, map<int, map<GPSTIME, EPHEM>> map_eph);
void DeleteZero(double *arr, int n);
int SPP(OBS Obs, SPP_CAL *SPP_Cal, map<int, map<GPSTIME, EPHEM>> map_eph, CALCULATION *Calculation);
void PreEdit(OBS *Obs, map<int, map<GPSTIME, EPHEM>> map_eph);

//RTK
int RTK(OBS Obs_1, OBS Obs_2, CALCULATION Calculation, RTK_CAL *RTK_Cal);
int SPV(CALCULATION *Calculation, OBS Obs, XYZ *RecPos, XYZ *RecVel, XYZ *sigma2_v);