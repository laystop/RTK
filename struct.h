#pragma once

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>
#include <map>

#define CGCS2000_a 6378137
#define CGCS2000_f 1/298.257222101
#define PI 3.14159265357989

//#define MAXRAWLEN 40960 //??????
#define GPSsatnum 32
#define BDSsatnum 63
/*
#define POLYCRC32 0xEDB88320u
#define MsgID_Obs 43
#define MsgID_GPSEph 7
#define MsgID_BDSEph 1696
#define MsgID_ion 8
*/
enum NAVSYS { GPS, BDS };

#define SpeedofLight 2.99792458e8
#define GPS_miu 3.986005e14
#define GPS_omegae_dot 7.2921151467e-5
#define BDS_miu 3.986004418e14
#define BDS_omegae_dot 7.2921150e-5

#define GPS_L1 1.57542e9
#define GPS_L2 1.2276e9
#define BDS_B1 1.561098e9
#define BDS_B3 1.26852e9

#define lambda_L1 SpeedofLight/GPS_L1
#define lambda_L2 SpeedofLight/GPS_L2
#define lambda_B1 SpeedofLight/BDS_B1
#define lambda_B3 SpeedofLight/BDS_B3

#define EarthRotation 7.292115e-5 

#define SatComMax 60 //��๲�����Ǹ���

/*ʱ��ϵͳ*/
struct GPSTIME {
	int Week;
	double SecofWeek;

	GPSTIME()
	{
		Week = 0;
		SecofWeek = 0.0;
	}

	bool operator< (const GPSTIME &other) const
	{
		return Week * 604800 + SecofWeek < (other.Week * 604800 + other.SecofWeek);
		//����С�ںţ�����map������ʽ������Ϊ�˴Ӵ�С���򣬰�<�ĳ�>
	}
};

struct COMMONTIME {
	int Year;
	int Month;
	int Day;
	int Hour;
	int Minute;
	double Second;
};

/*�ѿ�������ϵ����*/
struct XYZ {
	double x;
	double y;
	double z;
};

/*�������ϵ����*/
struct BLH {
	double b;
	double l;
	double h;
};

/*���ݽ���*/
struct SAT {
	unsigned short Prn;
	//NAVSYS Sys;
	double Psr[2];
	double Adr[2]; // BDS, Adr1=B1, Adr2=B3; GPS, Adr1=L1, Adr2=L2
	double Dop[2];
	double snr[2], LockTime[2];
	//float sigma_psr[2], sigma_adr[2];
};

struct OBS {
	SAT GPS_Sat[GPSsatnum], BDS_Sat[BDSsatnum];
	COMMONTIME ObsTime_ymd;
	GPSTIME ObsTime;
	int GPS_SatNum, BDS_SatNum, SatNum;
	int GPS_L1C_index, GPS_L1P_index, GPS_L1D_index, GPS_L2C_index, GPS_L2P_index, GPS_L2D_index;
	int BDS_B1C_index, BDS_B1P_index, BDS_B1D_index, BDS_B3C_index, BDS_B3P_index, BDS_B3D_index;
	XYZ StaPos_ref;
};

struct EPHEM {
	unsigned short PRN, iodc; //iode? ������aode��aodc����ʾ����������ʱ�䣺��Сʱ
	NAVSYS Sys;
	GPSTIME EphTime;
	double tow, toe, sqrt_A, deltaN, M0;
	double ecc, w, cuc, cus, crc, crs, cic, cis;
	double Iup0, Idown0, w0, wdot, toc, tgd, i0;
	double af0, af1, af2, URA, IDOT, omega, omega_0, omega_dot;
	unsigned long IODE, Zweek, week;

	bool operator< (const EPHEM &other) const
	{
		return PRN <other.PRN;
	}
};

/*
struct RAWDATA {
	//GPSTIME StartTime[3], EndTime[3]; //�ļ���ʼ������ʱ��
	OBS Obs[2];  //����վ
	//EPHEM GPSEph[GPSsatnum], BDSEph[BDSsatnum];
	map<int, map<GPSTIME, EPHEM>> map_eph;  //GPS 32, BDS 63, prn��һ��95����GPS:n, BDS:32+n��
};
*/

/*��������λ���ٶ�*/
struct SAT_POSnVEL {
	XYZ SatPos; //����λ��
	double v_x, v_y, v_z; //�����ٶ�
	double delta_tsv, delta_tsv_dot; //�����Ӳ������
	GPSTIME SignalTrTime; //���Ƿ���ʱ�䣨GPS BDS����GPST��
	double Psr_InonCorrect, TroCorrect; //����㡢���������
	int iono_flag;
	double Element; //�߶Ƚ�
	double sigma2;  //���ģ��
};

struct CALCULATION {
	SAT_POSnVEL GPS_POSnVEL[GPSsatnum], BDS_POSnVEL[BDSsatnum];
};

struct SPP_CAL {
	XYZ RecPos;
	BLH RecBlh;
	GPSTIME time;
	double delta_t_receive[2], delta_t_dot[2];
	double sigma, PDOP;
	int SatNum;
};

struct RTK_CAL {
	GPSTIME time;
	XYZ StaPos_ref, StaPos_flow_tmp, StaPos_flow;
	XYZ float_line, fixed_line; //����
	double Qxx[9], Qxx_fixed[9]; //Э��������
	double ratio;
	int SatNum;
};