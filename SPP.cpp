#include "stdafx.h"
#include "struct.h"
#include<iostream>
#include<cmath>
#include<map>
#include<algorithm>
#include<fstream>
#include<string>
#include<iomanip>

using namespace std;

/**************************************************************
DeleteZero
目的：删除矩阵中等于0的元素

参数：
arr   待处理的数组
n     数组维数
**************************************************************/
void DeleteZero(double *arr, int n)
{
	//直接从头到尾遍历整个数组，是0就将n减1，同时数组下表往前移
	for (int i = 0; i < n; i++)
	{
		if (abs(arr[i]) < 1e-10)
		{
			//将所有元素向前移
			int j = i;
			for (; j < n - 1; j++)
			{
				arr[j] = arr[j + 1];
			}
			//将要判断的元素的下标前移
			i = i - 1;
			//非0元素个数减1
			n--;
		}
	}
}

/**************************************************************
PreEdit
目的：对观测文件进行数据处理，如果观测到星历里不存在的卫星，就删除

参数：


**************************************************************/
void PreEdit(OBS *Obs, map<int, map<GPSTIME, EPHEM>> map_eph)
{
	int i, prn;
	for (i = 0; i < Obs->GPS_SatNum; i++)
	{
		prn = Obs->GPS_Sat[i].Prn;
		if (map_eph.find(prn) == map_eph.end())
		{
			Obs->GPS_Sat[i] = Obs->GPS_Sat[Obs->GPS_SatNum - 1];
			Obs->GPS_SatNum--;
			i--;
		}
	}
	for (i = 0; i < Obs->BDS_SatNum; i++)
	{
		prn = Obs->BDS_Sat[i].Prn + 32;
		if (map_eph.find(prn) == map_eph.end())
		{
			Obs->BDS_Sat[i] = Obs->BDS_Sat[Obs->BDS_SatNum - 1];
			Obs->BDS_SatNum--;
			i--;
		}
	}
}

/**************************************************************
SPP
目的：单点定位

参数：
RawData              原始数据
Rcv                  接收机计算结果
Calculation          计算结果

返回值：0=卫星数目不足，无法计算  1=正常
***************************************************************/
int SPP(OBS Obs, SPP_CAL *SPP_Cal, map<int, map<GPSTIME, EPHEM>> map_eph, CALCULATION *Calculation)
{
	int i, j, k, prn, n_GPS, n_BDS, n;
	double azimuth, element, iono_correct, tropo_correct, rho, rho_dot, position, position_new, tao;
	double v[30], B[150], B_T[150], P_value[30], P[900], BTP[150], w[30], BB_T[25], BB_T1[25], BBTP[150], BB_T2[150], BB_T3[150], x[5], v_T[30]; //最多接收30颗卫星
	double  Vv[30], VB[150], VB_T[150], VP[900], VBTP[150], Vw[30], VBB_T[25], VBB_T1[25], VBBTP[150], VBB_T2[150], Vx[5], Vv_T[30];
	double vtpv[1];
	double sigma, PDOP;
	XYZ RecTemp;
	double delta_t_receive[2], delta_t_dot[2];

	azimuth = 0.0;
	element = 0.0;
	RecTemp.x = SPP_Cal->RecPos.x;
	RecTemp.y = SPP_Cal->RecPos.y;
	RecTemp.z = SPP_Cal->RecPos.z;
	position = sqrt(pow(RecTemp.x, 2) + pow(RecTemp.y, 2) + pow(RecTemp.z, 2));
	delta_t_receive[0] = SPP_Cal->delta_t_receive[0];
	delta_t_receive[1] = SPP_Cal->delta_t_receive[1];
	delta_t_dot[0] = SPP_Cal->delta_t_dot[0];
	delta_t_dot[1] = SPP_Cal->delta_t_dot[1];
	memset(P, 0, sizeof(P)); //对高度角加权
	SPP_Cal->time = Obs.ObsTime;

	//计算卫星信号发射时间
	GPS_SignalTransmit(Obs, Calculation, map_eph);
	BDS_SignalTransmit(Obs, Calculation, map_eph);

	//单点定位
	for (j = 0; j < 10; j++)
	{
		IonoFree(&RecTemp, &Obs, Calculation);

		n_GPS = 0; //可用于计算的卫星数量（同时有星历、观测数据且星历未过期）
		for (i = 0; i < Obs.GPS_SatNum; i++) //GPS
		{
			prn = Obs.GPS_Sat[i].Prn;
			tao = Obs.ObsTime.SecofWeek - Calculation->GPS_POSnVEL[prn - 1].SignalTrTime.SecofWeek;
			if (Calculation->GPS_POSnVEL[prn - 1].iono_flag == 1) //双频伪距之差过大
				continue;
			if (GPS_SatPositionVelocity(Calculation->GPS_POSnVEL[prn - 1].SignalTrTime, prn, Calculation, map_eph, tao) == 2) //计算卫星位置速度
				continue;
			ComputeAzimuthElement(CGCS2000_a, CGCS2000_f, &RecTemp, &Calculation->GPS_POSnVEL[prn - 1].SatPos, &azimuth, &element); //计算高度角方位角（更新）

			if (abs(SPP_Cal->RecPos.x) < 1e-10 || element > 15 / 180.0*PI)
			{
				/////////
				if (abs(SPP_Cal->RecPos.x) < 1e-10 || azimuth > 0)
				{
					tropo_correct = Hopfield(&RecTemp, element); //计算对流层改正（更新）
					Calculation->GPS_POSnVEL[prn - 1].TroCorrect = tropo_correct;
					Calculation->GPS_POSnVEL[prn - 1].Element = element / PI*180.0;
					Calculation->GPS_POSnVEL[prn - 1].sigma2 = 0.04*0.04 + 0.03*0.03 / sin(element) / sin(element);
					//Calculation->GPS_POSnVEL[prn - 1].sigma2 = 0.003*0.003*(1 + 1.5*cos(element)*cos(element));

					rho = sqrt(pow((Calculation->GPS_POSnVEL[prn - 1].SatPos.x - RecTemp.x), 2) + pow((Calculation->GPS_POSnVEL[prn - 1].SatPos.y - RecTemp.y), 2) + pow((Calculation->GPS_POSnVEL[prn - 1].SatPos.z - RecTemp.z), 2));

					B[n_GPS * 5] = (RecTemp.x - Calculation->GPS_POSnVEL[prn - 1].SatPos.x) / rho;
					B[n_GPS * 5 + 1] = (RecTemp.y - Calculation->GPS_POSnVEL[prn - 1].SatPos.y) / rho;
					B[n_GPS * 5 + 2] = (RecTemp.z - Calculation->GPS_POSnVEL[prn - 1].SatPos.z) / rho;
					B[n_GPS * 5 + 3] = 1;
					B[n_GPS * 5 + 4] = 0;
					//P_value[n_GPS] = sin(element)*sin(element);
					//P_value[n_GPS] = 1 / (0.03*0.03 + 0.04*0.04 / sin(element) / sin(element));
					P_value[n_GPS] = 1 / (pow(0.13 + 0.56*exp(-element / PI * 180 / 10), 2));

					if (Calculation->GPS_POSnVEL[prn - 1].iono_flag == 2)
						w[n_GPS] = Calculation->GPS_POSnVEL[prn - 1].Psr_InonCorrect - rho - delta_t_receive[0] + SpeedofLight*Calculation->GPS_POSnVEL[prn - 1].delta_tsv - tropo_correct;
					else if (Obs.GPS_Sat[i].Psr[0] < 1)
					{
						//iono_correct = Klobuchar(&(Obs.ObsTime), &RecTemp, &(Calculation->GPS_POSnVEL[prn - 1].SatPos), azimuth, element, RawData);
						w[n_GPS] = Obs.GPS_Sat[i].Psr[1] - rho - delta_t_receive[0] + SpeedofLight*Calculation->GPS_POSnVEL[prn - 1].delta_tsv - tropo_correct;
					}
					else if (Obs.GPS_Sat[i].Psr[1] < 1)
					{
						//iono_correct = Klobuchar(&(Obs.ObsTime), &RecTemp, &(Calculation->GPS_POSnVEL[prn - 1].SatPos), azimuth, element, RawData);
						w[n_GPS] = Obs.GPS_Sat[i].Psr[0] - rho - delta_t_receive[0] + SpeedofLight*Calculation->GPS_POSnVEL[prn - 1].delta_tsv - tropo_correct;
					}
					else
					{
						//iono_correct = Klobuchar(&(Calculation->GPS_POSnVEL[prn - 1].SignalTrTime), &RecTemp, &(Calculation->GPS_POSnVEL[prn - 1].SatPos), azimuth, element, RawData);
						w[n_GPS] = (Obs.GPS_Sat[i].Psr[0] + Obs.GPS_Sat[i].Psr[1]) / 2 - rho - delta_t_receive[0] + SpeedofLight*Calculation->GPS_POSnVEL[prn - 1].delta_tsv - tropo_correct;
					}
					n_GPS++;
				}
			}
		}
		
		n_BDS = 0;		
		for (i = 0; i < Obs.BDS_SatNum; i++) //BDS
		{
			prn = Obs.BDS_Sat[i].Prn;
			tao = Obs.ObsTime.SecofWeek - Calculation->BDS_POSnVEL[prn - 1].SignalTrTime.SecofWeek;
			if (Calculation->BDS_POSnVEL[prn - 1].iono_flag == 1 || Calculation->BDS_POSnVEL[prn - 1].iono_flag == 3)//双频伪距之差过大
				continue;
			//if (Calculation->BDS_POSnVEL[prn - 1].iono_flag == 1)//双频伪距之差过大
			//	continue;
			if (BDS_SatPositionVelocity(Calculation->BDS_POSnVEL[prn - 1].SignalTrTime, prn, Calculation, map_eph, tao) == 2) //计算卫星位置速度
				continue;
			ComputeAzimuthElement(CGCS2000_a, CGCS2000_f, &RecTemp, &Calculation->BDS_POSnVEL[prn - 1].SatPos, &azimuth, &element); //计算高度角方位角（更新）

			if (abs(SPP_Cal->RecPos.x) < 1e-10 || element > 15 / 180.0*PI)
			{
				/////////
				if (abs(SPP_Cal->RecPos.x) < 1e-10 || azimuth > 0)
				{
					tropo_correct = Hopfield(&RecTemp, element); //计算对流层改正（更新）
					Calculation->BDS_POSnVEL[prn - 1].TroCorrect = tropo_correct;
					Calculation->BDS_POSnVEL[prn - 1].Element = element / PI*180.0;
					Calculation->BDS_POSnVEL[prn - 1].sigma2 = 0.04*0.04 + 0.03*0.03 / sin(element) / sin(element);
					//Calculation->BDS_POSnVEL[prn - 1].sigma2 = 0.003*0.003*(1 + 1.5*cos(element)*cos(element));

					rho = sqrt(pow((Calculation->BDS_POSnVEL[prn - 1].SatPos.x - RecTemp.x), 2) + pow((Calculation->BDS_POSnVEL[prn - 1].SatPos.y - RecTemp.y), 2) + pow((Calculation->BDS_POSnVEL[prn - 1].SatPos.z - RecTemp.z), 2));

					B[(n_GPS + n_BDS) * 5] = (RecTemp.x - Calculation->BDS_POSnVEL[prn - 1].SatPos.x) / rho;
					B[(n_GPS + n_BDS) * 5 + 1] = (RecTemp.y - Calculation->BDS_POSnVEL[prn - 1].SatPos.y) / rho;
					B[(n_GPS + n_BDS) * 5 + 2] = (RecTemp.z - Calculation->BDS_POSnVEL[prn - 1].SatPos.z) / rho;
					B[(n_GPS + n_BDS) * 5 + 3] = 0;
					B[(n_GPS + n_BDS) * 5 + 4] = 1;
					//P_value[n_GPS + n_BDS] = sin(element)*sin(element);
					//P_value[n_GPS + n_BDS] = 1 / (0.03*0.03 + 0.04*0.04 / sin(element) / sin(element));
					P_value[n_GPS + n_BDS] = 1 / (pow(0.13 + 0.56*exp(-element / PI * 180 / 10), 2));

					if (Calculation->BDS_POSnVEL[prn - 1].iono_flag == 2)
						w[n_GPS + n_BDS] = Calculation->BDS_POSnVEL[prn - 1].Psr_InonCorrect - rho - delta_t_receive[1] + SpeedofLight*Calculation->BDS_POSnVEL[prn - 1].delta_tsv - tropo_correct;
					else if (Obs.BDS_Sat[i].Psr[0] < 1)
					{
						//iono_correct = Klobuchar(&(Calculation->BDS_POSnVEL[prn - 1].SignalTrTime), &RecTemp, &(Calculation->BDS_POSnVEL[prn - 1].SatPos), azimuth, element, RawData);
						w[n_GPS + n_BDS] = Obs.BDS_Sat[i].Psr[1] - rho - delta_t_receive[1] + SpeedofLight*Calculation->BDS_POSnVEL[prn - 1].delta_tsv - tropo_correct;
					}
					else if (Obs.BDS_Sat[i].Psr[1] < 1)
					{
						//iono_correct = Klobuchar(&(Calculation->BDS_POSnVEL[prn - 1].SignalTrTime), &RecTemp, &(Calculation->BDS_POSnVEL[prn - 1].SatPos), azimuth, element, RawData);
						w[n_GPS + n_BDS] = Obs.BDS_Sat[i].Psr[0] - rho - delta_t_receive[1] + SpeedofLight*Calculation->BDS_POSnVEL[prn - 1].delta_tsv - tropo_correct;
					}
					else
					{
						//iono_correct = Klobuchar(&(Calculation->BDS_POSnVEL[prn - 1].SignalTrTime), &RecTemp, &(Calculation->BDS_POSnVEL[prn - 1].SatPos), azimuth, element, RawData);
						w[n_GPS + n_BDS] = (Obs.BDS_Sat[i].Psr[0] + Obs.BDS_Sat[i].Psr[1]) / 2 - rho - delta_t_receive[1] + SpeedofLight*Calculation->BDS_POSnVEL[prn - 1].delta_tsv - tropo_correct;
					}
					n_BDS++;
				}
			}
		}
		
		n = n_GPS + n_BDS; //可用于计算的全部卫星数量
		if (n_GPS != 0 && n_BDS != 0 && n > 5) //GPS和北斗都有
		{
			for (k = 0; k < n; k++)
			{
				//P[k*n + k] = P_value[k];
				P[k*n + k] = 1;
			}

			MatrixTrans(n, 5, B, B_T);
			MatrixMul(5, n, n, n, B_T, P, BTP);//
			MatrixMul(5, n, n, 5, BTP, B, BB_T);//
												//MatrixMul(5, n, n, 5, B_T, B, BB_T);
			MatrixInv(5, BB_T, BB_T1);
			MatrixMul(5, 5, 5, n, BB_T1, B_T, BB_T2);
			MatrixMul(5, n, n, n, BB_T2, P, BB_T3);//
												   //MatrixMul(5, n, n, 1, BB_T2, w, x);
			MatrixMul(5, n, n, 1, BB_T3, w, x);

			RecTemp.x = RecTemp.x + x[0];
			RecTemp.y = RecTemp.y + x[1];
			RecTemp.z = RecTemp.z + x[2];
			delta_t_receive[0] = delta_t_receive[0] + x[3];
			delta_t_receive[1] = delta_t_receive[1] + x[4];
			position_new = sqrt(pow((RecTemp.x), 2) + pow((RecTemp.y), 2) + pow((RecTemp.z), 2) + pow(delta_t_receive[0], 2) + pow(delta_t_receive[1], 2));

			if (fabs(position - position_new) > 1e-5)
			{
				position = position_new;
				continue;
			}
			else
			{
				MatrixMul(n, 5, 5, 1, B, x, v);
				MatrixSub(n, 1, n, 1, v, w, v);
				MatrixTrans(n, 1, v, v_T);
				MatrixMul(1, n, n, 1, v_T, v, vtpv);
				sigma = sqrt(vtpv[0] / (n - 5));
				PDOP = sqrt(BB_T1[0] * BB_T1[0] + BB_T1[6] * BB_T1[6] + BB_T1[12] * BB_T1[12]);

				SPP_Cal->RecPos.x = RecTemp.x;
				SPP_Cal->RecPos.y = RecTemp.y;
				SPP_Cal->RecPos.z = RecTemp.z;
				SPP_Cal->delta_t_receive[0] = delta_t_receive[0];
				SPP_Cal->delta_t_receive[1] = delta_t_receive[1];
				SPP_Cal->delta_t_dot[0] = delta_t_dot[0];
				SPP_Cal->delta_t_dot[1] = delta_t_dot[1];
				SPP_Cal->sigma = sigma;
				SPP_Cal->PDOP = PDOP;
				SPP_Cal->SatNum = n;
				SPP_Cal->time = Obs.ObsTime;
				XYZToBLH(CGCS2000_a, CGCS2000_f, &SPP_Cal->RecBlh, &SPP_Cal->RecPos);

				return 1;
			}
		}
		else if (n_BDS == 0 && n_GPS > 4) //单GPS
		{
			DeleteZero(B, 150);

			MatrixTrans(n, 4, B, B_T);
			MatrixMul(4, n, n, 4, B_T, B, BB_T);
			MatrixInv(4, BB_T, BB_T1);
			MatrixMul(4, 4, 4, n, BB_T1, B_T, BB_T2);
			MatrixMul(4, n, n, 1, BB_T2, w, x);

			RecTemp.x = RecTemp.x + x[0];
			RecTemp.y = RecTemp.y + x[1];
			RecTemp.z = RecTemp.z + x[2];
			delta_t_receive[0] = delta_t_receive[0] + x[3];
			position_new = sqrt(pow((RecTemp.x), 2) + pow((RecTemp.y), 2) + pow((RecTemp.z), 2));

			if (fabs(position - position_new) > 1)
			{
				position = position_new;
				continue;
			}
			else
			{
				MatrixMul(n, 4, 4, 1, B, x, v);
				MatrixSub(n, 1, n, 1, v, w, v);
				MatrixTrans(n, 1, v, v_T);
				MatrixMul(1, n, n, 1, v_T, v, vtpv);
				sigma = sqrt(vtpv[0] / (n - 4));

				SPP_Cal->RecPos.x = RecTemp.x;
				SPP_Cal->RecPos.y = RecTemp.y;
				SPP_Cal->RecPos.z = RecTemp.z;
				SPP_Cal->delta_t_receive[0] = delta_t_receive[0];
				SPP_Cal->delta_t_dot[0] = delta_t_dot[0];
				XYZToBLH(CGCS2000_a, CGCS2000_f, &SPP_Cal->RecBlh, &SPP_Cal->RecPos);

				return 1;
			}

		}
		else if (n_GPS == 0 && n_BDS > 4) //单北斗
		{
			DeleteZero(B, 150);

			MatrixTrans(n, 4, B, B_T);
			MatrixMul(4, n, n, 4, B_T, B, BB_T);
			MatrixInv(4, BB_T, BB_T1);
			MatrixMul(4, 4, 4, n, BB_T1, B_T, BB_T2);
			MatrixMul(4, n, n, 1, BB_T2, w, x);

			RecTemp.x = RecTemp.x + x[0];
			RecTemp.y = RecTemp.y + x[1];
			RecTemp.z = RecTemp.z + x[2];
			delta_t_receive[1] = delta_t_receive[1] + x[3];
			position_new = sqrt(pow((RecTemp.x), 2) + pow((RecTemp.y), 2) + pow((RecTemp.z), 2));

			if (fabs(position - position_new) > 1)
			{
				position = position_new;
				continue;
			}
			else
			{
				MatrixMul(n, 4, 4, 1, B, x, v);
				MatrixSub(n, 1, n, 1, v, w, v);
				MatrixTrans(n, 1, v, v_T);
				MatrixMul(1, n, n, 1, v_T, v, vtpv);
				sigma = sqrt(vtpv[0] / (n - 4));

				SPP_Cal->RecPos.x = RecTemp.x;
				SPP_Cal->RecPos.y = RecTemp.y;
				SPP_Cal->RecPos.z = RecTemp.z;
				SPP_Cal->delta_t_receive[1] = delta_t_receive[1];
				SPP_Cal->delta_t_dot[1] = delta_t_dot[1];
				XYZToBLH(CGCS2000_a, CGCS2000_f, &SPP_Cal->RecBlh, &SPP_Cal->RecPos);

				return 1;
			}
		}
		else
		{
			cout << "number of Sat too small" << endl;
			break;
		}
	}
	if (j == 100)
	{
		cout << "SPP failed!" << endl;
		return 0;
	}
}

/**************************************************************
GPS_SatPositionVelocity
目的：计算GPS卫星位置和速度

参数：
t             时刻（如信号发射时刻）
prn           卫星prn号
Calculation   GPS卫星位置、速度的计算结果
map_eph       星历
tao           信号传播时间

返回值：1=正常 2=星历过期 0=迭代失败
***************************************************************/
int GPS_SatPositionVelocity(GPSTIME t, int prn, CALCULATION *Calculation, map<int, map<GPSTIME, EPHEM>> map_eph, double &tao)
{
	double n0, t_k, n, M_k, e, E_k, delta, E_k_new, v_k, phi_k, delta_uk, delta_rk, delta_ik, uk, rk, ik, omega_k, xk, yk;
	double Ek_dot, phik_dot, uk_dot, rk_dot, ik_dot, omegak_dot, xk_dot, yk_dot, delta_t;
	int j;

	//离当前时刻最近的星历
	map<GPSTIME, EPHEM> ::iterator find_t;
	find_t = map_eph[prn].upper_bound(t);
	GPSTIME t_eph = find_t->first;

	//星历是否过期
	//delta_t = (t.Week - map_eph[prn][t_eph].week) * 604800 + (t.SecofWeek - map_eph[prn][t_eph].toe);
	delta_t = (t.Week - map_eph[prn][t_eph].week) * 604800;
	delta_t = (t.SecofWeek - map_eph[prn][t_eph].toe);
	delta_t = (t.Week - map_eph[prn][t_eph].week) * 604800 + (t.SecofWeek - map_eph[prn][t_eph].toe);
	if (delta_t <= 2 * 3600)
	{
		n0 = sqrt(GPS_miu / pow(pow(map_eph[prn][t_eph].sqrt_A, 2), 3));
		t_k = t.SecofWeek - map_eph[prn][t_eph].toe;
		if (t_k > 302400) t_k = t_k - 604800;
		else if (t_k < -302400) t_k = t_k + 604800;
		n = n0 + map_eph[prn][t_eph].deltaN;
		M_k = map_eph[prn][t_eph].M0 + n*t_k;
		E_k = M_k;
		delta = 1e-14;
		e = map_eph[prn][t_eph].ecc;
		for (j = 0; j < 10; j++)
		{
			E_k_new = M_k + sin(E_k)*e;
			if (abs(E_k_new - E_k) < delta) break;
			else E_k = E_k_new;
			if (j == 10) return 0;
		}
		v_k = atan2(sqrt(1 - e*e)*sin(E_k) / (1 - e*cos(E_k)), (cos(E_k) - e) / (1 - e*cos(E_k)));
		phi_k = v_k + map_eph[prn][t_eph].omega;
		delta_uk = map_eph[prn][t_eph].cus*sin(2 * phi_k) + map_eph[prn][t_eph].cuc*cos(2 * phi_k);
		delta_rk = map_eph[prn][t_eph].crs*sin(2 * phi_k) + map_eph[prn][t_eph].crc*cos(2 * phi_k);
		delta_ik = map_eph[prn][t_eph].cis*sin(2 * phi_k) + map_eph[prn][t_eph].cic*cos(2 * phi_k);
		uk = phi_k + delta_uk;
		rk = pow(map_eph[prn][t_eph].sqrt_A, 2)*(1 - e*cos(E_k)) + delta_rk;
		ik = map_eph[prn][t_eph].i0 + delta_ik + map_eph[prn][t_eph].IDOT*t_k;
		xk = rk*cos(uk);
		yk = rk*sin(uk);
		omega_k = map_eph[prn][t_eph].omega_0 + (map_eph[prn][t_eph].omega_dot - GPS_omegae_dot)*t_k - GPS_omegae_dot*map_eph[prn][t_eph].toe;
		Calculation->GPS_POSnVEL[prn - 1].SatPos.x = xk*cos(omega_k) - yk*cos(ik)*sin(omega_k);
		Calculation->GPS_POSnVEL[prn - 1].SatPos.y = xk*sin(omega_k) + yk*cos(ik)*cos(omega_k);
		Calculation->GPS_POSnVEL[prn - 1].SatPos.z = yk*sin(ik);

		EarthRotate(&Calculation->GPS_POSnVEL[prn - 1].SatPos, tao);

		Ek_dot = n / (1 - e*cos(E_k));
		phik_dot = sqrt((1 + e) / (1 - e))*pow(cos(v_k / 2) / cos(E_k / 2), 2)*Ek_dot;
		uk_dot = 2 * (map_eph[prn][t_eph].cus*cos(2 * phi_k) - map_eph[prn][t_eph].cuc*sin(2 * phi_k))*phik_dot + phik_dot;
		rk_dot = pow(map_eph[prn][t_eph].sqrt_A, 2)*e*sin(E_k)*Ek_dot + 2 * (map_eph[prn][t_eph].crs*cos(2 * phi_k) - map_eph[prn][t_eph].crc*sin(2 * phi_k))*phik_dot;
		ik_dot = map_eph[prn][t_eph].IDOT + 2 * (map_eph[prn][t_eph].cis*cos(2 * phi_k) - map_eph[prn][t_eph].cic*sin(2 * phi_k))*phik_dot;
		omegak_dot = map_eph[prn][t_eph].omega_dot - GPS_omegae_dot;
		xk_dot = rk_dot*cos(uk) - rk*uk_dot*sin(uk);
		yk_dot = rk_dot*sin(uk) + rk*uk_dot*cos(uk);
		Calculation->GPS_POSnVEL[prn - 1].v_x = cos(omega_k)*xk_dot - sin(omega_k)*cos(ik)*yk_dot - (xk*sin(omega_k) + yk*cos(omega_k)*cos(ik))*omegak_dot + yk*sin(omega_k)*sin(ik)*ik_dot;
		Calculation->GPS_POSnVEL[prn - 1].v_y = sin(omega_k)*xk_dot + cos(omega_k)*cos(ik)*yk_dot + (xk*cos(omega_k) - yk*sin(omega_k)*cos(ik))*omegak_dot + yk*cos(omega_k)*sin(ik)*ik_dot;
		Calculation->GPS_POSnVEL[prn - 1].v_z = sin(ik)*yk_dot + yk*cos(ik)*ik_dot;
		return 1;
	}
	else
		return 2;
}

/**************************************************************
GPS_SatClock
目的：计算GPS卫星钟差与钟速

参数：
t_sv          时刻（如信号发射时刻 周内秒）
prn           卫星prn号
Calculation   卫星位置、速度的计算结果
map_eph       观测数据、星历等数据

返回值：1=正常 0=迭代失败
***************************************************************/
int GPS_SatClock(GPSTIME t_sv, int prn, CALCULATION *Calculation, map<int, map<GPSTIME, EPHEM>> map_eph)
{
	double F, delta_tr, e, A, E_k, t_k, n0, n, M_k, E_k_new, t, delta_t, delta_t_new;
	double clkbias, clkdrift, clkdriftrate, Ek_dot, delta_tr_dot;
	int i, j;

	//离当前时刻最近的星历
	map<GPSTIME, EPHEM> ::iterator find_t;
	find_t = map_eph[prn].upper_bound(t_sv);
	GPSTIME t_eph = find_t->first;

	F = -4.442807633e-10;
	e = map_eph[prn][t_eph].ecc;
	A = map_eph[prn][t_eph].sqrt_A*map_eph[prn][t_eph].sqrt_A;
	clkbias = map_eph[prn][t_eph].af0;
	clkdrift = map_eph[prn][t_eph].af1;
	clkdriftrate = map_eph[prn][t_eph].af2;
	delta_t = 0;
	for (i = 0; i < 10; i++)
	{
		t = t_sv.SecofWeek - delta_t;
		/*计算E_k*/
		n0 = sqrt(GPS_miu / pow(A, 3));
		t_k = t_sv.SecofWeek - map_eph[prn][t_eph].toe;
		if (t_k > 302400) t_k = t_k - 604800;
		else if (t_k < -302400) t_k = t_k + 604800;
		n = n0 + map_eph[prn][t_eph].deltaN;
		M_k = map_eph[prn][t_eph].M0 + n*t_k;
		E_k = M_k;
		e = map_eph[prn][t_eph].ecc;
		for (j = 0; j < 10; j++)
		{
			E_k_new = M_k + sin(E_k)*e;
			if (abs(E_k_new - E_k) < 1e-14) break;
			else E_k = E_k_new;
			if (j == 9) return 0;
		}
		delta_tr = F*e*sqrt(A)*sin(E_k);
		delta_t_new = clkbias + clkdrift*(t - t_eph.SecofWeek)
			+ clkdriftrate*pow(t - t_eph.SecofWeek, 2) + delta_tr;
		//-RawData->GPSEph[prn - 1].tgd[0];
		if (abs(delta_t - delta_t_new) < 1e-10) break;
		else delta_t = delta_t_new;
		if (i == 9) return 0;
	}
	Calculation->GPS_POSnVEL[prn - 1].delta_tsv = delta_t;

	Ek_dot = n / (1 - e*cos(E_k));
	delta_tr_dot = F*e*sqrt(A)*cos(E_k)*Ek_dot;
	Calculation->GPS_POSnVEL[prn - 1].delta_tsv_dot = clkdrift + 2 * clkdriftrate*(t - delta_t - t_eph.SecofWeek) + delta_tr_dot;
	return 1;
}

/**************************************************************
BDS_SatPositionVelocity
目的：计算BDS卫星位置和速度

参数：
t             时刻（如信号发射时刻）
prn           卫星prn号
Calculation   卫星位置、速度、钟差、钟速的计算结果
RawData       观测数据、星历等数据
tao           信号传播时间

返回值：1=正常 0=迭代失败
***************************************************************/
int BDS_SatPositionVelocity(GPSTIME t, int prn, CALCULATION *Calculation, map<int, map<GPSTIME, EPHEM>> map_eph, double &tao)
{
	double n0, t_k, n, M_k, e, E_k, delta, E_k_new, v_k, phi_k, delta_uk, delta_rk, delta_ik, uk, rk, ik, omega_k, xk, yk, x_GK, y_GK, z_GK;
	double Ek_dot, phik_dot, uk_dot, rk_dot, ik_dot, omegak_dot, xk_dot, yk_dot, delta_t, vx_GK, vy_GK, vz_GK, rotate_1, rotate_2;
	int j;
	//星历是否过期
	//BDST = GPST - 1356week - 14s;
	//delta_t = (t->Week - RawData->BDSEph[prn - 1].week) * 604800 + (t->SecofWeek - RawData->BDSEph[prn - 1].toe);

	//离当前时刻最近的星历
	map<GPSTIME, EPHEM> ::iterator find_t;
	find_t = map_eph[prn + 32].upper_bound(t);
	GPSTIME t_eph = (find_t)->first;

	//delta_t = (t.Week - map_eph[prn + 32][t_eph].week) * 604800 + (t.SecofWeek - map_eph[prn + 32][t_eph].toe);
	delta_t = (t.Week - map_eph[prn + 32][t_eph].week) * 604800 + (t.SecofWeek - map_eph[prn + 32][t_eph].toe - 14);
	if (delta_t <= 3600)
	{
		n0 = sqrt(BDS_miu / pow(pow(map_eph[prn + 32][t_eph].sqrt_A, 2), 3));
		//t_k = t.SecofWeek - map_eph[prn + 32][t_eph].toe;
		t_k = t.SecofWeek - map_eph[prn + 32][t_eph].toe - 14;
		if (t_k > 302400) t_k = t_k - 604800;
		else if (t_k < -302400) t_k = t_k + 604800;
		n = n0 + map_eph[prn + 32][t_eph].deltaN;
		M_k = map_eph[prn + 32][t_eph].M0 + n*t_k;
		E_k = M_k;
		delta = 1e-14;
		e = map_eph[prn + 32][t_eph].ecc;
		for (j = 0; j < 10; j++)
		{
			E_k_new = M_k + sin(E_k)*e;
			if (abs(E_k_new - E_k) < delta) break;
			else E_k = E_k_new;
			if (j == 10) return 0;
		}
		v_k = atan2(sqrt(1 - e*e)*sin(E_k) / (1 - e*cos(E_k)), (cos(E_k) - e) / (1 - e*cos(E_k)));
		phi_k = v_k + map_eph[prn + 32][t_eph].omega;
		delta_uk = map_eph[prn + 32][t_eph].cus*sin(2 * phi_k) + map_eph[prn + 32][t_eph].cuc*cos(2 * phi_k);
		delta_rk = map_eph[prn + 32][t_eph].crs*sin(2 * phi_k) + map_eph[prn + 32][t_eph].crc*cos(2 * phi_k);
		delta_ik = map_eph[prn + 32][t_eph].cis*sin(2 * phi_k) + map_eph[prn + 32][t_eph].cic*cos(2 * phi_k);
		uk = phi_k + delta_uk;
		rk = pow(map_eph[prn + 32][t_eph].sqrt_A, 2)*(1 - e*cos(E_k)) + delta_rk;
		ik = map_eph[prn + 32][t_eph].i0 + delta_ik + map_eph[prn + 32][t_eph].IDOT*t_k;
		xk = rk*cos(uk);
		yk = rk*sin(uk);

		Ek_dot = n / (1 - e*cos(E_k));
		phik_dot = sqrt((1 + e) / (1 - e))*pow(cos(v_k / 2) / cos(E_k / 2), 2)*Ek_dot;
		uk_dot = 2 * (map_eph[prn + 32][t_eph].cus*cos(2 * phi_k) - map_eph[prn + 32][t_eph].cuc*sin(2 * phi_k))*phik_dot + phik_dot;
		rk_dot = pow(map_eph[prn + 32][t_eph].sqrt_A, 2)*e*sin(E_k)*Ek_dot + 2 * (map_eph[prn + 32][t_eph].crs*cos(2 * phi_k) - map_eph[prn + 32][t_eph].crc*sin(2 * phi_k))*phik_dot;
		ik_dot = map_eph[prn + 32][t_eph].IDOT + 2 * (map_eph[prn + 32][t_eph].cis*cos(2 * phi_k) - map_eph[prn + 32][t_eph].cic*sin(2 * phi_k))*phik_dot;
		xk_dot = rk_dot*cos(uk) - rk*uk_dot*sin(uk);
		yk_dot = rk_dot*sin(uk) + rk*uk_dot*cos(uk);

		if (prn == 1|| prn == 2|| prn == 3|| prn == 4|| prn == 5|| prn == 59|| prn == 60|| prn == 61) //GEO
		{
			omega_k = map_eph[prn + 32][t_eph].omega_0 + map_eph[prn + 32][t_eph].omega_dot*t_k - BDS_omegae_dot*map_eph[prn + 32][t_eph].toe;
			omegak_dot = map_eph[prn + 32][t_eph].omega_dot;
			x_GK = xk*cos(omega_k) - yk*cos(ik)*sin(omega_k);
			y_GK = xk*sin(omega_k) + yk*cos(ik)*cos(omega_k);
			z_GK = yk*sin(ik);
			rotate_1 = BDS_omegae_dot*t_k;
			rotate_2 = -5.0 / 180.0 * PI;
			double Rz[9] = { cos(rotate_1), sin(rotate_1), 0, -sin(rotate_1), cos(rotate_1), 0, 0, 0, 1 };
			double Rx[9] = { 1,0,0,0,cos(rotate_2),sin(rotate_2),0,-sin(rotate_2),cos(rotate_2) };
			double Rzx[9];
			MatrixMul(3, 3, 3, 3, Rz, Rx, Rzx);
			Calculation->BDS_POSnVEL[prn - 1].SatPos.x = Rzx[0] * x_GK + Rzx[1] * y_GK + Rzx[2] * z_GK;
			Calculation->BDS_POSnVEL[prn - 1].SatPos.y = Rzx[3] * x_GK + Rzx[4] * y_GK + Rzx[5] * z_GK;
			Calculation->BDS_POSnVEL[prn - 1].SatPos.z = Rzx[6] * x_GK + Rzx[7] * y_GK + Rzx[8] * z_GK;

			vx_GK = cos(omega_k)*xk_dot - sin(omega_k)*cos(ik)*yk_dot - (xk*sin(omega_k) + yk*cos(omega_k)*cos(ik))*omegak_dot + yk*sin(omega_k)*sin(ik)*ik_dot;
			vy_GK = sin(omega_k)*xk_dot + cos(omega_k)*cos(ik)*yk_dot + (xk*cos(omega_k) - yk*sin(omega_k)*cos(ik))*omegak_dot + yk*cos(omega_k)*sin(ik)*ik_dot;
			vz_GK = sin(ik)*yk_dot + yk*cos(ik)*ik_dot;

			Calculation->BDS_POSnVEL[prn - 1].v_x = BDS_omegae_dot*(-sin(rotate_1))*x_GK + cos(rotate_1)*vx_GK
				+ BDS_omegae_dot*cos(rotate_1)*cos(rotate_2)*y_GK + sin(rotate_1)*cos(rotate_2)*vy_GK
				+ BDS_omegae_dot*cos(rotate_1)*sin(rotate_2)*z_GK + sin(rotate_1)*sin(rotate_2)*vz_GK;
			Calculation->BDS_POSnVEL[prn - 1].v_y = -BDS_omegae_dot*cos(rotate_1)*x_GK - sin(rotate_1)*vx_GK
				- BDS_omegae_dot*sin(rotate_1)*cos(rotate_2)*y_GK + cos(rotate_1)*cos(rotate_2)*vy_GK
				- BDS_omegae_dot*sin(rotate_1)*sin(rotate_2)*z_GK + cos(rotate_1)*sin(rotate_2)*vz_GK;
			Calculation->BDS_POSnVEL[prn - 1].v_z = -sin(rotate_2)*vy_GK + cos(rotate_2)*vz_GK;
		}
		else
		{
			omega_k = map_eph[prn + 32][t_eph].omega_0 + (map_eph[prn + 32][t_eph].omega_dot - BDS_omegae_dot)*t_k - BDS_omegae_dot*map_eph[prn + 32][t_eph].toe;
			omegak_dot = map_eph[prn + 32][t_eph].omega_dot - BDS_omegae_dot;
			Calculation->BDS_POSnVEL[prn - 1].SatPos.x = xk*cos(omega_k) - yk*cos(ik)*sin(omega_k);
			Calculation->BDS_POSnVEL[prn - 1].SatPos.y = xk*sin(omega_k) + yk*cos(ik)*cos(omega_k);
			Calculation->BDS_POSnVEL[prn - 1].SatPos.z = yk*sin(ik);

			Calculation->BDS_POSnVEL[prn - 1].v_x = cos(omega_k)*xk_dot - sin(omega_k)*cos(ik)*yk_dot - (xk*sin(omega_k) + yk*cos(omega_k)*cos(ik))*omegak_dot + yk*sin(omega_k)*sin(ik)*ik_dot;
			Calculation->BDS_POSnVEL[prn - 1].v_y = sin(omega_k)*xk_dot + cos(omega_k)*cos(ik)*yk_dot + (xk*cos(omega_k) - yk*sin(omega_k)*cos(ik))*omegak_dot + yk*cos(omega_k)*sin(ik)*ik_dot;
			Calculation->BDS_POSnVEL[prn - 1].v_z = sin(ik)*yk_dot + yk*cos(ik)*ik_dot;
		}

		EarthRotate(&Calculation->BDS_POSnVEL[prn - 1].SatPos, tao);
		return 1;
	}
	else
		return 2;
}

/**************************************************************
BDS_SatClock
目的：计算BDS卫星钟差与钟速

参数：
t_sv          时刻（如信号发射时刻）
prn           卫星prn号
Calculation   卫星位置、速度的计算结果
RawData       观测数据、星历等数据

返回值：1=正常 0=迭代失败
***************************************************************/
int BDS_SatClock(GPSTIME t_sv, int prn, CALCULATION *Calculation, map<int, map<GPSTIME, EPHEM>> map_eph)
{
	double F, delta_tr, e, A, E_k, t_k, n0, n, M_k, E_k_new, t, delta_t, delta_t_new;
	double clkbias, clkdrift, clkdriftrate, Ek_dot, delta_tr_dot;
	int i, j;

	prn = prn + 32;
	//离当前时刻最近的星历
	map<GPSTIME, EPHEM> ::iterator find_t;
	GPSTIME t_eph;
	find_t = map_eph[prn].upper_bound(t_sv);
	t_eph = find_t->first;


	F = -4.442807633e-10;
	e = map_eph[prn][t_eph].ecc;
	A = map_eph[prn][t_eph].sqrt_A*map_eph[prn][t_eph].sqrt_A;
	clkbias = map_eph[prn][t_eph].af0;
	clkdrift = map_eph[prn][t_eph].af1;
	clkdriftrate = map_eph[prn][t_eph].af2;
	delta_t = 0;
	for (i = 0; i < 10; i++)
	{
		//t = t_sv.SecofWeek - delta_t - 14;
		t = t_sv.SecofWeek - delta_t;
		/*计算E_k*/
		n0 = sqrt(BDS_miu / pow(A, 3));
		//t_k = t_sv.SecofWeek - map_eph[prn][t_eph].toe - 14;
		t_k = t_sv.SecofWeek - map_eph[prn][t_eph].toe;
		if (t_k > 302400) t_k = t_k - 604800;
		else if (t_k < -302400) t_k = t_k + 604800;
		n = n0 + map_eph[prn][t_eph].deltaN;
		M_k = map_eph[prn][t_eph].M0 + n*t_k;
		E_k = M_k;
		e = map_eph[prn][t_eph].ecc;
		for (j = 0; j < 10; j++)
		{
			E_k_new = M_k + sin(E_k)*e;
			if (abs(E_k_new - E_k) < 1e-10) break;
			else E_k = E_k_new;
			if (j == 9) return 0;
		}
		delta_tr = F*e*sqrt(A)*sin(E_k);
		delta_t_new = clkbias + clkdrift*(t - t_eph.SecofWeek)
			+ clkdriftrate*pow(t - t_eph.SecofWeek, 2) + delta_tr;
		if (abs(delta_t - delta_t_new) < 1e-10) break;
		else delta_t = delta_t_new;
		if (i == 9) return 0;
	}
	Calculation->BDS_POSnVEL[prn - 32 - 1].delta_tsv = delta_t;

	Ek_dot = n / (1 - e*cos(E_k));
	delta_tr_dot = F*e*sqrt(A)*cos(E_k)*Ek_dot;
	Calculation->BDS_POSnVEL[prn - 32 - 1].delta_tsv_dot = clkdrift + 2 * clkdriftrate*(t - delta_t - t_eph.SecofWeek) + delta_tr_dot;
	return 1;
}

/**************************************************************
GPS_SignalTransmit
目的：计算GPS信号发射时间

参数：
RawData        原始数据
Calculation   计算结果

返回值：0=计算失败 1=成功
***************************************************************/
int GPS_SignalTransmit(OBS Obs, CALCULATION *Calculation, map<int, map<GPSTIME, EPHEM>> map_eph)
{
	double t_tr, P, delta_t;
	int i, j, prn;

	for (i = 0; i < Obs.GPS_SatNum; i++)
	{
		prn = Obs.GPS_Sat[i].Prn;
		delta_t = 0;
		//计算伪距（单频双频）
		if (Obs.GPS_Sat[i].Psr[0]>1 && Obs.GPS_Sat[i].Psr[1] > 1)
			P = (GPS_L1*GPS_L1* Obs.GPS_Sat[i].Psr[0] - GPS_L2*GPS_L2* Obs.GPS_Sat[i].Psr[1]) / (GPS_L1*GPS_L1 - GPS_L2*GPS_L2);
		else if (Obs.GPS_Sat[i].Psr[1] < 1)
		{
			P = Obs.GPS_Sat[i].Psr[0];
		}
		else if (Obs.GPS_Sat[i].Psr[0] < 1)
		{
			P = Obs.GPS_Sat[i].Psr[1];
		}
		else
			return 0;

		Calculation->GPS_POSnVEL[prn - 1].SignalTrTime.Week = Obs.ObsTime.Week;
		for (j = 0; j < 2; j++)
		{
			t_tr = Obs.ObsTime.SecofWeek - P / SpeedofLight - delta_t;
			Calculation->GPS_POSnVEL[prn - 1].SignalTrTime.SecofWeek = t_tr;
			GPS_SatClock(Calculation->GPS_POSnVEL[prn - 1].SignalTrTime, prn, Calculation, map_eph);
			delta_t = Calculation->GPS_POSnVEL[prn - 1].delta_tsv;
		}
		Calculation->GPS_POSnVEL[prn - 1].SignalTrTime.SecofWeek = t_tr;
	}
	return 1;
}

/**************************************************************
BDS_SignalTransmit
目的：计算BDS信号发射时间

参数：
RawData        原始数据
Calculation   计算结果

返回值：0=计算失败 1=成功
***************************************************************/
int BDS_SignalTransmit(OBS Obs, CALCULATION *Calculation, map<int, map<GPSTIME, EPHEM>> map_eph)
{
	double t_tr, P, delta_t;
	int i, j, prn;
	for (i = 0; i < Obs.BDS_SatNum; i++)
	{
		prn = Obs.BDS_Sat[i].Prn;
		delta_t = 0;

		P = (BDS_B1*BDS_B1* Obs.BDS_Sat[i].Psr[0] - BDS_B3*BDS_B3* Obs.BDS_Sat[i].Psr[1]) / (BDS_B1*BDS_B1 - BDS_B3*BDS_B3);
		/*
		//计算伪距（单频双频）
		if (Obs.BDS_Sat[i].Psr[0]>1 && Obs.BDS_Sat[i].Psr[1] > 1)
			P = (BDS_B1*BDS_B1* Obs.BDS_Sat[i].Psr[0] - BDS_B3*BDS_B3* Obs.BDS_Sat[i].Psr[1]) / (BDS_B1*BDS_B1 - BDS_B3*BDS_B3);
		else if (Obs.BDS_Sat[i].Psr[1] < 1)
		{
			P = Obs.BDS_Sat[i].Psr[0];
		}
		else if (Obs.BDS_Sat[i].Psr[0] < 1)
		{
			P = Obs.BDS_Sat[i].Psr[1];
		}
		else
			return 0;
		*/

		Calculation->BDS_POSnVEL[prn - 1].SignalTrTime.Week = Obs.ObsTime.Week;
		for (j = 0; j < 2; j++)
		{
			t_tr = Obs.ObsTime.SecofWeek - P / SpeedofLight - delta_t;
			Calculation->BDS_POSnVEL[prn - 1].SignalTrTime.SecofWeek = t_tr;
			BDS_SatClock(Calculation->BDS_POSnVEL[prn - 1].SignalTrTime, prn, Calculation, map_eph);
			delta_t = Calculation->BDS_POSnVEL[prn - 1].delta_tsv;
		}
		Calculation->BDS_POSnVEL[prn - 1].SignalTrTime.SecofWeek = t_tr;
	}
	return 1;
}



/**************************************************************
SPV
目的：单点测速

参数：
Obs                  原始观测数据
RecPos               接收机位置
RecVel               接收机速度
Calculation          卫星位置速度计算结果
sigma2_v             方差阵

返回值：0=卫星数目不足，无法计算  1=正常
***************************************************************/
int SPV(CALCULATION *Calculation, OBS Obs, XYZ *RecPos, XYZ *RecVel, XYZ *sigma2_v)
{
	int i, prn, n_GPS, n;
	double Sat_X, Sat_Y, Sat_Z, rho, rho_dot;
	double B[60], w[15], D[15], B_T[60], BTP[60], BTPB[16], N[16], NBT[60], NBTP[60], x[4];
	double Bx[15], v[15], vT[15], vTP[15], vTPv[1], sigma0;
	double Q[225], P[225];

	if (Obs.GPS_SatNum > 4)
	{
		n_GPS = 0;
		memset(Q, 0, sizeof(Q));
		for (i = 0; i < Obs.GPS_SatNum; i++)
		{
			prn = Obs.GPS_Sat[i].Prn;
			Sat_X = Calculation->GPS_POSnVEL[prn - 1].SatPos.x;
			Sat_Y = Calculation->GPS_POSnVEL[prn - 1].SatPos.y;
			Sat_Z = Calculation->GPS_POSnVEL[prn - 1].SatPos.z;

			rho = sqrt(pow((Sat_X - RecPos->x), 2) + pow((Sat_Y - RecPos->y), 2) + pow((Sat_Z - RecPos->z), 2));
			B[n_GPS * 4] = (RecPos->x - Calculation->GPS_POSnVEL[prn - 1].SatPos.x) / rho;
			B[n_GPS * 4 + 1] = (RecPos->y - Calculation->GPS_POSnVEL[prn - 1].SatPos.y) / rho;
			B[n_GPS * 4 + 2] = (RecPos->z - Calculation->GPS_POSnVEL[prn - 1].SatPos.z) / rho;
			B[n_GPS * 4 + 3] = 1;
			rho_dot = B[n_GPS * 4] * Calculation->GPS_POSnVEL[prn - 1].v_x + B[n_GPS * 4 + 1] * Calculation->GPS_POSnVEL[prn - 1].v_y + B[n_GPS * 4 + 2] * Calculation->GPS_POSnVEL[prn - 1].v_z;
			rho_dot = -rho_dot;
			D[n_GPS] = Obs.GPS_Sat[i].Dop[0] * SpeedofLight / GPS_L1;
			w[n_GPS] = -D[n_GPS] - rho_dot + SpeedofLight*Calculation->GPS_POSnVEL[prn - 1].delta_tsv_dot;
			Q[i*Obs.GPS_SatNum + i] = Calculation->GPS_POSnVEL[prn - 1].sigma2;

			n_GPS++;
		}

		n = n_GPS;
		MatrixInv(n, Q, P);

		MatrixTrans(n, 4, B, B_T);
		MatrixMul(4, n, n, n, B_T, P, BTP);
		MatrixMul(4, n, n, 4, BTP, B, BTPB);
		MatrixInv(4, BTPB, N);
		MatrixMul(4, 4, 4, n, N, B_T, NBT);
		MatrixMul(4, n, n, n, NBT, P, NBTP);
		MatrixMul(4, n, n, 1, NBTP, w, x);

		RecVel->x = x[0];
		RecVel->y = x[1];
		RecVel->z = x[2];

		//方差
		MatrixMul(n, 4, 4, 1, B, x, Bx);
		MatrixSub(n, 1, n, 1, Bx, w, v);
		MatrixTrans(n, 1, v, vT);
		MatrixMul(1, n, n, n, vT, P, vTP);
		MatrixMul(1, n, n, 1, vTP, v, vTPv);
		sigma0 = sqrt(vTPv[0] / (n - 4));
		sigma2_v->x = sigma0*N[0];
		sigma2_v->y = sigma0*N[5];
		sigma2_v->z = sigma0*N[10];

		return 1;
	}

	else
	{
		cout << "Compute Velocity: number of Sat too small" << endl;
		return 0;
	}
}