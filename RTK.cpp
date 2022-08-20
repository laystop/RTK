#include "stdafx.h"
#include "struct.h"
#include<iostream>
#include<cmath>
#include<string>
#include<fstream>

using namespace std;

/*Obs_1流动站，Obs_2基准站*/
int RTK(OBS Obs_1, OBS Obs_2, CALCULATION Calculation, RTK_CAL *RTK_Cal)
{
	//找两个基站的共视卫星、高度角最高的卫星
	int GPS_prn_common[SatComMax/2], BDS_prn_common[SatComMax/2];   //共视卫星prn
	int GPS_num_common = 0, BDS_num_common = 0;						//共视卫星个数
	int center_GPS_prn, center_BDS_prn;								//基准卫星prn
	double center_GPS_element = 0, center_BDS_element = 0;			//基准卫星高度角

	double Q_GPS[SatComMax*SatComMax/4], Q_BDS[SatComMax*SatComMax/4]; //协方差
	double Q[4 * SatComMax * 4 * SatComMax];
	double P[4 * SatComMax * 4 * SatComMax];
	double sigma2_P_phi = 90000; //伪距和相位方差之比

	//最小二乘
	double AT[4 * SatComMax * 2 * SatComMax], ATP[2 * SatComMax * 4 * SatComMax],
		ATPA[2 * SatComMax * 2 * SatComMax], N[2 * SatComMax * 2 * SatComMax],
		NAT[2 * SatComMax * 4 * SatComMax], NATP[2 * SatComMax * 4 * SatComMax],
		NATPL[2 * SatComMax], result[2 * SatComMax];
	double Q_float[2 * SatComMax * 2 * SatComMax];

	memset(result, 0, sizeof(result));
	RTK_Cal->time = Obs_1.ObsTime;

	//找基准卫星
	for (int i = 0; i < Obs_1.GPS_SatNum; i++)
	{
		for (int j = 0; j < Obs_2.GPS_SatNum; j++)
		{
			if (Obs_1.GPS_Sat[i].Prn == Obs_2.GPS_Sat[j].Prn&&Calculation.GPS_POSnVEL[Obs_1.GPS_Sat[i].Prn - 1].Element > 10 && Calculation.GPS_POSnVEL[Obs_2.GPS_Sat[j].Prn - 1].Element > 10)
			{
				GPS_prn_common[GPS_num_common] = Obs_1.GPS_Sat[i].Prn;
				GPS_num_common++;
				if (Calculation.GPS_POSnVEL[Obs_1.GPS_Sat[i].Prn - 1].Element > center_GPS_element)
				{
					center_GPS_prn = Obs_1.GPS_Sat[i].Prn;
					center_GPS_element = Calculation.GPS_POSnVEL[Obs_1.GPS_Sat[i].Prn - 1].Element;
				}
			}
		}
	}
	for (int i = 0; i < Obs_1.BDS_SatNum; i++)
	{
		for (int j = 0; j < Obs_2.BDS_SatNum; j++)
		{
			if (Obs_1.BDS_Sat[i].Prn == Obs_2.BDS_Sat[j].Prn&&Calculation.BDS_POSnVEL[Obs_1.BDS_Sat[i].Prn - 1].Element > 10 && Calculation.BDS_POSnVEL[Obs_2.BDS_Sat[j].Prn - 1].Element > 10)
			{
				BDS_prn_common[BDS_num_common] = Obs_1.BDS_Sat[i].Prn;
				BDS_num_common++;
				if (Calculation.BDS_POSnVEL[Obs_1.BDS_Sat[i].Prn - 1].Element > center_BDS_element)
				{
					center_BDS_prn = Obs_1.BDS_Sat[i].Prn;
					center_BDS_element = Calculation.BDS_POSnVEL[Obs_1.BDS_Sat[i].Prn - 1].Element;
				}
			}
		}
	}

	//最小二乘
	for (int times = 0; times < 20; times++)
	{
		//双频双系统函数模型建立
		//方程个数：GPS_num_common-1 + BDS_num_common-1
		double L[4 * SatComMax], A[4 * SatComMax * 2 * SatComMax];
		double pho_bk, pho_bj, pho_ak, pho_aj;
		memset(A, 0, sizeof(A));

		//GPS
		for (int i = 0; i < GPS_num_common - 1; i++)
		{
			if (GPS_prn_common[i] != center_GPS_prn)
			{
				pho_bk = sqrt(pow((RTK_Cal->StaPos_flow_tmp.x - Calculation.GPS_POSnVEL[GPS_prn_common[i] - 1].SatPos.x), 2)
					+ pow((RTK_Cal->StaPos_flow_tmp.y - Calculation.GPS_POSnVEL[GPS_prn_common[i] - 1].SatPos.y), 2)
					+ pow((RTK_Cal->StaPos_flow_tmp.z - Calculation.GPS_POSnVEL[GPS_prn_common[i] - 1].SatPos.z), 2));
				pho_bj = sqrt(pow((RTK_Cal->StaPos_flow_tmp.x - Calculation.GPS_POSnVEL[center_GPS_prn - 1].SatPos.x), 2)
					+ pow((RTK_Cal->StaPos_flow_tmp.y - Calculation.GPS_POSnVEL[center_GPS_prn - 1].SatPos.y), 2)
					+ pow((RTK_Cal->StaPos_flow_tmp.z - Calculation.GPS_POSnVEL[center_GPS_prn - 1].SatPos.z), 2));
				pho_ak = sqrt(pow((RTK_Cal->StaPos_ref.x - Calculation.GPS_POSnVEL[GPS_prn_common[i] - 1].SatPos.x), 2)
					+ pow((RTK_Cal->StaPos_ref.y - Calculation.GPS_POSnVEL[GPS_prn_common[i] - 1].SatPos.y), 2)
					+ pow((RTK_Cal->StaPos_ref.z - Calculation.GPS_POSnVEL[GPS_prn_common[i] - 1].SatPos.z), 2));
				pho_aj = sqrt(pow((RTK_Cal->StaPos_ref.x - Calculation.GPS_POSnVEL[center_GPS_prn - 1].SatPos.x), 2)
					+ pow((RTK_Cal->StaPos_ref.y - Calculation.GPS_POSnVEL[center_GPS_prn - 1].SatPos.y), 2)
					+ pow((RTK_Cal->StaPos_ref.z - Calculation.GPS_POSnVEL[center_GPS_prn - 1].SatPos.z), 2));
				int index_1, index_2;                   //GPS_prn_common[i]卫星在obs1和obs2中存储的位置
				int index_center_1, index_center_2;     //center_GPS_prn卫星在obs1和obs2中存储的位置
				for (int j = 0; j < Obs_1.GPS_SatNum; j++)
				{
					if (Obs_1.GPS_Sat[j].Prn == GPS_prn_common[i])
						index_1 = j;
					if (Obs_1.GPS_Sat[j].Prn == center_GPS_prn)
						index_center_1 = j;
				}
				for (int j = 0; j < Obs_2.GPS_SatNum; j++)
				{
					if (Obs_2.GPS_Sat[j].Prn == GPS_prn_common[i])
						index_2 = j;
					if (Obs_2.GPS_Sat[j].Prn == center_GPS_prn)
						index_center_2 = j;
				}

				L[i] = lambda_L1*(Obs_1.GPS_Sat[index_1].Adr[0] - Obs_2.GPS_Sat[index_2].Adr[0]
					- Obs_1.GPS_Sat[index_center_1].Adr[0] + Obs_2.GPS_Sat[index_center_2].Adr[0])
					- pho_bk + pho_bj + pho_ak - pho_aj;   //GPS L1 相位
				L[i + GPS_num_common - 1] = 
					lambda_L2*(Obs_1.GPS_Sat[index_1].Adr[1] - Obs_2.GPS_Sat[index_2].Adr[1]
					- Obs_1.GPS_Sat[index_center_1].Adr[1] + Obs_2.GPS_Sat[index_center_2].Adr[1])
					- pho_bk + pho_bj + pho_ak - pho_aj;   //GPS L2 相位
				L[i + 2 * (GPS_num_common - 1) + 2 * (BDS_num_common - 1)] =
					Obs_1.GPS_Sat[index_1].Psr[0] - Obs_2.GPS_Sat[index_2].Psr[0]
					- Obs_1.GPS_Sat[index_center_1].Psr[0] + Obs_2.GPS_Sat[index_center_2].Psr[0]
					- pho_bk + pho_bj + pho_ak - pho_aj;   //GPS L1 伪距
				L[i + 3 * (GPS_num_common - 1) + 2 * (BDS_num_common - 1)] =
					Obs_1.GPS_Sat[index_1].Psr[1] - Obs_2.GPS_Sat[index_2].Psr[1]
					- Obs_1.GPS_Sat[index_center_1].Psr[1] + Obs_2.GPS_Sat[index_center_2].Psr[1]
					- pho_bk + pho_bj + pho_ak - pho_aj;   //GPS L2 伪距

				////////
				double A1, A2, A3;
				A1 = -(Calculation.GPS_POSnVEL[GPS_prn_common[i] - 1].SatPos.x - RTK_Cal->StaPos_flow_tmp.x) / pho_bk
					+ (Calculation.GPS_POSnVEL[center_GPS_prn - 1].SatPos.x - RTK_Cal->StaPos_flow_tmp.x) / pho_bj;
				A2 = -(Calculation.GPS_POSnVEL[GPS_prn_common[i] - 1].SatPos.y - RTK_Cal->StaPos_flow_tmp.y) / pho_bk
					+ (Calculation.GPS_POSnVEL[center_GPS_prn - 1].SatPos.y - RTK_Cal->StaPos_flow_tmp.y) / pho_bj;
				A3 = -(Calculation.GPS_POSnVEL[GPS_prn_common[i] - 1].SatPos.z - RTK_Cal->StaPos_flow_tmp.z) / pho_bk
					+ (Calculation.GPS_POSnVEL[center_GPS_prn - 1].SatPos.z - RTK_Cal->StaPos_flow_tmp.z) / pho_bj;
				//L1
				A[i*(2 * GPS_num_common + 2 * BDS_num_common - 1)] = A1;
				A[i*(2 * GPS_num_common + 2 * BDS_num_common - 1) + 1] = A2;
				A[i*(2 * GPS_num_common + 2 * BDS_num_common - 1) + 2] = A3;
				A[i*(2 * GPS_num_common + 2 * BDS_num_common - 1) + i + 3] = lambda_L1;
				//L2
				A[(i + GPS_num_common - 1)*(2 * GPS_num_common + 2 * BDS_num_common - 1)] = A1;
				A[(i + GPS_num_common - 1)*(2 * GPS_num_common + 2 * BDS_num_common - 1) + 1] = A2;
				A[(i + GPS_num_common - 1)*(2 * GPS_num_common + 2 * BDS_num_common - 1) + 2] = A3;
				A[(i + GPS_num_common - 1)*(2 * GPS_num_common + 2 * BDS_num_common - 1) + i + 3 + GPS_num_common - 1] = lambda_L2;
				//P1
				A[(i + 2 * (GPS_num_common - 1) + 2 * (BDS_num_common - 1))*(2 * GPS_num_common + 2 * BDS_num_common - 1)] = A1;
				A[(i + 2 * (GPS_num_common - 1) + 2 * (BDS_num_common - 1))*(2 * GPS_num_common + 2 * BDS_num_common - 1) + 1] = A2;
				A[(i + 2 * (GPS_num_common - 1) + 2 * (BDS_num_common - 1))*(2 * GPS_num_common + 2 * BDS_num_common - 1) + 2] = A3;
				//P2
				A[(i + 3 * (GPS_num_common - 1) + 2 * (BDS_num_common - 1))*(2 * GPS_num_common + 2 * BDS_num_common - 1)] = A1;
				A[(i + 3 * (GPS_num_common - 1) + 2 * (BDS_num_common - 1))*(2 * GPS_num_common + 2 * BDS_num_common - 1) + 1] = A2;
				A[(i + 3 * (GPS_num_common - 1) + 2 * (BDS_num_common - 1))*(2 * GPS_num_common + 2 * BDS_num_common - 1) + 2] = A3;

				//随机模型
				for (int j = 0; j < GPS_num_common - 1; j++)
				{
					Q_GPS[i*(GPS_num_common - 1) + j] = Calculation.GPS_POSnVEL[center_GPS_prn - 1].sigma2;
				}
				Q_GPS[i*(GPS_num_common - 1) + i] = Q_GPS[i*(GPS_num_common - 1) + i] + Calculation.GPS_POSnVEL[GPS_prn_common[i] - 1].sigma2;
			}
			else
			{
				for (int j = i + 1; j < GPS_num_common; j++)
				{
					GPS_prn_common[j - 1] = GPS_prn_common[j];
				}
				GPS_prn_common[GPS_num_common - 1] = center_GPS_prn;
				i--;
			}
		}

		//BDS
		for (int i = 0; i < BDS_num_common - 1; i++)
		{
			if (BDS_prn_common[i] != center_BDS_prn)
			{
				pho_bk = sqrt(pow((RTK_Cal->StaPos_flow_tmp.x - Calculation.BDS_POSnVEL[BDS_prn_common[i] - 1].SatPos.x), 2)
					+ pow((RTK_Cal->StaPos_flow_tmp.y - Calculation.BDS_POSnVEL[BDS_prn_common[i] - 1].SatPos.y), 2)
					+ pow((RTK_Cal->StaPos_flow_tmp.z - Calculation.BDS_POSnVEL[BDS_prn_common[i] - 1].SatPos.z), 2));
				pho_bj = sqrt(pow((RTK_Cal->StaPos_flow_tmp.x - Calculation.BDS_POSnVEL[center_BDS_prn - 1].SatPos.x), 2)
					+ pow((RTK_Cal->StaPos_flow_tmp.y - Calculation.BDS_POSnVEL[center_BDS_prn - 1].SatPos.y), 2)
					+ pow((RTK_Cal->StaPos_flow_tmp.z - Calculation.BDS_POSnVEL[center_BDS_prn - 1].SatPos.z), 2));
				pho_ak = sqrt(pow((RTK_Cal->StaPos_ref.x - Calculation.BDS_POSnVEL[BDS_prn_common[i] - 1].SatPos.x), 2)
					+ pow((RTK_Cal->StaPos_ref.y - Calculation.BDS_POSnVEL[BDS_prn_common[i] - 1].SatPos.y), 2)
					+ pow((RTK_Cal->StaPos_ref.z - Calculation.BDS_POSnVEL[BDS_prn_common[i] - 1].SatPos.z), 2));
				pho_aj = sqrt(pow((RTK_Cal->StaPos_ref.x - Calculation.BDS_POSnVEL[center_BDS_prn - 1].SatPos.x), 2)
					+ pow((RTK_Cal->StaPos_ref.y - Calculation.BDS_POSnVEL[center_BDS_prn - 1].SatPos.y), 2)
					+ pow((RTK_Cal->StaPos_ref.z - Calculation.BDS_POSnVEL[center_BDS_prn - 1].SatPos.z), 2));
				int index_1, index_2;                   //BDS_prn_common[i]卫星在obs1和obs2中存储的位置
				int index_center_1, index_center_2;     //center_BDS_prn卫星在obs1和obs2中存储的位置
				for (int j = 0; j < Obs_1.BDS_SatNum; j++)
				{
					if (Obs_1.BDS_Sat[j].Prn == BDS_prn_common[i])
						index_1 = j;
					if (Obs_1.BDS_Sat[j].Prn == center_BDS_prn)
						index_center_1 = j;
				}
				for (int j = 0; j < Obs_2.BDS_SatNum; j++)
				{
					if (Obs_2.BDS_Sat[j].Prn == BDS_prn_common[i])
						index_2 = j;
					if (Obs_2.BDS_Sat[j].Prn == center_BDS_prn)
						index_center_2 = j;
				}

				L[i + 2 * (GPS_num_common - 1)] =
					lambda_B1*(Obs_1.BDS_Sat[index_1].Adr[0] - Obs_2.BDS_Sat[index_2].Adr[0]
					- Obs_1.BDS_Sat[index_center_1].Adr[0] + Obs_2.BDS_Sat[index_center_2].Adr[0])
					- pho_bk + pho_bj + pho_ak - pho_aj;   //BDS B1 相位
				L[i + 2 * (GPS_num_common - 1) + BDS_num_common - 1] =
					lambda_B3*(Obs_1.BDS_Sat[index_1].Adr[1] - Obs_2.BDS_Sat[index_2].Adr[1]
					- Obs_1.BDS_Sat[index_center_1].Adr[1] + Obs_2.BDS_Sat[index_center_2].Adr[1])
					- pho_bk + pho_bj + pho_ak - pho_aj;   //BDS B3 相位
				L[i + 4 * (GPS_num_common - 1) + 2 * (BDS_num_common - 1)] =
					Obs_1.BDS_Sat[index_1].Psr[0] - Obs_2.BDS_Sat[index_2].Psr[0]
					- Obs_1.BDS_Sat[index_center_1].Psr[0] + Obs_2.BDS_Sat[index_center_2].Psr[0]
					- pho_bk + pho_bj + pho_ak - pho_aj;   //BDS B1 伪距
				L[i + 4 * (GPS_num_common - 1) + 3 * (BDS_num_common - 1)] =
					Obs_1.BDS_Sat[index_1].Psr[1] - Obs_2.BDS_Sat[index_2].Psr[1]
					- Obs_1.BDS_Sat[index_center_1].Psr[1] + Obs_2.BDS_Sat[index_center_2].Psr[1]
					- pho_bk + pho_bj + pho_ak - pho_aj;   //BDS B3 伪距

				double A1, A2, A3;
				A1 = -(Calculation.BDS_POSnVEL[BDS_prn_common[i] - 1].SatPos.x - RTK_Cal->StaPos_flow_tmp.x) / pho_bk
					+ (Calculation.BDS_POSnVEL[center_BDS_prn - 1].SatPos.x - RTK_Cal->StaPos_flow_tmp.x) / pho_bj;
				A2 = -(Calculation.BDS_POSnVEL[BDS_prn_common[i] - 1].SatPos.y - RTK_Cal->StaPos_flow_tmp.y) / pho_bk
					+ (Calculation.BDS_POSnVEL[center_BDS_prn - 1].SatPos.y - RTK_Cal->StaPos_flow_tmp.y) / pho_bj;
				A3 = -(Calculation.BDS_POSnVEL[BDS_prn_common[i] - 1].SatPos.z - RTK_Cal->StaPos_flow_tmp.z) / pho_bk
					+ (Calculation.BDS_POSnVEL[center_BDS_prn - 1].SatPos.z - RTK_Cal->StaPos_flow_tmp.z) / pho_bj;
				//B1
				A[(i + 2 * (GPS_num_common - 1))*(2 * GPS_num_common + 2 * BDS_num_common - 1)] = A1;
				A[(i + 2 * (GPS_num_common - 1))*(2 * GPS_num_common + 2 * BDS_num_common - 1) + 1] = A2;
				A[(i + 2 * (GPS_num_common - 1))*(2 * GPS_num_common + 2 * BDS_num_common - 1) + 2] = A3;
				A[(i + 2 * (GPS_num_common - 1))*(2 * GPS_num_common + 2 * BDS_num_common - 1) + i + 3 + 2 * (GPS_num_common - 1)] = lambda_B1;
				//B3
				A[(i + 2 * (GPS_num_common - 1) + BDS_num_common - 1)*(2 * GPS_num_common + 2 * BDS_num_common - 1)] = A1;
				A[(i + 2 * (GPS_num_common - 1) + BDS_num_common - 1)*(2 * GPS_num_common + 2 * BDS_num_common - 1) + 1] = A2;
				A[(i + 2 * (GPS_num_common - 1) + BDS_num_common - 1)*(2 * GPS_num_common + 2 * BDS_num_common - 1) + 2] = A3;
				A[(i + 2 * (GPS_num_common - 1) + BDS_num_common - 1)*(2 * GPS_num_common + 2 * BDS_num_common - 1) + i + 3 + 2 * (GPS_num_common - 1) + BDS_num_common - 1] = lambda_B3;
				//P1
				A[(i + 4 * (GPS_num_common - 1) + 2 * (BDS_num_common - 1))*(2 * GPS_num_common + 2 * BDS_num_common - 1)] = A1;
				A[(i + 4 * (GPS_num_common - 1) + 2 * (BDS_num_common - 1))*(2 * GPS_num_common + 2 * BDS_num_common - 1) + 1] = A2;
				A[(i + 4 * (GPS_num_common - 1) + 2 * (BDS_num_common - 1))*(2 * GPS_num_common + 2 * BDS_num_common - 1) + 2] = A3;
				//P3
				A[(i + 4 * (GPS_num_common - 1) + 3 * (BDS_num_common - 1))*(2 * GPS_num_common + 2 * BDS_num_common - 1)] = A1;
				A[(i + 4 * (GPS_num_common - 1) + 3 * (BDS_num_common - 1))*(2 * GPS_num_common + 2 * BDS_num_common - 1) + 1] = A2;
				A[(i + 4 * (GPS_num_common - 1) + 3 * (BDS_num_common - 1))*(2 * GPS_num_common + 2 * BDS_num_common - 1) + 2] = A3;

				//随机模型
				for (int j = 0; j < BDS_num_common - 1; j++)
				{
					Q_BDS[i*(BDS_num_common - 1) + j] = Calculation.BDS_POSnVEL[center_BDS_prn - 1].sigma2;
				}
				Q_BDS[i*(BDS_num_common - 1) + i] = Q_BDS[i*(BDS_num_common - 1) + i] + Calculation.BDS_POSnVEL[BDS_prn_common[i] - 1].sigma2;
			}
			else
			{
				for (int j = i + 1; j < BDS_num_common; j++)
				{
					BDS_prn_common[j - 1] = BDS_prn_common[j];
				}
				BDS_prn_common[BDS_num_common - 1] = center_BDS_prn;
				i--;
			}
		}

		//随机模型建立 P矩阵行和列都是 4*(GPS_num_common+BDS_num_common-2)
		int Q_roll = 4 * (GPS_num_common + BDS_num_common - 2);
		int Qg_roll = GPS_num_common - 1;
		int Qb_roll = BDS_num_common - 1;
		memset(Q, 0, sizeof(Q));

		for (int i = 0; i < Qg_roll; i++)
		{
			for (int j = 0; j < Qg_roll; j++) //i行j列
			{
				Q[i*Q_roll + j] = Q_GPS[i*Qg_roll + j]; //f1相位
				Q[(i + Qg_roll)*Q_roll + (j + Qg_roll)] = Q_GPS[i*Qg_roll + j]; //f2相位
				Q[(i + 2 * Qg_roll + 2 * Qb_roll)*Q_roll + (j + 2 * Qg_roll + 2 * Qb_roll)] = sigma2_P_phi*Q_GPS[i*Qg_roll + j]; //f1伪距
				Q[(i + 3 * Qg_roll + 2 * Qb_roll)*Q_roll + (j + 3 * Qg_roll + 2 * Qb_roll)] = sigma2_P_phi*Q_GPS[i*Qg_roll + j]; //f2伪距
			}
		}
		for (int i = 0; i < Qb_roll; i++)
		{
			for (int j = 0; j < Qb_roll; j++) //i行j列
			{
				Q[(i + 2 * Qg_roll)*Q_roll + (j + 2 * Qg_roll)] = Q_BDS[i*Qb_roll + j]; //b1相位
				Q[(i + 2 * Qg_roll + Qb_roll)*Q_roll + (j + 2 * Qg_roll + Qb_roll)] = Q_BDS[i*Qb_roll + j]; //b3相位
				Q[(i + 4 * Qg_roll + 2 * Qb_roll)*Q_roll + (j + 4 * Qg_roll + 2 * Qb_roll)] = sigma2_P_phi*Q_BDS[i*Qb_roll + j]; //b1伪距
				Q[(i + 4 * Qg_roll + 3 * Qb_roll)*Q_roll + (j + 4 * Qg_roll + 3 * Qb_roll)] = sigma2_P_phi*Q_BDS[i*Qb_roll + j]; //b3伪距
			}
		}
		MatrixInv(4 * (Qg_roll + Qb_roll), Q, P);

		//解算
		MatrixTrans(4 * (Qg_roll + Qb_roll), 3 + 2 * (Qg_roll + Qb_roll), A, AT);
		MatrixMul(3 + 2 * (Qg_roll + Qb_roll), 4 * (Qg_roll + Qb_roll), 4 * (Qg_roll + Qb_roll), 4 * (Qg_roll + Qb_roll), AT, P, ATP);
		MatrixMul(3 + 2 * (Qg_roll + Qb_roll), 4 * (Qg_roll + Qb_roll), 4 * (Qg_roll + Qb_roll), 3 + 2 * (Qg_roll + Qb_roll), ATP, A, ATPA);
		MatrixInv(3 + 2 * (Qg_roll + Qb_roll), ATPA, N);
		MatrixMul(3 + 2 * (Qg_roll + Qb_roll), 3 + 2 * (Qg_roll + Qb_roll), 3 + 2 * (Qg_roll + Qb_roll), 4 * (Qg_roll + Qb_roll), N, AT, NAT);
		MatrixMul(3 + 2 * (Qg_roll + Qb_roll), 4 * (Qg_roll + Qb_roll), 4 * (Qg_roll + Qb_roll), 4 * (Qg_roll + Qb_roll), NAT, P, NATP);
		MatrixMul(3 + 2 * (Qg_roll + Qb_roll), 4 * (Qg_roll + Qb_roll), 4 * (Qg_roll + Qb_roll), 1, NATP, L, NATPL);


		if ((pow(NATPL[0], 2) + pow(NATPL[1], 2) + pow(NATPL[2], 2)) < 1e-5)
		{
			//输出结果
			RTK_Cal->SatNum = GPS_num_common + BDS_num_common;
			RTK_Cal->StaPos_flow = RTK_Cal->StaPos_flow_tmp;
			RTK_Cal->float_line.x = RTK_Cal->StaPos_flow.x - RTK_Cal->StaPos_ref.x;
			RTK_Cal->float_line.y = RTK_Cal->StaPos_flow.y - RTK_Cal->StaPos_ref.y;
			RTK_Cal->float_line.z = RTK_Cal->StaPos_flow.z - RTK_Cal->StaPos_ref.z;
			//精度评定
			memcpy(Q_float, N, sizeof(N));  // 3+2*(Qg_roll+Qb_roll)
			memset(RTK_Cal->Qxx, 0, sizeof(RTK_Cal->Qxx));
			RTK_Cal->Qxx[0] = N[0];
			RTK_Cal->Qxx[1] = N[1];
			RTK_Cal->Qxx[2] = N[2];
			RTK_Cal->Qxx[3] = N[2 * (Qg_roll + Qb_roll) + 3];
			RTK_Cal->Qxx[4] = N[2 * (Qg_roll + Qb_roll) + 3 + 1];
			RTK_Cal->Qxx[5] = N[2 * (Qg_roll + Qb_roll) + 3 + 2];
			RTK_Cal->Qxx[6] = N[2 * (2 * (Qg_roll + Qb_roll) + 3)];
			RTK_Cal->Qxx[7] = N[2 * (2 * (Qg_roll + Qb_roll) + 3) + 1];
			RTK_Cal->Qxx[8] = N[2 * (2 * (Qg_roll + Qb_roll) + 3) + 2];

			//lambda固定模糊度
			double* fa = NULL, *Qa = NULL, Qb[9] = { 0 }, *Qab = NULL;
			double *F = NULL, *s = NULL;
			int n, m = 2;
			double info;
			n = 2 * (Qg_roll + Qb_roll); //模糊度维数
			fa = mat(n, 1);
			Qa = mat(n, n);
			Qab = mat(n, 3);
			memset(Qa, 0, sizeof(Qa));
			memset(Qab, 0, sizeof(Qab));
			memset(Qb, 0, sizeof(Qb));
			for (int i = 0; i < n; i++)
			{
				fa[i] = result[i + 3];
				for (int j = 0; j < n; j++)
				{
					Qa[n*i + j] = N[(3 + n)*(3 + i) + 3 + j];
				}
				for (int j = 0; j < 3; j++)
				{
					Qab[3 * i + j] = N[(3 + n)*(3 + i) + j];
				}
			}
			Qb[0] = N[0];
			Qb[4] = N[n + 3 + 1];
			Qb[8] = N[2 * (n + 3) + 2];

			// 调用Lambda进行计算 
			// n： 模糊度维数
			// m： 模糊度候选解个数，设为2
			// fa：模糊度浮点解向量，n*1
			// Qa：模糊度方差协方差矩阵，n*n
			// F： 模糊度固定解向量，n*m
			// s： 模糊度残差二次型，1*m
			F = mat(n, m); s = mat(1, m);
			info = lambda(n, m, fa, Qa, F, s);
			double fix[2 * SatComMax];

			//更新
			double Qba[2 * SatComMax * 3], QaN[2 * SatComMax * 2 * SatComMax];
			double delta_a[2 * SatComMax];
			double Q1[2 * SatComMax * 3], Q2[3], Q3[9];
			for (int i = 0; i < n; i++)
			{
				delta_a[i] = fa[i] - F[i];
			}
			MatrixTrans(n, 3, Qab, Qba); //Qba:3*n
			MatrixInv(n, Qa, QaN);
			MatrixMul(3, n, n, n, Qba, QaN, Q1);
			MatrixMul(3, n, n, 1, Q1, delta_a, Q2);	

			//协因数矩阵更新
			MatrixMul(3, n, n, 3, Q1, Qab, Q3);
			MatrixSub(3, 3, 3, 3, RTK_Cal->Qxx, Q3, RTK_Cal->Qxx_fixed);

			RTK_Cal->fixed_line.x = RTK_Cal->float_line.x - Q2[0];
			RTK_Cal->fixed_line.y = RTK_Cal->float_line.y - Q2[1];
			RTK_Cal->fixed_line.z = RTK_Cal->float_line.z - Q2[2];
			RTK_Cal->ratio = s[1] / s[0];
			if (info == -1)
			{
				RTK_Cal->ratio = 0;
			}

			/*
			ofstream outfile("C:\\gnss\\RTK\\R.txt", ios::out);
			//输出
			for (int i = 0; i < Qg_roll; i++)
			{
				outfile << fa[i] << "  ";
			}
			outfile << endl;
			for (int i = 0; i < Qg_roll; i++)
			{
				outfile << F[i] << "  ";
			}
			*/

			break;
		}

		//迭代
		else
		{
			RTK_Cal->StaPos_flow_tmp.x = RTK_Cal->StaPos_flow_tmp.x + NATPL[0];
			RTK_Cal->StaPos_flow_tmp.y = RTK_Cal->StaPos_flow_tmp.y + NATPL[1];
			RTK_Cal->StaPos_flow_tmp.z = RTK_Cal->StaPos_flow_tmp.z + NATPL[2];			
			for (int i = 0; i < 3 + 2 * (Qg_roll + Qb_roll); i++)
			{
				result[i] = result[i] + NATPL[i];
			}
		}

		if (times == 19)
		{
			cout << "最小二乘失败" << endl;
			system("PAUSE");
		}
	}

	return 0;
}