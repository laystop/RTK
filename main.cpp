// RTK.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "struct.h"
#include<iostream>
#include<fstream>
#include<string>
#include<map>
#include<iomanip>
#include<sstream>

using namespace std;

int main()
{
	string filename[3];
	//filename[0] = "C:\\gnss\\RTK\\计算参考结果\\20210412_Zero_OEM7\\202104122230-7200.21O";
	//filename[1] = "C:\\gnss\\RTK\\计算参考结果\\20210412_Zero_OEM7\\202104122230-7190.21O";
	//filename[2] = "C:\\gnss\\RTK\\计算参考结果\\20210412_Zero_OEM7\\brdm1020.21p";

	filename[0] = "C:\\gnss\\RTK\\计算参考结果\\20210515_Short_OEM7\\CENT_20210515.21O";
	filename[1] = "C:\\gnss\\RTK\\计算参考结果\\20210515_Short_OEM7\\SGG_20210515.21O";
	filename[2] = "C:\\gnss\\RTK\\计算参考结果\\20210515_Short_OEM7\\brdm1350.21p";

	OBS Obs[2];
	map<int, map<GPSTIME, EPHEM>> map_eph; //GPS 32, BDS 63, prn号一共95个（GPS:n, BDS:32+n）
	double delta_t; //时间同步
	SPP_CAL SPP_Cal[2];
	CALCULATION Calculation[2];
	RTK_CAL RTK_Cal;
	double RDOP;

	//读观测文件
	ifstream fin_0(filename[0]), fin_1(filename[1]);
	ReadHeader_RENIX(fin_0,&Obs[0]); 
	ReadHeader_RENIX(fin_1, &Obs[1]);
	//读星历
	ifstream fin_2(filename[2]);
	ReadHeader_RENIX(fin_2, &Obs[0]);
	ReadEph_RENIX(fin_2, map_eph);
	//输出
	ofstream outfile("C:\\gnss\\RTK\\515.txt", ios::out);
	double n, e, u, nf, ef, uf;

	while (!fin_0.eof() && !fin_1.eof())
	{
		InitObs(&Obs[0]);
		InitObs(&Obs[1]);

		ReadObs_RENIX(fin_0, &Obs[0]); //读流动站观测文件
		ReadObs_RENIX(fin_1, &Obs[1]); //读基准站观测文件

		//时间同步		
		delta_t = (Obs[0].ObsTime.Week - Obs[1].ObsTime.Week) * 604800 + Obs[0].ObsTime.SecofWeek - Obs[1].ObsTime.SecofWeek;
		while (abs(delta_t) > 1e-5)
		{
			if (fin_0.eof() || fin_1.eof()) //文件尾
			{
				return 0;
			}
			if (delta_t > 0) //基准站落后
			{
				InitObs(&Obs[1]);
				ReadObs_RENIX(fin_1, &Obs[1]);
			}
			else if (delta_t < 0) //流动站落后
			{
				InitObs(&Obs[0]);
				ReadObs_RENIX(fin_0, &Obs[0]);
			}
			delta_t = (Obs[0].ObsTime.Week - Obs[1].ObsTime.Week) * 604800 + Obs[0].ObsTime.SecofWeek - Obs[1].ObsTime.SecofWeek;
		}


		//循环，逐一计算挨个历元
		memset(&Calculation[0], 0, sizeof(CALCULATION));
		memset(&Calculation[1], 0, sizeof(CALCULATION));
		memset(&SPP_Cal[0], 0, sizeof(SPP_CAL));
		memset(&SPP_Cal[1], 0, sizeof(SPP_CAL));
		//单点定位
		PreEdit(&Obs[0], map_eph);
		SPP(Obs[0], &SPP_Cal[0], map_eph, &Calculation[0]);
		PreEdit(&Obs[1], map_eph);
		//SPP(Obs[1], &SPP_Cal[1], map_eph, &Calculation[1]);		

		RTK_Cal.StaPos_ref.x = -2267810.196;
		RTK_Cal.StaPos_ref.y = 5009356.572;
		RTK_Cal.StaPos_ref.z = 3221000.818;

		//RTK
		RTK_Cal.time = SPP_Cal[0].time;
		RTK_Cal.StaPos_flow_tmp = SPP_Cal[0].RecPos; //基准站和流动站分别是1,0
		//RTK_Cal.StaPos_ref = Obs[1].StaPos_ref;    //基准站直接从文件头读
		RTK(Obs[0], Obs[1], Calculation[0], &RTK_Cal);

		if (RTK_Cal.ratio > 2) //固定解
		{
			XYZToNEU(RTK_Cal.fixed_line, RTK_Cal.StaPos_ref, &nf, &ef, &uf);
			XYZToNEU(RTK_Cal.float_line, RTK_Cal.StaPos_ref, &n, &e, &u);
			RDOP = sqrt(RTK_Cal.Qxx_fixed[0] + RTK_Cal.Qxx_fixed[4] + RTK_Cal.Qxx_fixed[8]);
			//输出
			/*outfile << fixed << setprecision(5) << RTK_Cal.time.Week << "   " << RTK_Cal.time.SecofWeek << setw(12)
				<< RTK_Cal.float_line.x << setw(12) << RTK_Cal.float_line.y << setw(12) << RTK_Cal.float_line.z << setw(12)
				<< RTK_Cal.fixed_line.x << setw(12) << RTK_Cal.fixed_line.y << setw(12) << RTK_Cal.fixed_line.z
				<< setw(12) << n << setw(12) << e << setw(12) << u << setw(12) << RTK_Cal.ratio << endl;*/
			outfile << fixed << setprecision(5) << RTK_Cal.time.Week << "   " << RTK_Cal.time.SecofWeek << setw(12)
				<< ef << setw(12) << nf << setw(12) << uf << setw(12) << e << setw(12) << n << setw(12) << u << setw(12) << RTK_Cal.ratio << endl;
		}
		else if (RTK_Cal.ratio < 2)
		{
			XYZToNEU(RTK_Cal.float_line, RTK_Cal.StaPos_ref, &n, &e, &u);
			//RDOP = sqrt(RTK_Cal.Qxx[0] + RTK_Cal.Qxx[4] + RTK_Cal.Qxx[8]);
			//输出
			/*outfile << fixed << setprecision(5) << RTK_Cal.time.Week << "   " << RTK_Cal.time.SecofWeek << setw(12)
				<< RTK_Cal.float_line.x << setw(12) << RTK_Cal.float_line.y << setw(12) << RTK_Cal.float_line.z
				<< setw(12) << RTK_Cal.float_line.x << setw(12) << RTK_Cal.float_line.y << setw(12) << RTK_Cal.float_line.z
				<< setw(12) << n << setw(12) << e << setw(12) << u << setw(12) << RTK_Cal.ratio << endl;*/
			outfile << fixed << setprecision(5) << RTK_Cal.time.Week << "   " << RTK_Cal.time.SecofWeek << setw(12)
				<< e << setw(12) << n << setw(12) << u << setw(12) << e << setw(12) << n << setw(12) << u << setw(12) << RTK_Cal.ratio << endl;
		}


	}

	fin_0.close();
	fin_1.close();
	fin_2.close();
	outfile.close();
    return 0;
}