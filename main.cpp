// RTK.cpp : �������̨Ӧ�ó������ڵ㡣
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
	//filename[0] = "C:\\gnss\\RTK\\����ο����\\20210412_Zero_OEM7\\202104122230-7200.21O";
	//filename[1] = "C:\\gnss\\RTK\\����ο����\\20210412_Zero_OEM7\\202104122230-7190.21O";
	//filename[2] = "C:\\gnss\\RTK\\����ο����\\20210412_Zero_OEM7\\brdm1020.21p";

	filename[0] = "C:\\gnss\\RTK\\����ο����\\20210515_Short_OEM7\\CENT_20210515.21O";
	filename[1] = "C:\\gnss\\RTK\\����ο����\\20210515_Short_OEM7\\SGG_20210515.21O";
	filename[2] = "C:\\gnss\\RTK\\����ο����\\20210515_Short_OEM7\\brdm1350.21p";

	OBS Obs[2];
	map<int, map<GPSTIME, EPHEM>> map_eph; //GPS 32, BDS 63, prn��һ��95����GPS:n, BDS:32+n��
	double delta_t; //ʱ��ͬ��
	SPP_CAL SPP_Cal[2];
	CALCULATION Calculation[2];
	RTK_CAL RTK_Cal;
	double RDOP;

	//���۲��ļ�
	ifstream fin_0(filename[0]), fin_1(filename[1]);
	ReadHeader_RENIX(fin_0,&Obs[0]); 
	ReadHeader_RENIX(fin_1, &Obs[1]);
	//������
	ifstream fin_2(filename[2]);
	ReadHeader_RENIX(fin_2, &Obs[0]);
	ReadEph_RENIX(fin_2, map_eph);
	//���
	ofstream outfile("C:\\gnss\\RTK\\515.txt", ios::out);
	double n, e, u, nf, ef, uf;

	while (!fin_0.eof() && !fin_1.eof())
	{
		InitObs(&Obs[0]);
		InitObs(&Obs[1]);

		ReadObs_RENIX(fin_0, &Obs[0]); //������վ�۲��ļ�
		ReadObs_RENIX(fin_1, &Obs[1]); //����׼վ�۲��ļ�

		//ʱ��ͬ��		
		delta_t = (Obs[0].ObsTime.Week - Obs[1].ObsTime.Week) * 604800 + Obs[0].ObsTime.SecofWeek - Obs[1].ObsTime.SecofWeek;
		while (abs(delta_t) > 1e-5)
		{
			if (fin_0.eof() || fin_1.eof()) //�ļ�β
			{
				return 0;
			}
			if (delta_t > 0) //��׼վ���
			{
				InitObs(&Obs[1]);
				ReadObs_RENIX(fin_1, &Obs[1]);
			}
			else if (delta_t < 0) //����վ���
			{
				InitObs(&Obs[0]);
				ReadObs_RENIX(fin_0, &Obs[0]);
			}
			delta_t = (Obs[0].ObsTime.Week - Obs[1].ObsTime.Week) * 604800 + Obs[0].ObsTime.SecofWeek - Obs[1].ObsTime.SecofWeek;
		}


		//ѭ������һ���㰤����Ԫ
		memset(&Calculation[0], 0, sizeof(CALCULATION));
		memset(&Calculation[1], 0, sizeof(CALCULATION));
		memset(&SPP_Cal[0], 0, sizeof(SPP_CAL));
		memset(&SPP_Cal[1], 0, sizeof(SPP_CAL));
		//���㶨λ
		PreEdit(&Obs[0], map_eph);
		SPP(Obs[0], &SPP_Cal[0], map_eph, &Calculation[0]);
		PreEdit(&Obs[1], map_eph);
		//SPP(Obs[1], &SPP_Cal[1], map_eph, &Calculation[1]);		

		RTK_Cal.StaPos_ref.x = -2267810.196;
		RTK_Cal.StaPos_ref.y = 5009356.572;
		RTK_Cal.StaPos_ref.z = 3221000.818;

		//RTK
		RTK_Cal.time = SPP_Cal[0].time;
		RTK_Cal.StaPos_flow_tmp = SPP_Cal[0].RecPos; //��׼վ������վ�ֱ���1,0
		//RTK_Cal.StaPos_ref = Obs[1].StaPos_ref;    //��׼վֱ�Ӵ��ļ�ͷ��
		RTK(Obs[0], Obs[1], Calculation[0], &RTK_Cal);

		if (RTK_Cal.ratio > 2) //�̶���
		{
			XYZToNEU(RTK_Cal.fixed_line, RTK_Cal.StaPos_ref, &nf, &ef, &uf);
			XYZToNEU(RTK_Cal.float_line, RTK_Cal.StaPos_ref, &n, &e, &u);
			RDOP = sqrt(RTK_Cal.Qxx_fixed[0] + RTK_Cal.Qxx_fixed[4] + RTK_Cal.Qxx_fixed[8]);
			//���
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
			//���
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