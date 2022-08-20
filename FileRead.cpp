#include "stdafx.h"
#include "struct.h"
#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<map>

using namespace std;


/*************************************
ReadHeader_RENIX
目的：读 RENIX3.04 文件头

参数：


返回值：1=观测文件，2=星历，3=无法成功读取
**************************************/
int ReadHeader_RENIX(ifstream &fin, OBS *obs)
{
	string oneline, sys, oneline2;
	int flag = 0;
	int obsnum = 0;
	int lines = 1;
	bool first_flag = true;
	int pos;
	string strs[30];
	while (getline(fin, oneline))
	{
		if (oneline.find("OBSERVATION DATA", 0) != string::npos)
		{
			flag = 1;
			continue;
		}

		if (oneline.find("NAV DATA", 0) != string::npos)
		{
			flag = 2;
			continue;
		}

		if (oneline.find("SYS / # / OBS TYPES", 0) != string::npos)
		{
			sys = oneline[0];

			obsnum = atoi((oneline.substr(4, 2)).c_str());
			lines = floor(obsnum / 13) + 1;
			if (obsnum % 13 == 0)
			{
				lines--;
			}
			int nl = lines;
			while (nl > 1)
			{
				getline(fin, oneline2);
				oneline = oneline + oneline2;
				nl--;
			}

			if (sys == "G")
			{
				while (lines > 0)
				{
					pos = oneline.find("SYS / # / OBS TYPES", 0);
					oneline.erase(pos, 19);
					lines--;
				}
				oneline.erase(0, 7);
				int index_temp = -1;
				oneline2 = Remove_space(oneline);
				istringstream iss(oneline2);
				string str_one;
				while (getline(iss,str_one,' '))
				{
					index_temp++;
					if (str_one == "C1C")
						obs->GPS_L1C_index = index_temp;
					if (str_one == "L1C")
						obs->GPS_L1P_index = index_temp;
					if (str_one == "D1C")
						obs->GPS_L1D_index = index_temp;

					if (str_one == "C2W")
						obs->GPS_L2C_index = index_temp;
					if (str_one == "L2W")
						obs->GPS_L2P_index = index_temp;
					if (str_one == "D2W")
						obs->GPS_L2D_index = index_temp;
				}


			}
			else if (sys == "C")
			{
				while (lines > 0)
				{
					pos = oneline.find("SYS / # / OBS TYPES", 0);
					oneline.erase(pos, 19);
					lines--;
				}
				oneline.erase(0, 7);
				int index_temp = -1;
				oneline2 = Remove_space(oneline);
				istringstream iss(oneline2);
				string str_one;
				while (getline(iss, str_one, ' '))
				{
					index_temp++;
					if (str_one == "C2I")
						obs->BDS_B1C_index = index_temp;
					if (str_one == "L2I")
						obs->BDS_B1P_index = index_temp;
					if (str_one == "D2I")
						obs->BDS_B1D_index = index_temp;

					if (str_one == "C6I")
						obs->BDS_B3C_index = index_temp;
					if (str_one == "L6I")
						obs->BDS_B3P_index = index_temp;
					if (str_one == "D6I")
						obs->BDS_B3D_index = index_temp;
				}
			}
		}

		if (oneline.find("END OF HEADER", 0) != string::npos)
		{
			return flag;
		}
	}

	
}


/*************************************
ReadObs_RENIX
目的：读 RENIX3.04 观测文件

参数：
fin  文件流
Obs  观测数据

返回值：0=正常
**************************************/
int ReadObs_RENIX(ifstream &fin, OBS *Obs)
{
	string oneline;
	stringstream ss;
	char skip;
	int LineNum; //行数

	Obs->GPS_SatNum = 0;
	Obs->BDS_SatNum = 0;

	while (getline(fin, oneline))
	{
		if (oneline[0] == '>') //找到观测历元开头
		{
			ss.str(oneline);
			ss >> skip;
			ss >> Obs->ObsTime_ymd.Year;
			ss >> Obs->ObsTime_ymd.Month;
			ss >> Obs->ObsTime_ymd.Day;
			ss >> Obs->ObsTime_ymd.Hour;
			ss >> Obs->ObsTime_ymd.Minute;
			ss >> Obs->ObsTime_ymd.Second;
			ss >> Obs->SatNum;
			ss >> Obs->SatNum;
			CommonTimeToGPSTime(Obs->ObsTime_ymd.Year, Obs->ObsTime_ymd.Month, Obs->ObsTime_ymd.Day,
				Obs->ObsTime_ymd.Hour, Obs->ObsTime_ymd.Minute, Obs->ObsTime_ymd.Second,
				Obs->ObsTime.Week, Obs->ObsTime.SecofWeek);
			LineNum = -1;
			if (Obs->ObsTime_ymd.Second == 12)
				int i = 0;
		}

		else if (oneline[0] == 'G') //GPS
		{
			Obs->GPS_Sat[Obs->GPS_SatNum].Prn = atoi((oneline.substr(1, 2)).c_str()); //prn
			Obs->GPS_Sat[Obs->GPS_SatNum].Psr[0] = atof((oneline.substr(16 * Obs->GPS_L1C_index + 3, 15)).c_str()); //L1伪距
			Obs->GPS_Sat[Obs->GPS_SatNum].Adr[0] = atof((oneline.substr(16 * Obs->GPS_L1P_index + 3, 15)).c_str());//L1相位
			Obs->GPS_Sat[Obs->GPS_SatNum].Dop[0] = atof((oneline.substr(16 * Obs->GPS_L1D_index + 3, 15)).c_str());//L1多普勒
			Obs->GPS_Sat[Obs->GPS_SatNum].Psr[1] = atof((oneline.substr(16 * Obs->GPS_L2C_index + 3, 15)).c_str()); //L2伪距
			Obs->GPS_Sat[Obs->GPS_SatNum].Adr[1] = atof((oneline.substr(16 * Obs->GPS_L2P_index + 3, 15)).c_str());//L2相位
			Obs->GPS_Sat[Obs->GPS_SatNum].Dop[1] = atof((oneline.substr(16 * Obs->GPS_L2D_index + 3, 15)).c_str());//L2多普勒

			//如果双频观测数据缺失
			if (abs(Obs->GPS_Sat[Obs->GPS_SatNum].Psr[0]) < 1e-5 || abs(Obs->GPS_Sat[Obs->GPS_SatNum].Psr[1]) < 1e-5)
			{
				Obs->GPS_SatNum = Obs->GPS_SatNum - 1;
			}
			Obs->GPS_SatNum++;
		}
		else if (oneline[0] == 'C') //BDS
		{
			Obs->BDS_Sat[Obs->BDS_SatNum].Prn = atoi((oneline.substr(1, 2)).c_str()); //prn
			Obs->BDS_Sat[Obs->BDS_SatNum].Psr[0] = atof((oneline.substr(16 * Obs->BDS_B1C_index + 3, 15)).c_str());//B1伪距
			Obs->BDS_Sat[Obs->BDS_SatNum].Adr[0] = atof((oneline.substr(16 * Obs->BDS_B1P_index + 3, 15)).c_str());//B1相位
			Obs->BDS_Sat[Obs->BDS_SatNum].Dop[0] = atof((oneline.substr(16 * Obs->BDS_B1D_index + 3, 15)).c_str());//B1多普勒
			Obs->BDS_Sat[Obs->BDS_SatNum].Psr[1] = atof((oneline.substr(16 * Obs->BDS_B3C_index + 3, 15)).c_str());//B3伪距
			Obs->BDS_Sat[Obs->BDS_SatNum].Adr[1] = atof((oneline.substr(16 * Obs->BDS_B3P_index + 3, 15)).c_str());//B3相位
			Obs->BDS_Sat[Obs->BDS_SatNum].Dop[1] = atof((oneline.substr(16 * Obs->BDS_B3D_index + 3, 15)).c_str());//B3多普勒

			if (abs(Obs->BDS_Sat[Obs->BDS_SatNum].Psr[0]) < 1e-5 || abs(Obs->BDS_Sat[Obs->BDS_SatNum].Psr[1]) < 1e-5)
			{
				Obs->BDS_SatNum = Obs->BDS_SatNum - 1;
			}
			Obs->BDS_SatNum++;
		}
		LineNum++;
		if (LineNum == Obs->SatNum)
			break;
	}
	return 0;
}

/*************************************
ReadEph_RENIX
目的：读 RENIX3.04 星历文件

参数：
fin      文件流
map_eph  星历数据

返回值：0=正常
**************************************/
int ReadEph_RENIX(ifstream &fin, map<int, map<GPSTIME, EPHEM>> &map_eph)
{
	string oneline;
	stringstream ss;
	string sys_prn;
	int prn;
	int y, m, d, h, min;
	double s;
	GPSTIME t;

	while (getline(fin, oneline))
	{
		ss.clear();
		ss.str(oneline);
		ss >> sys_prn;   //prn
		prn = atoi(sys_prn.substr(1, 2).c_str());
		if (sys_prn[0] == 'C') prn = prn + 32;
		else if (sys_prn[0] == 'G') prn = prn;
		else continue;  //其它系统的跳过

		ss >> y >> m >> d >> h >> min >> s; //时间
		CommonTimeToGPSTime(y, m, d, h, min, s, t.Week, t.SecofWeek);
		ss.clear(ios::goodbit);

		map_eph[prn][t].af0 = atof(oneline.substr(23, 19).c_str());
		map_eph[prn][t].af1 = atof(oneline.substr(42, 19).c_str());
		map_eph[prn][t].af2 = atof(oneline.substr(61, 19).c_str());

		getline(fin, oneline); //第二行
		map_eph[prn][t].IODE = atof(oneline.substr(4, 19).c_str());
		map_eph[prn][t].crs = atof(oneline.substr(23, 19).c_str());
		map_eph[prn][t].deltaN = atof(oneline.substr(42, 19).c_str());
		map_eph[prn][t].M0 = atof(oneline.substr(61, 19).c_str());
		
		getline(fin, oneline); //第三行
		map_eph[prn][t].cuc = atof(oneline.substr(4, 19).c_str());
		map_eph[prn][t].ecc = atof(oneline.substr(23, 19).c_str());
		map_eph[prn][t].cus = atof(oneline.substr(42, 19).c_str());
		map_eph[prn][t].sqrt_A = atof(oneline.substr(61, 19).c_str());

		getline(fin, oneline); //第四行
		map_eph[prn][t].toe = atof(oneline.substr(4, 19).c_str());
		map_eph[prn][t].cic = atof(oneline.substr(23, 19).c_str());
		map_eph[prn][t].omega_0 = atof(oneline.substr(42, 19).c_str());
		map_eph[prn][t].cis = atof(oneline.substr(61, 19).c_str());

		getline(fin, oneline); //第五行
		map_eph[prn][t].i0 = atof(oneline.substr(4, 19).c_str());
		map_eph[prn][t].crc = atof(oneline.substr(23, 19).c_str());
		map_eph[prn][t].omega = atof(oneline.substr(42, 19).c_str());
		map_eph[prn][t].omega_dot = atof(oneline.substr(61, 19).c_str());

		getline(fin, oneline); //第六行
		map_eph[prn][t].IDOT = atof(oneline.substr(4, 19).c_str());
		//map_eph[prn][t].crc = atof(oneline.substr(23, 19).c_str());
		map_eph[prn][t].week = atof(oneline.substr(42, 19).c_str());
		//map_eph[prn][t].omega_dot = atof(oneline.substr(61, 19).c_str());

		getline(fin, oneline); //第七行
		//map_eph[prn][t].i0 = atof(oneline.substr(4, 19).c_str());
		//map_eph[prn][t].crc = atof(oneline.substr(23, 19).c_str());
		map_eph[prn][t].tgd = atof(oneline.substr(42, 19).c_str());
		map_eph[prn][t].iodc = atof(oneline.substr(61, 19).c_str());

		//getline(fin, oneline); //第八行
		//map_eph[prn][t].i0 = atof(oneline.substr(4, 19).c_str());
		//map_eph[prn][t].crc = atof(oneline.substr(23, 19).c_str());
		if (prn > 32)
		{
			map_eph[prn][t].week = map_eph[prn][t].week + 1356;
			//map_eph[prn][t].toe = map_eph[prn][t].toe + 14;
		}
	}

	return 0;
}

string Remove_space(const string& src)
{
	string result = "";
	for (int i = 0; src[i] != '\0'; i++)
	{
		if (src[i] != ' ')
			result.append(1, src[i]);
		else
			if (src[i + 1] != ' ')
				result.append(1, src[i]);
	}
	return result;
}

/*除了文件头里数据位置以外的变量全部清零*/
void InitObs(OBS *Obs)
{
	memset(&(Obs->GPS_Sat), 0, sizeof(SAT)*GPSsatnum);
	memset(&(Obs->BDS_Sat), 0, sizeof(SAT)*BDSsatnum);
	Obs->SatNum = 0;
	Obs->GPS_SatNum = 0;
	Obs->BDS_SatNum = 0;
	memset(&(Obs->ObsTime), 0, sizeof(GPSTIME));
	memset(&(Obs->ObsTime_ymd), 0, sizeof(COMMONTIME));
	memset(&(Obs->StaPos_ref), 0, sizeof(XYZ));
}