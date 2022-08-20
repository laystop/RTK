#include "stdafx.h"
#include <iostream>
#include<cmath>
using namespace std;

/*************************************
BDST -> GPST
BDST + 1356week +14second = GPST
*************************************/


/*************************************
CommonTimeToMJDTime
目的：通用时到简化儒略日的转换算法

参数：
year     年（通用时）
month    月（通用时）
day      日（通用时）
hour     时（通用时）
minute   分（通用时）
second   秒（通用时）
days     整数天（简化儒略日）
fracday  小数天（简化儒略日）

返回值：1=正常 0=输入时间为负
**************************************/
int CommonTimeToMJDTime(int &year, int &month, int &day, int &hour, int &minute, double &second, int &days, double &fracday)
{
	double ut;
	if (year < 0 || month < 0 || day < 0 || hour < 0 || minute < 0 || second < 0)
	{
		cout << "Error input CommonTime parameter!";
		return 0;
	}
	else
	{
		if (month <= 2)
		{
			year = year - 1;
			month = month + 12;
		}
		days = floor(365.25*year) + floor(30.6001*(month + 1)) + day - 679019;
		ut = hour + minute / 60.0 + second / 3600.0;
		days = days + floor(ut / 24.0);
		fracday = ut / 24.0 - floor(ut / 24.0);

		return 1;
	}
}

/*************************************
MJDTimeToCommonTime
目的：简化儒略日到通用时的转换算法

参数：
year     年（通用时）
month    月（通用时）
day      日（通用时）
hour     时（通用时）
minute   分（通用时）
second   秒（通用时）
days     整数天（简化儒略日）
fracday  小数天（简化儒略日）

返回值：1=正常 0=输入时间为负
**************************************/
int MJDTimeToCommonTime(int &days, double &fracday, int &year, int &month, int &day, int &hour, int &minute, double &second)
{
	int a, b, c, d, e;
	double D;
	if (days < 0 || fracday < 0)
	{
		cout << "Error input MJDTime parameter!";
		return 0;
	}
	else
	{
		a = floor(days + fracday + 2400001);
		b = a + 1537;
		c = floor((b - 122.1) / 365.25);
		d = floor(365.25*c);
		e = floor((b - d) / 30.6001);
		D = b - d - floor(30.6001*e) + days + fracday + 2400001 - floor(days + fracday + 2400001);
		month = e - 1 - 12 * floor(e / 14.0);
		year = c - 4715 - floor((7 + month) / 10.0);
		day = floor(D);
		hour = floor((D - day)*24.0);
		minute = floor(((D - day)*24.0 - hour)*60.0);
		second = (((D - day)*24.0 - hour)*60.0 - minute)*60.0;

		return 1;
	}
}

/*************************************
GPSTimeToMJDTime
目的：GPS时到简化儒略日的转换算法

参数：
week           周（GPS）
secondofweek   周内秒（GPS）
days           整数天（简化儒略日）
fracday        小数天（简化儒略日）

返回值：1=正常 0=输入时间为负
**************************************/
int GPSTimeToMJDTime(int &week, double &secondofweek, int &days, double &fracday)
{
	if (week < 0 || secondofweek < 0)
	{
		cout << "Error input GPSTime parameter!";
		return 0;
	}
	else
	{
		days = 44244 + week * 7 + floor(secondofweek / 86400.0);
		fracday = secondofweek / 86400.0 - floor(secondofweek / 86400.0);

		return 1;
	}
}

/*************************************
MJDTimeToGPSTime
目的：简化儒略日到GPS时的转换算法

参数：
days           整数天（简化儒略日）
fracday        小数天（简化儒略日）
week           周（GPS）
secondofweek   周内秒（GPS）

返回值：1=正常 0=输入时间为负
**************************************/
int MJDTimeToGPSTime(int &days, double &fracday, int &week, double &secondofweek)
{
	if (days < 0 || fracday < 0)
	{
		cout << "Error input MJDTime parameter!";
		return 0;
	}
	else
	{
		week = floor(((days + fracday) - 44244) / 7.0);
		secondofweek = ((days + fracday) - 44244 - week * 7) * 86400;

		return 1;
	}
}

/*************************************
CommonTimeToGPSTime
目的：通用时到GPS时的转换算法

参数：
year     年（通用时）
month    月（通用时）
day      日（通用时）
hour     时（通用时）
minute   分（通用时）
second   秒（通用时）
week           周（GPS）
secondofweek   周内秒（GPS）

返回值：1=正常 0=输入时间为负
**************************************/
int CommonTimeToGPSTime(int &year, int &month, int &day, int &hour, int &minute, double &second, int &week, double &secondofweek)
{
	int days;
	double fracday;
	if (CommonTimeToMJDTime(year, month, day, hour, minute, second, days, fracday) == 1)
	{
		if (MJDTimeToGPSTime(days, fracday, week, secondofweek) == 1)
		{
			return 1;
		}
	}
	else
		return 0;
}