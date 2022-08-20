#include "stdafx.h"
#include <iostream>
using namespace std;

/*************************************
MatrixAdd
目的：矩阵加法

参数：
line_a  矩阵a行数
row_a   矩阵a列数
lin_b   矩阵b行数
row_b   矩阵b列数
a,b     输入矩阵
c       输出矩阵
返回值：1=正常，0=致命错误
**************************************/
int MatrixAdd(int line_a, int row_a, int line_b, int row_b, double a[], double b[], double c[])
{
	int i, j, line, row;
	if (line_a != line_b || row_a != row_b)
	{
		cout << "Error dimension in MatrixAdd!" << endl;
		return 0;
	}
	else
	{
		line = line_a;
		row = row_a;
		for (i = 0; i < line; i++)
		{
			for (j = 0; j < row; j++)
			{
				c[i*row + j] = a[i*row + j] + b[i*row + j];
			}
		}
		return 1;
	}
}

/*************************************
MatrixSub
目的：矩阵减法

参数：
line_a  矩阵a行数
row_a   矩阵a列数
lin_b   矩阵b行数
row_b   矩阵b列数
a,b     输入矩阵
c       输出矩阵
返回值：1=正常，0=致命错误
**************************************/
int MatrixSub(int line_a, int row_a, int line_b, int row_b, double a[], double b[], double c[])
{
	int i, j, line, row;
	if (line_a != line_b || row_a != row_b)
	{
		cout << "Error dimension in MatrixSub!" << endl;
		return 0;
	}
	else
	{
		line = line_a;
		row = row_a;
		for (i = 0; i < line; i++)
		{
			for (j = 0; j < row; j++)
			{
				c[i*row + j] = a[i*row + j] - b[i*row + j];
			}
		}
		return 1;
	}
}

/*************************************
MatrixMul
目的：矩阵乘法

参数：
a,b     输入矩阵
c       输出矩阵
line_a  矩阵a行数
row_a   矩阵a列数
lin_b   矩阵b行数
row_b   矩阵b列数
返回值：1=正常，0=致命错误
**************************************/
int MatrixMul(int line_a, int row_a, int line_b, int row_b, double a[], double b[], double c[])
{
	int i, j, k, m, n, q;
	if (row_a != line_b)
	{
		cout << "Error dimension in MatrixMul!" << endl;
		return 0;
	}
	else
	{
		m = line_a;
		n = row_a;
		q = row_b;
		for (i = 0; i < m; i++)
		{
			for (j = 0; j < q; j++)
			{
				c[i*q + j] = 0;
				for (k = 0; k < n; k++)
				{
					c[i*q + j] = c[i*q + j] + a[i*n + k] * b[k*q + j];
				}
			}
		}
		return 1;
	}
}

/*************************************
MatrixInv
目的：矩阵求逆

参数：
n  矩阵行数和列数
a  输入矩阵
b  输出矩阵，b=inv(a)
返回值：1=正常，0=致命错误
**************************************/
int MatrixInv(int n, double a[], double b[])
{
	int i, j, k, l, u, v, is[100], js[100];   /* matrix dimension <= 10 */
	double d, p;

	if (n <= 0)
	{
		printf("Error dimension in MatrixInv!\n");
		exit(EXIT_FAILURE);
	}

	/* 将输入矩阵赋值给输出矩阵b，下面对b矩阵求逆，a矩阵不变 */
	for (i = 0; i<n; i++)
	{
		for (j = 0; j<n; j++)
		{
			b[i*n + j] = a[i*n + j];
		}
	}

	for (k = 0; k<n; k++)
	{
		d = 0.0;
		for (i = k; i<n; i++)   /* 查找右下角方阵中主元素的位置 */
		{
			for (j = k; j<n; j++)
			{
				l = n*i + j;
				p = fabs(b[l]);
				if (p>d)
				{
					d = p;
					is[k] = i;
					js[k] = j;
				}
			}
		}

		if (d<DBL_EPSILON)   /* 主元素接近于0，矩阵不可逆 */
		{
			printf("Divided by 0 in MatrixInv!\n");
			return 0;
		}

		if (is[k] != k)  /* 对主元素所在的行与右下角方阵的首行进行调换 */
		{
			for (j = 0; j<n; j++)
			{
				u = k*n + j;
				v = is[k] * n + j;
				p = b[u];
				b[u] = b[v];
				b[v] = p;
			}
		}

		if (js[k] != k)  /* 对主元素所在的列与右下角方阵的首列进行调换 */
		{
			for (i = 0; i<n; i++)
			{
				u = i*n + k;
				v = i*n + js[k];
				p = b[u];
				b[u] = b[v];
				b[v] = p;
			}
		}

		l = k*n + k;
		b[l] = 1.0 / b[l];  /* 初等行变换 */
		for (j = 0; j<n; j++)
		{
			if (j != k)
			{
				u = k*n + j;
				b[u] = b[u] * b[l];
			}
		}
		for (i = 0; i<n; i++)
		{
			if (i != k)
			{
				for (j = 0; j<n; j++)
				{
					if (j != k)
					{
						u = i*n + j;
						b[u] = b[u] - b[i*n + k] * b[k*n + j];
					}
				}
			}
		}
		for (i = 0; i<n; i++)
		{
			if (i != k)
			{
				u = i*n + k;
				b[u] = -b[u] * b[l];
			}
		}
	}

	for (k = n - 1; k >= 0; k--)  /* 将上面的行列调换重新恢复 */
	{
		if (js[k] != k)
		{
			for (j = 0; j<n; j++)
			{
				u = k*n + j;
				v = js[k] * n + j;
				p = b[u];
				b[u] = b[v];
				b[v] = p;
			}
		}
		if (is[k] != k)
		{
			for (i = 0; i<n; i++)
			{
				u = i*n + k;
				v = is[k] + i*n;
				p = b[u];
				b[u] = b[v];
				b[v] = p;
			}
		}
	}

	return 1;
}

/*************************************
MatrixTrans
目的：矩阵转置

参数：
a  输入矩阵
line_a  矩阵a行数
row_a   矩阵a列数
b       输出矩阵，b=trans(a)
line_b  矩阵b行数
row_b   矩阵b列数
返回值：1=正常，0=致命错误
**************************************/
//int MatrixTrans(int line_a, int row_a, int line_b, int row_b, double a[], double b[])
int MatrixTrans(int line_a, int row_a, double a[], double b[])
{
	int i, j;
	/*
	if (line_a != row_b || row_a != line_b)
	{
		cout << "Error dimension in MatrixTrans!" << endl;
		return 0;
	}
	else
	{
		for (i = 0; i < line_a; i++)
		{
			for (j = 0; j < row_a; j++)
			{
				b[j*line_a + i] = a[i*row_a + j];
			}
		}
		return 1;
	}
	*/
	for (i = 0; i < line_a; i++)
	{
		for (j = 0; j < row_a; j++)
		{
			b[j*line_a + i] = a[i*row_a + j];
		}
	}
	return 1;
}