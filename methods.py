#include <iostream>

using namespace std;

#include <iostream>
#include <iomanip>

using namespace std;

double g(double* x)
{
	return 2 * x[0] + 3 * x[1] - 1;
}

double f1(double* x)
{
	return x[0] * x[0] + 8 * x[1] * x[1] - x[0] * x[1] + x[0];
}

double rk;
double (*f)(double* x);
int m;
int n;


double g_plus(double x)
{
	if (x > 0)
		return x;
	else return 0;
}


double func(double* x, int n)
{
	return f(x) + (rk / 2.0) * pow(g(x), 2);
}

double p(double* x)
{	
	return (rk / 2.0) * pow(g(x), 2);
}

void multiplyMatrixByVector(double** m, double* v, double* res, double n)
{
	for (int i = 0; i < n; i++)
		res[i] = 0;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			res[i] += m[i][j] * v[j];
}


void nablaF(double* x, double* res, double (*func)(double* x, int n), int n, double deltax)
{
	//res[0] = 2 * x[0] - x[1] + 1;
	//res[1] = 8 * x[1] - x[0];

	double* p1 = new double[n];
	double* p2 = new double[n];
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < n; i++)
		{
			if (i == j)
			{
				p1[i] = x[i] + deltax;
				p2[i] = x[i] - deltax;
			}
			else
			{
				p1[i] = x[i];
				p2[i] = x[i];
			}
		}
		res[j] = (func(p1, n) - func(p2, n)) / (2 * deltax);
	}
	delete[]p1;
	delete[]p2;
}

void hessian(double** hessian, double* x, double (*func)(double* x, int n), double deltax0, double deltax1, int n)
{
	/*	hessian[0][0] = 2;
		hessian[0][1] = -1;
		hessian[1][0] = -1;
		hessian[1][1] = 8;*/

	double* p1 = new double[n];
	double* p2 = new double[n];
	double* p3 = new double[n];
	double* p4 = new double[n];

	for (int i = 0; i < n; i++)
	{
		p1[i] = x[i] + deltax0;
		p2[i] = x[i] - deltax0;
		p3[i] = x[i] + deltax0;
		p4[i] = x[i] - deltax0;


		for (int j = 0; j < n; j++)
		{
			if (i == j)
			{
				p1[j] += deltax1;
				p2[j] += deltax1;
				p3[j] -= deltax1;
				p4[j] -= deltax1;
			}
			else
			{
				p1[j] = x[j] + deltax1;
				p2[j] = x[j] + deltax1;
				p3[j] = x[j] - deltax1;
				p4[j] = x[j] - deltax1;
			}
			for (int k = 0; k < n; k++)
			{
				if (k != i && k != j)
				{
					p1[k] = x[k];
					p2[k] = x[k];
					p3[k] = x[k];
					p4[k] = x[k];
				}
			}

			hessian[i][j] = (func(p1, n) - func(p2, n) - func(p3, n) + func(p4, n)) / (4 * deltax0 * deltax1);

			if (i == j)
			{
				p1[j] -= deltax1;
				p2[j] -= deltax1;
				p3[j] += deltax1;
				p4[j] += deltax1;
			}
		}

	}

	delete[]p1;
	delete[]p2;
	delete[]p3;
	delete[]p4;

}

void reversematrix(double** matrix, double** revmatr, int n)
{
	double buff = matrix[0][0];
	double m = 1.0 / (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]);
	revmatr[0][0] = matrix[1][1] * m;
	revmatr[0][1] = -matrix[0][1] * m;
	revmatr[1][0] = -matrix[1][0] * m;
	revmatr[1][1] = buff * m;
}

double* newton_ravson(double* x0, double (*func)(double* x, int n), double eps1, double eps2, int m, int n)
{
	int k = 0;
	bool cond = false;
	double* xk = new double[n];
	double* xkplus1 = new double[n];
	double* nablaxk = new double[n];
	double* dk = new double[n];
	double** hess = new double* [n];
	double** hessreverse = new double* [n];
	for (int i = 0; i < n; i++)
	{
		hess[i] = new double[n];
		hessreverse[i] = new double[n];
	}

	for (int i = 0; i < n; i++)
		xk[i] = x0[i];

	while (true)
	{
		nablaF(xk, nablaxk, func, n, 0.05);

		double normanablaxk = 0;
		for (int i = 0; i < n; i++)
			normanablaxk += pow(nablaxk[i], n);
		normanablaxk = sqrt(normanablaxk);

		if (normanablaxk <= eps1)
		{
			delete[]xkplus1;
			delete[]dk;
			delete[]hess;
			delete[]hessreverse;
			delete[]nablaxk;
			//cout << "Количество шагов = " << k << endl;
			return xk;
		}

		if (k >= m)
		{
			delete[]xkplus1;
			delete[]dk;
			delete[]hess;
			delete[]hessreverse;
			delete[]nablaxk;
			//cout << "Количество шагов = " << k << endl;
			return xk;
		}


		hessian(hess, xk, func, 0.05, 0.05, n);
		reversematrix(hess, hessreverse, n);

		double delta1 = 0, delta2 = 0;
		delta1 = hessreverse[0][0];
		delta2 = hessreverse[0][0] * hessreverse[1][1] - hessreverse[0][1] * hessreverse[1][0];

		if (delta1 > 0 && delta2 > 0)
		{
			multiplyMatrixByVector(hessreverse, nablaxk, dk, n);
		}
		else
		{
			for (int i = 0; i < n; i++)
				dk[i] = nablaxk[i];
		}


		double tk;
		double a = 0, b = 2, epsilon = 0.0000000000001, delta = 0.000000000000000000001;
		double* p1 = new double[n];
		double* p2 = new double[n];
		do
		{
			double x0 = (a + b - delta) / 2.0;
			double y0 = (a + b + delta) / 2.0;

			for (int i = 0; i < n; i++)
			{
				p1[i] = xk[i] - x0 * dk[i];
				p2[i] = xk[i] - y0 * dk[i];
			}

			if (func(p1, n) < func(p2, n))
				b = y0;
			else if (func(p1, n) > func(p2, n))
				a = x0;
			else if (func(p1, n) == func(p2, n))
			{
				a = x0;
				b = y0;
			}
		} while (abs(b - a) >= 2 * epsilon);
		delete[]p1;
		delete[]p2;

		tk = (a + b) / 2.0;
		for (int i = 0; i < n; i++)
			xkplus1[i] = xk[i] - tk * dk[i];

		//Диагностический вывод
		//cout << "f(xk)   = " << func(xk, n) << endl;
		//cout << "f(xk+1) = " << func(xkplus1, n) << endl;
		//

		double normaxkplus1minusxk = 0;
		for (int i = 0; i < n; i++)
			normaxkplus1minusxk += pow(xkplus1[i] - xk[i], 2);
		normaxkplus1minusxk = sqrt(normaxkplus1minusxk);

		if (normaxkplus1minusxk < eps2 && abs(func(xkplus1, n) - func(xk, n)) < eps2 && cond)
		{
			delete[]xk;
			delete[]dk;
			delete[]hess;
			delete[]hessreverse;
			delete[]nablaxk;
			//cout << "Количество шагов = " << k + 1 << endl;
			return xkplus1;
		}

		if (normaxkplus1minusxk < eps2 && abs(func(xkplus1, n) - func(xk, n)) < eps2)
		{
			cond = true;
		}
		else
			cond = false;

		for (int i = 0; i < n; i++)
			xk[i] = xkplus1[i];
		k += 1;
	}
}


double* strafniy_functs(double* x0, int n, double r0, double c, double eps)
{
	int k = 0;
	double* xk = new double[n];
	for (int i = 0; i < n; i++)
		xk[i] = x0[i];

	while (true)
	{
		double* xmin;
		xmin = newton_ravson(xk, func, 0.01, 0.01, 100, n);
		if (p(xmin) <= eps)
		{
			delete xk;
			cout << "Количество шагов: " << k << endl;
			return xmin;
		}
		else
		{
			rk = c * rk;
			for (int i = 0; i < n; i++)
				xk[i] = xmin[i];
		}
		delete[] xmin;
		k++;
	}

	cout << "Количество шагов: " << k;
	return 0;
}

int main()
{
	setlocale(LC_ALL, "rus");
	cout << fixed << setprecision(6);
	double* xmin;
	double x0[] = { 3, 1 };
	n = 2;
	m = 1;
	rk = 0.01;
	double c = 20;
	double eps = 0.01;
	f = f1;
	xmin = strafniy_functs(x0, n, rk, c, eps);
	for (int i = 0; i < n; i++)
		cout  << xmin[i] << " ";
	cout << "\nmin f(x) = " << f(xmin) << endl;
	delete[]xmin;

	return 0;
}
