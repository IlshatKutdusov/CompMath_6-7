#include <iostream>
#include "locale.h"
#include <math.h>
#include <iomanip>

# define M_PI 3.14159265358979323846

using namespace std;
const int rr = 10000; // Округление. Количество нулей в числе = количество знаков после запятой
const double x0 = 1.4, xn = 1.6; // Границы отрезка
const int n = 14; // Кол-во шагов
const double x[] = { x0, (x0 + xn) / 2.0, xn }; // xi (узловые точки)
double y[3]; // yi (узловые точки)
const double hx = (xn - x0) / n; // Для построения ИФЛ и АП
const double ht = 2.0 / n; // Для I и II ИФН
double a[3]; // коэффициенты P2
/*
* Содержание работы:
1. Изучить принципы построения интерполяционной формулы Лагранжа, I и II интерполяционных формул Ньютона и аппроксимационного полинома.
2. На конкретном примере усвоить порядок построения указанных полиномов с помощью ЭВМ.
3. Составить программу на любом языке программирования, реализующую процесс построения указанных полиномов второго порядка для системы из трех равноотстоящих узловых точек.
4. Сделать вывод о точности построения полиномов.
5. Составить отчет о проделанной работе.
*/

double f(double x) // Исходная функция
{
	///return round((8.8 / sqrt(M_PI)) * exp((2.0 * x + 3.0) / (x * x + 1.0)) * 1000) / 1000;
	//return 19 / M_PI * exp(pow(x + 1, 2));
	return x * x;
}
void uzlovye_tochki() // Рассчет узловых точек
{
	y[0] = f(x[0]);
	y[1] = f(x[1]);
	y[2] = f(x[2]);
	cout << "-------------------------------------------------" << endl;
	cout << "| xi |" << setw(8) << "x0 = " << x[0] << setw(4) << "|" << setw(8) << "x1 = " <<
		x[1] << setw(4) << "|" << setw(8) << "x2 = " << x[2] << setw(4) << "|" << endl;
	cout << "-------------------------------------------------" << endl;
	cout << "| yi |" << setw(9) << y[0] << setw(4) << "|" << setw(11) << y[1] << setw(4) <<
		"|" << setw(9) << y[2] << setw(4) << "|" << endl;
}
double L2(double _x) // Многочлен Лагранжа
{
	double a[3];
	a[0] = y[0] / ((x[0] - x[1]) * (x[0] - x[2]));
	a[1] = y[1] / ((x[1] - x[0]) * (x[1] - x[2]));
	a[2] = y[2] / ((x[2] - x[0]) * (x[2] - x[1]));
	double result = a[0] * (_x - x[1]) * (_x - x[2]) + a[1] * (_x - x[0]) * (_x - x[2]) +
		a[2] * (_x - x[0]) * (_x - x[1]);
	return round(result * rr) / rr;
}
double N1(double t) // I интерполяционная формула
{
	double dy[] = { y[1] - y[0], y[2] - y[1] };
	double result = y[0] + t * dy[0] + t * (t - 1) / 2 * (dy[1] - dy[0]);
	return round(result * rr) / rr;
}
double N2(double t) // II интерполяционная формула
{
	double dy[] = { y[1] - y[0], y[2] - y[1] };
	double result = y[2] + t * dy[1] + t * (t + 1) / 2 * (dy[1] - dy[0]);
	return round(result * rr) / rr;
}
void koefs_P2()
{

	double sumx[5], sumy, sumxy, sumxxy; // sumx[5], 5 - просто для удобства пользования нумерацией, 0 не используем
	double det[4]; // Определители
	sumx[1] = x[0] + x[1] + x[2];
	sumx[2] = (x[0] * x[0]) + (x[1] * x[1]) + (x[2] * x[2]);
	sumx[3] = (x[0] * x[0] * x[0]) + (x[1] * x[1] * x[1]) + (x[2] * x[2] * x[2]);
	sumx[4] = (x[0] * x[0] * x[0] * x[0]) + (x[1] * x[1] * x[1] * x[1]) + (x[2] * x[2] *
		x[2] * x[2]);
	sumy = y[0] + y[1] + y[2];
	sumxy = x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
	sumxxy = x[0] * x[0] * y[0] + x[1] * x[1] * y[1] + x[2] * x[2] * y[2];
	det[0] = 3 * sumx[2] * sumx[4] + sumx[1] * sumx[3] * sumx[2] + sumx[2] * sumx[1] *
		sumx[3] - sumx[2] * sumx[2] * sumx[2] - sumx[3] * sumx[3] * 3 - sumx[4] * sumx[1] * sumx[1];
	det[1] = sumy * sumx[2] * sumx[4] + sumx[1] * sumx[3] * sumxxy + sumx[2] * sumxy *
		sumx[3] - sumxxy * sumx[2] * sumx[2] - sumx[3] * sumx[3] * sumy - sumx[4] * sumxy * sumx[1];
	det[2] = 3 * sumxy * sumx[4] + sumy * sumx[3] * sumx[2] + sumx[2] * sumx[1] * sumxxy -
		sumx[2] * sumxy * sumx[2] - sumxxy * sumx[3] * 3 - sumx[4] * sumx[1] * sumy;
	det[3] = 3 * sumx[2] * sumxxy + sumx[1] * sumxy * sumx[2] + sumy * sumx[1] * sumx[3] -
		sumx[2] * sumx[2] * sumy - sumx[3] * sumxy * 3 - sumxxy * sumx[1] * sumx[1];
	a[0] = det[1] / det[0];
	a[1] = det[2] / det[0];
	a[2] = det[3] / det[0];
}
double P2(double _x) // Аппроксимационный многочлен
{
	double result = a[2] * _x * _x + a[1] * _x + a[0];
	return round(result * rr) / rr;
}
void tab_func()
{
	double _x = x[0], t1 = 0.0, t2 = -2.0;
	for (int i = 0; i <= n; i++) {
		_x = round(_x * 10000) / 10000;
		t1 = round(t1 * 10000) / 10000;
		t2 = round(t2 * 10000) / 10000;
		cout << setw(3) << i << setw(3) << "|" << setw(6) << _x << setw(2) << "|" <<

			setw(5) << t1 << setw(2)

			<< "|" << setw(5) << t2 << setw(3) << "|" << setw(8) << f(_x) << setw(2)
			<< "|" << setw(9) << L2(_x) << setw(2) << "|" << setw(13) <<

			round(abs(f(_x) - L2(_x)) * rr) / rr << setw(8)

			<< "|" << setw(9) << N1(t1) << setw(2) << "|" << setw(13) <<

			round(abs(f(_x) - N1(t1)) * rr) / rr << setw(8)

			<< "|" << setw(9) << N2(t2) << setw(2) << "|" << setw(13) <<

			round(abs(f(_x) - N2(t2)) * rr) / rr << setw(8)

			<< "|" << setw(9) << P2(_x) << setw(2) << "|" << setw(13) <<

			round(abs(f(_x) - P2(_x)) * rr) / rr << setw(8)

			<< "|" << endl;
		_x += hx;
		t1 += ht;
		t2 += ht;
	}
}
void show_table()
{
	cout << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------\n";
	cout << " i   |   xi   |  ti   |  t'i  |  f(xi)   |  L2(xi)  | |(f(xi) - L2(xi))| |  N1(xi)  | |(f(xi) - N1(xi))| |  N2(xi)  | |(f(xi) - N2(xi))| |  P2(xi) | |(f(xi) - P2(xi))| |\n";
	cout << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------\n";
}
int main()
{

	setlocale(0, "");
	cout << "Лабораторная работа по Вычислительной математике No6-7\n" << endl;
	cout << "Выполнили студенты группы 4310 Кутдусов Ильшат и Мулюков Роберт\n" <<
		endl;
	cout << "ИНТЕРПОЛЯЦИЯ И АППРОКСИМАЦИЯ ФУНКЦИЙ" << endl;
	cout << "Вариант No20:" << endl;
	cout << "(8.8/sqrt(pi))*e^((2*x+3)/(x^2+1)), x ∈ [0, 1]" << endl;
	uzlovye_tochki();
	koefs_P2();
	show_table();
	tab_func();
	system("pause");
	return 0;
}