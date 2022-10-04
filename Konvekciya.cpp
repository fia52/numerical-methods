#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <string>
#include <math.h>
#include <float.h>
#include <assert.h>
#include "discpp.h"

using namespace std;
bool RevMode = false;
const int jmax = 21, nmax = 100, maxex = 100;
const double u = 0.1, tmax = 5, alpha = 1.0e-05, L = 1.0, scale = 100.0, tol = 1.0e-05, s = 0.5, eps = 1.e-03, check = -87.0, PI = atan(1.) * 4.;
void DIFFextra(double[], double[], double);
double RMS(double[], double[], int);
void BTCS(double);
void Plot(int, double, double[], double[], double[], double, double, double);

int main()
{
	double cur = 0.8;
	BTCS(cur);
	//cin.get();
	//FTCS(td2, s2);
	//cin.get();
	//FTCS(td3, s3);
	//return EXIT_SUCCESS;
}

void BTCS(double cur)
{
	int i, j, jmap, ajm, n, m = 20;
	double delx, delt, t, miny, maxy, ystep, aa, bb, cc, td[jmax], d[3], x[jmax], te[jmax];
	jmap = jmax - 1; ajm = jmap;
	delx = L / ajm; delt = cur * delx / u;
	for (i = 0; i < jmax; i++) td[i] = 0.0; 
	n = 0; 
	t = 0;
	//
	for (j = 0; j < jmax; j++)
	{
		x[j] = delx * j;
		if (x[j] >= 0. && x[j] <= 0.1) td[j] = sin(10 * PI * x[j]);
		if (x[j] > 0.1 && x[j] <= L) td[j] = 0.;
	}
	//
	aa = 1.; bb = -cur; cc = cur;
	do
	{
		td[0] = 0.; td[jmax - 1] = 0.;
		d[1] = td[0];
		for (j = 1; j < jmap; j++)
		{
			d[2] = d[1];
			d[1] = td[j];
			d[0] = td[j];
			td[j] = aa * d[0] + bb * d[1] + cc * d[2];
		}
		t = t + delt; n = n + 1;
		if (n >= nmax) break;
		cout << fixed << setprecision(2);
		cout << "t= " << setw(7) << t << " ";
		for (i = 0; i < jmax; i++) cout << " " << setw(7) << td[i];
		cout << endl;
	} while (t < tmax);
	DIFFextra(te, x, t);
	cout << endl << endl;
	cout << "t= " << setw(7) << t << " ";
	for (i = 0; i < jmax; i++) cout << setw(8) << te[i];
	cout << endl << "x= " << setw(8) << " ";
	for (i = 0; i < jmax; i++) cout << setw(8) << x[i];
	cout << endl << endl;
	cout << setprecision(4) << " Error estimate, RMS = " << RMS(te, td, ajm)
		<< "; Delta x = " << delx << "; Delta t = " << delt << "; Parameter cur =  " << cur << endl;
	miny = td[0];
	for (i = 0; i < m; ++i)miny = fmin(miny, td[i]);
	for (i = 0; i < m; ++i)miny = fmin(miny, te[i]);
	maxy = td[0];
	for (i = 0; i < m; ++i)maxy = fmax(maxy, td[i]);
	for (i = 0; i < m; ++i)maxy = fmax(maxy, te[i]);
	ystep = (maxy - miny) / jmax;
	Plot(jmax, delx, x, td, te, ystep, miny, maxy);
}
void DIFFextra(double te[], double x[], double time)
{
	int j;
	for (j = 0; j < jmax; j++)
	{
		if (x[j] >= 0. && x[j] <= u * time) te[j] = 0.;
		if (x[j] > u * time && x[j] <= u * time + 0.1) te[j] = sin(10 * PI * (x[j] - u * time));
		if (x[j] > u * time + 0.1 && x[j] <= 1.0) te[j] = 0.;
	}
}
double RMS(double te[], double td[], int ajm)
{
	int j;
	double rms, dmp, avs, sum = 0.0;
	for (j = 0; j < jmax; ++j)
	{
		dmp = te[j] - td[j];
		sum = sum + dmp * dmp;
	}
	//estimate the errors
	avs = sum / (1.0 + ajm);
	rms = sqrt(avs);
	return rms;
}
void Plot(int m, double xstep, double xout[], double yout[], double textra[], double ystep, double miny, double maxy)
{
	Dislin g;
	int n = m - 1;
	char cdev[] = "XWIN";
	g.metafl(cdev);
	g.setpag("da41");
	if (RevMode)g.scrmod("revers");
	g.disini();
	g.pagera();
	g.hwfont();
	g.axspos(450, 1800);
	g.axslen(2200, 1400);
	g.name("X-axis", "x");
	g.name("Y-axis", "y");
	g.labdig(-1, "x");
	g.ticks(10, "xy");
	g.titlin("extra(blue), our(yellow)", 2);
	g.graf(xout[0], xout[n], xout[0], xstep*20, miny, maxy, textra[0], ystep*2.6);
	g.title();
	g.color("yellow");
	g.curve(xout, yout, m);
	g.color("blue");
	g.curve(xout, textra, m);
	g.dash();
	g.xaxgit();
	g.disfin();
}
