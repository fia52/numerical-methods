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
const double one = 1.0, six = 6.0, tmax = 2999.0, alpha = 1.0e-05, L = 1.0, scale = 100.0, tol = 1.0e-05, s = 0.5, eps = 1.e-03, check = -87.0;
void DIFFextra(double[], double[], double);
double RMS(double[], double[], int);
void FTCS(double[], double);
void Plot(int, double[], double[], double[], double, double, double);

int main()
{
	double td1[jmax] = { NAN }, td2[jmax] = {NAN}, td3[jmax] = { NAN };;
	double s1 = 0.6, s2 = 0.1, s3 = 1./6.;
	FTCS(td1, s1);
	cin.get();
	FTCS(td2, s2);
	cin.get();
	FTCS(td3, s3);
	return EXIT_SUCCESS;
}

void FTCS(double td[], double s)
{
	int i, j, jmap, ajm, n, m = 20;
	double delx, delt, t, miny, maxy, ystep;
	jmap = jmax - 1; ajm = jmap;
	delx = L / ajm; delt = delx * delx * s / alpha;
	double tn[jmax], dum[jmax], x[jmax], te[jmax];
	for (i = 0; i < jmax; i++) tn[i] = 0.0;
	n = 0; t = 0;
	//
	for (j = 0; j < jmax; j++) x[j] = delx * j;
	//
	do
	{
		tn[0] = 1.0; tn[jmax - 1] = 1.0;
		if (t < tol)
		{
			tn[0] = 0.5;
			tn[jmax - 1] = 0.5;
		}
		td[0] = scale * tn[0]; td[jmax - 1] = scale * tn[jmax - 1];
		for (j = 1; j < jmap; j++) dum[j] = (1. - 2. * s) * tn[j] + s * (tn[j - 1] + tn[j + 1]);
		for (j = 1; j < jmap; j++)tn[j] = dum[j];
		for (j = 1; j < jmap; j++)td[j] = scale * tn[j];
		t = t + delt; n = n + 1;
		if (n >= nmax) break;
		cout << fixed << setprecision(2);
		cout << "t= " << setw(7) << t << " ";
		for (i = 0; i < jmax; i++) cout << " " << setw(7) << td[i];
		cout << endl;
	} 
	while (t < tmax);
	DIFFextra(te, x, t);
	cout << endl << endl;
	cout << "t= " << setw(7) << t << " ";
	for (i = 0; i < jmax; i++) cout << setw(8) << te[i];
	cout << endl << "x= " << setw(8) << " ";
	for (i = 0; i < jmax; i++) cout << setw(8) << x[i];
	cout << endl << endl;
	cout << setprecision(4) << " Error estimate, RMS = " << RMS(te, td, ajm)
		<< "; Delta x = " << delx << "; Delta t = " << delt << "; Parameter S =  " << s << endl;
	miny = td[0];
	for (i = 0; i < m; ++i)miny = fmin(miny, td[i]);
	for (i = 0; i < m; ++i)miny = fmin(miny, te[i]);
	maxy = td[0];
	for (i = 0; i < m; ++i)maxy = fmax(maxy, td[i]);
	for (i = 0; i < m; ++i)maxy = fmax(maxy, te[i]);
	ystep = (maxy - miny) / 20.0;
	Plot(m, x, td, te, ystep, miny, maxy);
}
void DIFFextra(double te[], double x[], double time)
{
	int j, m;
	double dxm, dam, dtm, PI = atan(1.) * 4.;
	for (j = 0; j < jmax; j++)
	{
		te[j] = 100.;
		for (m = 0; m < maxex - 1; m++)
		{
			dam = (2. * (m + 1) - 1.);
			dxm = dam * PI * x[j];
			dtm = -alpha * dam * dam * PI * PI * time;
			if (dtm < check) dtm = check;
			dtm = exp(dtm);
			if (dtm < eps) break;
			te[j] = te[j] - 400.0 / dam / PI * sin(dam * PI * x[j])*dtm; 
		}
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
void Plot(int m, double xout[], double yout[], double textra[], double ystep, double miny, double maxy)
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
	g.graf(xout[0], xout[n], xout[0], 1, miny, maxy, textra[0], ystep);
	g.title();
	g.color("yellow");
	g.curve(xout, yout, m);
	g.color("blue");
	g.curve(xout, textra, m);
	g.dash();
	g.xaxgit();
	g.disfin();
}
