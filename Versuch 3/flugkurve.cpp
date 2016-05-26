// P3_2015.cpp : Definiert den Einstiegspunkt für die Konsolenanwendung.
//

// #include "stdafx.h" remove since it is not found
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

using namespace std;

int Romberg_3(double(*f)(double), double a, double b, double tol, int L_in, int &L_out, int &fcnt, double &Q, double &E);
int spline(int n, double *x, double *y, double marg_0, double marg_n, int marg_cond, double *b, double *c, double *d);
double spval(int n, double x0, double *a, double *b, double *c, double *d, double *x, double *ausg);
void datenlesen(vector<double> &t, vector<double>&x, vector<double> &y, vector<double> &z);

int n;
vector<double> t, x, y, z;
vector<double> xb, xc, xd, yb, yc, yd, zb, zc, zd;

double laenge(double x0){
	double sx[3], sy[3], sz[3];
	spval(n, x0, &x[0], &xb[0], &xc[0], &xd[0], &t[0], sx);
	spval(n, x0, &y[0], &yb[0], &yc[0], &yd[0], &t[0], sy);
	spval(n, x0, &z[0], &zb[0], &zc[0], &zd[0], &t[0], sz);

	return sqrt(sx[0]*sx[0] + sy[0]*sy[0] + sz[0]*sz[0]);
}

int main(void)
{
	fstream file1;
	file1.open("Flug_testen.txt", ios::out);
	file1 << "Heute ist schönes Wetter" << endl;
	file1.close();

	datenlesen(t, x, y, z);
	// an dieser Stelle das Hauptprogramm einfuegen !!

	n = t.size();
	xb.resize(n);
	xc.resize(n);
	xd.resize(n);
	yb.resize(n);
	yc.resize(n);
	yd.resize(n);
	zb.resize(n);
	zc.resize(n);
	zd.resize(n);

	spline(n, &t[0], &x[0], 0, 0, 2, &xb[0], &xc[0], &xd[0]);
	spline(n, &t[0], &y[0], 0, 0, 2, &yb[0], &yc[0], &yd[0]);
	spline(n, &t[0], &z[0], 0, 0, 2, &zb[0], &zc[0], &zd[0]);

	double ta, tb;
	cout << "ta eingeben: ";
	cin >> ta;
	ta = ta*3600;
	cout << "tb eingeben: ";
	cin >> tb;
	tb = tb*3600;

	if(ta >= tb){
		cout << "Fehler: ta muss groesser als tb sein! ta und tb wurden auf die globalen Intervallgrenzen gesetzt" << endl;
		ta = t[0];
		tb = t[n-1];
	}

	int L_out, fcnt;
	double Q, E;

	Romberg_3(laenge, ta, tb, 1e-14, 20, L_out, fcnt, Q, E);

}

/* ----------------------------- Modul datenlesen ------------------------- */
void datenlesen(vector<double> &t, vector<double>&x, vector<double> &y, vector<double> &z)
{
	double tb;
	fstream file;
	char inputFilename[] = "Flug_txyz.txt";
	file.open(inputFilename, ios::in);
	if (!file) {
		cerr << "Eingabe-File " << inputFilename << " kann nicht geöffnet werden!" << endl;
		cout << " Gib irgendeine Zahl ein (Modul:datenlesen)" << endl;
		cin >> tb;
		exit(1);
	}

	double tt, tx, ty, tz;
	int helpcount = 0;
	file >> tt >> tx >> ty >> tz;
	while (!file.eof())
	{
		helpcount++;
		t.push_back(tt);
		x.push_back(tx);
		y.push_back(ty);
		z.push_back(tz);
		file >> tt >> tx >> ty >> tz;
	}
	file.close();
	//cout << "Laenge: " << t.size() << endl;
}
/* ----------------------------- Ende datenlesen ------------------------- */


/* ----------------------------- Modul Tridiagonal --------------------------- */
int trdiag(int n, double *lower, double *diag, double *upper, double *b, int rep)
{
	register int i;

	if (n < 2) return 1;
	if (rep == 0)
	{
		for (i = 1; i<n; i++)
		{
			if (fabs(diag[i - 1]) < DBL_EPSILON) return 2;
			lower[i] /= diag[i - 1];
			diag[i] -= lower[i] * upper[i - 1];
		}
	}
	if (fabs(diag[n - 1]) < DBL_EPSILON) return 2;
	for (i = 1; i<n; i++)
		b[i] -= lower[i] * b[i - 1];
	b[n - 1] /= diag[n - 1];
	for (i = n - 2; i >= 0; i--)
		b[i] = (b[i] - upper[i] * b[i + 1]) / diag[i];
	return 0;
}
/* ----------------------------- Ende Tridiagonal ---------------------------- */

/* ----------------------------- Modul Spline -------------------------------- */

int spline(int n, double *x, double *y, double marg_0, double marg_n, int marg_cond, double *b, double *c, double *d)
/********************************************************************************
* spline berechnet zu den vorgegebenen Wertepaaren                              *
*      (x[i], y[i], i=0(1)n-1)                                                  *
* die Koeffizienten eines kubischen Splines.                                    *
* Die Art der Randbedingung wird durch die Variable marg_cond                   *
* festgelegt.                                                                   *
*===============================================================================*
* Eingabeparameter:                                                             *
* -----------------                                                             *
* Name		Typ/Laenge              Bedeutung									*
* ------------------------------------------------------------------------------*
* n			int/---		Anzahl der Stuetzstellen; n>2                           *
* x			double/[n]	x-Koordinaten der Wertepaare                            *
* y			double/[n]	y-Koordinaten der Wertepaare                            *
* marg_0	double/---	Randbedingung in x[0]                                   *
* marg_n	double/---	Randbedingung in x[n-1]                                 *
* marg_cond int/---		= 0 : not-a-knot-Bedingung                              *
*						       marg_0,marg_n ohne Bedeutung                     *
*						= 1 : marg_0, marg_n sind 1. Ableitungen                *
*						= 2 : marg_0, marg_n sind 2. Ableitungen                *
*						= 3 : marg_0, marg_n sind 3. Ableitungen				*
*                                                                               *
* Ausgabeparameter:                                                             *
* -----------------                                                             *
* Name		Typ/Laenge			Bedeutung                                       *
* ------------------------------------------------------------------------------*
* b			double/[n]	Splinekoeffizienten:			                        *
* c			double/[n]	s(x)=a[i]+b[i]*t+c[i]*t^2+d[i]*t^3			            *
* d	        double/[n]	t=x-x[i]; x aus [ x[i], x[i+1] ]			            *
*						(a entspricht y)					                    *
*                       Die Splinekoeffizienten brauchten nur mit [n-1]         *
*						vereinbart werden. Es wird jedoch bei der Rechnung      *
*						[n] als Zwischenspeicher benoetigt.                     *
*																	        	*
* Rueckgabewert:                                                                *
* -------------                                                                 *
* = -i : Monotoniefehler: x[i-1] >= x[i]                                        *
* =  0 : kein Fehler                                                            *
* =  1 : falscher Wert fuer marg_cond                                           *
* =  2 : n < 3                                                                  *
* =  3 : es steht nicht genuegend Speicher fuer die notwendigen Hilfsfelder     *
*        zur Verfuegung                                                         *
* >  3 : Fehler in trdiag                                                       *
********************************************************************************/
{
	int i, error;
	double *h, *a;

	if (n < 3) return (2);
	/* es wird Speicherplatz fuer das Hilfsfeld zur Verfuegung gestellt */
	a = (double *)malloc(n * sizeof (double));
	if (a == NULL) return 3;
	h = (double *)malloc((--n) * sizeof (double));
	if (h == NULL) { free(a); return 3; }

	for (i = 0; i <= n - 1; i++)
	{
		h[i] = x[i + 1] - x[i];
		if (h[i] <= 0.) return (-i);  /* Ueberpruefung der Monotonie */
	}
	/* Aufstellen des Gleichungssystems */
	for (i = 0; i <= n - 2; i++)
	{
		a[i] = 3. * ((y[i + 2] - y[i + 1]) / h[i + 1] - (y[i + 1] - y[i]) / h[i]);
		b[i] = h[i];
		c[i] = h[i + 1];
		d[i] = 2. * (h[i] + h[i + 1]);
	}
	/* In Abhaengigkeit von der Randbedingung werden Werte neu bestezt */
	switch (marg_cond)
	{
	case 0:
	{
			  if (n == 2)
			  {
				  a[0] /= 3.;
				  d[0] *= 0.5;
			  }
			  else
			  {
				  a[0] *= h[1] / (h[0] + h[1]);
				  a[n - 2] *= h[n - 2] / (h[n - 1] + h[n - 2]);
				  d[0] -= h[0];
				  d[n - 2] -= h[n - 1];
				  c[0] -= h[0];
				  b[n - 2] -= h[n - 1];
			  }
			  break;
	}
	case 1:
	{
			  a[0] -= 1.5 * ((y[1] - y[0]) / h[0] - marg_0);
			  a[n - 2] -= 1.5 * (marg_n - (y[n] - y[n - 1]) / h[n - 1]);
			  d[0] -= h[0] * 0.5;
			  d[n - 2] -= h[n - 1] * 0.5;
			  break;
	}
	case 2:
	{
			  a[0] -= h[0] * marg_0 * 0.5;
			  a[n - 2] -= h[n - 1] * marg_n * 0.5;
			  break;
	}
	case 3:
	{
			  a[0] += marg_0 * h[0] * h[0] * 0.5;
			  a[n - 2] -= marg_n * h[n - 1] * h[n - 1] * 0.5;
			  d[0] += h[0];
			  d[n - 2] += h[n - 1];
			  break;
	}
	default: {free(a); free(h); return 1; }
	}
	/* Berechnen der Koeffizienten */
	switch (n - 1)
	{
	case 1:
	{
			  c[1] = a[0] / d[0];
			  break;
	}
	default:
	{
			   error = trdiag(n - 1, b, d, c, a, 0);
			   if (error != 0) { free(a); free(h); return error + 3; }

			   for (i = 0; i <= n - 2; i++)  /* Ueberschreiben des Loesungs- */
				   c[i + 1] = a[i];     /* vektors auf c                */
	}
	}
	/* In Abhaengigkeit von der Randbedingung wird der erste und der letzte Wert von c neu besetzt. */
	switch (marg_cond)
	{
	case 0:
	{
			  if (n == 2)
				  c[0] = c[2] = c[1];
			  else
			  {
				  c[0] = c[1] + h[0] * (c[1] - c[2]) / h[1];
				  c[n] = c[n - 1] + h[n - 1] * (c[n - 1] - c[n - 2]) / h[n - 2];
			  }
			  break;
	}
	case 1:
	{
			  c[0] = 1.5 * ((y[1] - y[0]) / h[0] - marg_0);
			  c[0] = (c[0] - c[1] * h[0] * 0.5) / h[0];
			  c[n] = -1.5 * ((y[n] - y[n - 1]) / h[n - 1] - marg_n);
			  c[n] = (c[n] - c[n - 1] * h[n - 1] * 0.5) / h[n - 1];
			  break;
	}
	case 2:
	{
			  c[0] = marg_0 * 0.5;
			  c[n] = marg_n * 0.5;
			  break;
	}
	case 3:
	{
			  c[0] = c[1] - marg_0 * h[0] * 0.5;
			  c[n] = c[n - 1] + marg_n * h[n - 1] * 0.5;
	}
	}
	for (i = 0; i <= n - 1; i++)
	{
		b[i] = (y[i + 1] - y[i]) / h[i] - h[i] * (c[i + 1] + 2. * c[i]) / 3.;
		d[i] = (c[i + 1] - c[i]) / (3. * h[i]);
	}
	free(h);
	free(a);
	return 0;
}
/* ----------------------------- Ende Spline --------------------------------- */

/* ----------------------------- Modul SPVAL --------------------------------- */
double spval(int n, double x0, double *a, double *b, double *c, double *d, double *x, double *ausg)
/********************************************************************************
* spval berechnet Funktionswerte eines kubischen Splines und                    *
* seine nicht-trivialen Ableitungen.                                            *
*===============================================================================*
* Eingabeparameter:                                                             *
* -----------------                                                             *
* Name		Typ/Laenge		Bedeutung										    *
* ------------------------------------------------------------------------------*
* n		int/---			Anzahl der Stuetzstellen						        *
* x0	double/---		Stelle, an der der Funktionswert berechnet wird			*
* a		double/[n-1]	)														*
* b		double/[n-1]	) Splinekoeffizienten                                   *
* c		double/[n-1]	)														*
* d		double/[n-1]	)														*
* x		double/[n-1]	Stuetzstellen											*
*                                                                               *
* Ausgabeparameter:                                                             *
* -----------------                                                             *
* Name		Typ/Laenge	Bedeutung								                *
* ------------------------------------------------------------------------------*
* ausg		double/[3]	Ableitungen an der Stelle x0					        *
*				ausg[0] enthaelt die 1. Ableitung								*
*				ausg[1] enthaelt die 2. Ableitung							    *
*				ausg[2] enthaelt die 3. Ableitung							    *
*																				*
* Rueckgabewert:                                                                *
* -------------                                                                 *
* Funktionswert der Splinefunktion an der Stelle                                *
*               x0,     falls x[0]<=x0<=x[n-1]                                  *
*               x[0],   falls x0<x[0]                                           *
*               x[n-1], falls x[n-1]<x0                                         *
********************************************************************************/
{
	int i, k, m;
	//cout << "Hier ist spval..." << endl;
	//cout << "n = " << n << endl;
	if (x0 < x[0]) x0 = x[0];
	if (x0 > x[n - 1]) x0 = x[n - 1];

	i = 0, k = n;
	while (m = (i + k) >> 1, m != i)	/* das Intervall, in dem x0 */
	{
		/* liegt, wird gesucht.     */
		if (x0 < x[m]) k = m;
		else i = m;
	}
	if (i >= n - 1) i = n - 2;
	x0 -= x[i];
	//cout <<"in SPVAL: i = ["<<i<<"]"<<endl;

	ausg[0] = (3. * d[i] * x0 + 2. * c[i]) * x0 + b[i];
	ausg[1] = 6. * d[i] * x0 + 2. * c[i];
	ausg[2] = 6. * d[i];
	return (((d[i] * x0 + c[i]) * x0 + b[i]) * x0 + a[i]);
}

/* ----------------------------- Ende SPVAL ---------------------------------- */

/* ----------------------------- Modul Romberg_3 ---------------------------------- */
int Romberg_3(double(*f)(double), double a, double b, double tol, int L_in, int &L_out, int &fcnt, double &Q, double &E){
	/* N�herungsweise Berechnung des Integrals �ber f von a bis b: I(f;a,b)
	% mit dem Romberg-Schema mit max. L_in Stufen.
	%
	% Eingabe:
	% f    : Funktion der Form f(x); Integrand
	% a    : reelle Zahl; untere Integrationsgrenze
	% b    : reelle Zahl; obere Integrationsgrenze
	% tol  : reelle Zahl; Genauigkeitstoleranz: |I(f;a,b)-Q|<=tol*(1+|Q|)
	% L_in : Integer; max. Anzahl Stufen im Romberg-Schema
	% Ausgabe:
	% rc    : Integer; rc=0; Q wurde in gew�nschter Genauigkeit berechnet.
	%                  rc=1; Q ist nicht genau genug. Es w�rden evtl. mehr Stufen (als L_in) ben�tigt.
	%                  rc=2; Eingabefehler: tol <=0
	%
	% L_out : Integer; Anzahl genutzter Stufen
	% fcnt  : Integer; Anzahl Funktionsauswertungen des Integranden
	% Q     : reelle Zahl; N�herungswert f�r das Integral I(f;a,b)
	% E     : reelle Zahl; Fehlersch�tzung zu Q
	%------- R. Reuter ----- 20.05.2015 --------------------------------------*/
	double Q1 = 0., h = 0., ah = 0.;
	int n = 10; // Anzahl initialer Teilintervalle f�r ST-Formel
	int rc = 1; // Return Code
	double Fak = 1.0;
	int Lmax = L_in; // max. Anzahl Stufen
	int L = 0;
	double *QQ = new double[Lmax + 1]; // Array zum Speichern der jeweils aktuellen Zeile
	// des Romberg-Schemas
	fcnt = 0;
	if (tol <= 0.){
		rc = 2;
		return rc;
	}

	h = (b - a) / n;
	Q1 = 0.5 * (f(a) + f(b));
	fcnt = fcnt + 2;
	for (int k = 1; k<n; k++)
		Q1 = Q1 + f(a + k*h);
	fcnt = fcnt + n - 1;
	Q1 = h*Q1; // ST-Formel zu h
	QQ[1] = Q1;
	for (L = 2; L <= Lmax; L++){
		Q = 0.0;
		ah = a + 0.5*h;

		for (int k = 0; k<n; k++)
			Q = Q + f(ah + k*h);
		fcnt = fcnt + n;
		Q = 0.5*(QQ[1] + h*Q); // ST-Formel zu h/2
		Fak = 1.0;

		for (int LL = 1; LL<L; LL++){
			Fak = Fak*4.;
			Q1 = QQ[LL];
			QQ[LL] = Q;
			E = (Q - Q1) / (Fak - 1.0);

			Q = Q + E;
		}

		QQ[L] = Q;

		if (fabs(E) <= tol*(1 + fabs(Q))){
			rc = 0;
			break;
		}

		h = 0.5*h;
		n = 2 * n;
	}
	L_out = L;
	return rc;
}
/* ----------------------------- Ende Romberg_3 ---------------------------------- */
