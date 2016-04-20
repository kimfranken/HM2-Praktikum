// Romberg_Test.cpp : Definiert den Einstiegspunkt f�r die Konsolenanwendung.
//

//#include "stdafx.h"
#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <float.h>
#include <iostream>

using namespace std;
const double PI = 4.*atan(1.);


int Romberg(double(*f)(double), double a, double b, double tol, int L_in, int &L_out, int &fcnt, double &Q);

double mein_sin(double x){
	double tmp;
	tmp = sin(x);
	return tmp;//*tmp;
}


int main()
{
	double Q = 0.;
	double a = 0., b = PI;
	int L_in = 10, L_out = 0, rc = 0, k, fcnt;
	double I = 2.0;//PI/2;
	double eps = 1e-6;
	for (k = 0; k<11; k++)
	{
		eps = eps / 10.;
		rc = Romberg(mein_sin, a, b, eps, L_in, L_out, fcnt, Q);
		if (rc == 2){
			cout << "Fehler in Romberg aufgetreten. Fehler Nr. 2"<<endl;
			//exit(1);
		}
		else if (rc == 1)
			cout << "Fehler in Romberg aufgetreten. Fehler Nr. 1. Max. Stufe erreicht!"<<endl;
		else if (rc == 0)
			cout << " I-Q = " << I - Q << " zu eps = " << eps << " L_out = " << L_out << " fcnt = " << fcnt<<endl;
	}
	cout << "Irgeneine ganze Zahl eingeben!" << endl;
	int NIX;
	cin >> NIX;
	return 0;
}

int Romberg(double(*f)(double), double a, double b, double tol, int L_in, int &L_out, int &fcnt, double &Q){
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
	% Q     : relle Zahl; N�herungswert f�r das Integral I(f;a,b)
	% L_out : Integer; Anzahl genutzter Stufen
	%
	%------- R. Reuter ----- 24.03.2016 --------------------------------------*/
	double Q1 = 0., E = 0., h = 0., ah = 0.;
	int n = 10; // Anzahl initialer Teilintervalle f�r ST-Formel
	int rc = 1; // Return Code
	double Fak = 1.0;
	int Lmax = L_in; // max. Anzahl Stufen
	int L = 0;
	double *QQ = new double[Lmax + 1]; // Array zum Speichern der jeweils aktuellen Zeile
	// des Romberg-Schemas

	if (tol <= 0.){
		rc = 2;
		return rc;
	}

	h = (b - a) / n;
	Q1 = 0.5 * (f(a) + f(b));
	fcnt = 2;
	for (int k = 1; k<n; k++)
		Q1 = Q1 + f(a + k*h);
	fcnt = fcnt + (n - 1);
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
