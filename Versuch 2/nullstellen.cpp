#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <float.h>
#include <iostream>

using namespace std;
const double PI = 4.*atan(1.);
int fkt_nummer, d, N;
double x0;

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

double funktion(double x){
  switch (fkt_nummer) {
    case 1:
      return ln(x); /* ln oder log -> nachgucken */
    case 2:
      return: (x-2)*(x-2);
    case 3:
      return: cosh(x);
    case 4:
      return: sqrt(x);
  }
}

double euklid(double x){
  return (x0*x0 - 2*x0*x + x*x + funktion(x0)*funktion(x0) - 2*funktion(x0)*funktion(x) + funktion(x)*funktion(x) - d*d);
}

double euklid_ab(double x){
  return (- 2*x0 + 2*x - 2*funktion(x0)*funktion_ab1(x) + 2*funktion(x)*funktion_ab1(x) - d*d);
}

double newton(double xn, double(*f)(double), double(*f_ab)(double))
{
  double x;

  for (int i = 0; i < 20; i++)
  {
    x = xn - f(xn)/f_ab(xn);
    xn = x;
  }

}

int main()
{
  cout << "(1) ln(x)"<< endl;
  cout << "(2) (x-2)^2"<< endl;
  cout << "(3) cosh(x)"<< endl;
  cout << "(4) sqrt(x)"<< endl;
  cout << "Funktion durch Eingabe der jeweiligen Nummer wählen:" << endl;
  cin >> fkt_nummer;
  cout << "Abstand d eingeben:" << endl;
  cin >> d;
  cout << "Startwert x0 eingeben:" << endl;
  cin >> x0;
  cout << "Anzahl N eingeben:" << endl;
  cin >> N;

  y0 = funktion(x0);

	return 0;
}
