#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <float.h>
#include <iostream>

using namespace std;
const double PI = 4.*atan(1.);
int fkt_nummer, d, N, fcnt;
double x0, xk, bogenstuecklaenge;

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

double f(double x){
  switch (fkt_nummer) {
    case 1:
      return log(x); /* ln oder log -> nachgucken */
    case 2:
      return (x-2)*(x-2);
    case 3:
      return cosh(x);
    case 4:
      return sqrt(x);
  }
}

double d1f(double x){
  switch (fkt_nummer) {
    case 1:
      return 1/x;
    case 2:
      return 2*x - 4;
    case 3:
      return sinh(x);
    case 4:
      return 1/(2*sqrt(x));
  }
}

double euklid(double x){
  return (x0*x0 - 2*x0*x + x*x + f(x0)*f(x0) - 2*f(x0)*f(x) + f(x)*f(x) - d*d);
}

double euklid_ab(double x){
  return (- 2*x0 + 2*x - 2*f(x0)*d1f(x) + 2*f(x)*d1f(x));
}

double bogenlaenge_integrant(double x){
  return sqrt(1 + d1f(x)*d1f(x));
}

double newton(double xn, double(*f)(double), double(*f_ab)(double)){
  double x;

  for (int i = 0; i < 100; i++)
  {
    x = xn - f(xn)/f_ab(xn);
    xn = x;
    if (abs(xn - x) <= 10e-10){
    	return x;
    }
  }

}

double romberg_newton(double x){
	/*
		Gibt (Q-gewünschte Bogenstuecklänge) zurück, um dies in newton zu verwenden.
	*/
	double tol = 10e-14, Q, rc;
	int L_out;

	rc = Romberg(bogenlaenge_integrant, xk, x, tol, 20, L_out, fcnt, Q);
	return Q - bogenstuecklaenge;
}

int main()
{
	// 2.1 a)
  double y0, xn, yn;
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

	// 2.1 b)
  y0 = f(x0);
	cout << "(x0,y0): (" << x0 << "," << y0 << ")" << endl;

  xn = newton(x0 + d, euklid, euklid_ab);
  yn = f(xn);

	cout << "(xn,yn): (" << xn << "," << yn << ")" << endl;

  /*
  cout << xn << endl;
  cout << yn << endl;
  cout << "d:" << sqrt((x0-xn)*(x0-xn) + (y0-yn)*(y0-yn)) << endl; // Testen ob d richtig ist
  */

  xk = x0;
  int L_out;
  double rc = Romberg(bogenlaenge_integrant, xk, xn, 10e-14, 20, L_out, fcnt, bogenstuecklaenge);
	cout << "Bogenlänge:" << bogenstuecklaenge << endl;

  bogenstuecklaenge = bogenstuecklaenge/N;
  cout << "Bogenstücklänge:" << bogenstuecklaenge << endl;

	// 2.1 c)
	double X[N+1], Y[N+1];
	X[0] = x0;
	Y[0] = y0;

	for (int i = 0; i < N; i++) {
		// x_quer berechnen mit Pythagoras und Steigung
		double x_quer = X[i] + bogenstuecklaenge/sqrt(1 + d1f(X[i])*d1f(X[i]));
		double y_quer = f(x_quer);

		// x_dach berechnen mit Pythagoras und Steigung
		double m = (f(x_quer) - f(X[i])) / (x_quer - X[i]);
		double x_dach = X[i] + bogenstuecklaenge/sqrt(1 + m*m); // evtl noch falsch!

		X[i+1] = newton(x_dach, romberg_newton, bogenlaenge_integrant);
		Y[i+1] = f(X[i+1]);
		xk = X[i + 1]; // für romberg_newton setzten
	}

	// 2.1 e)
	double px, py;
	px = (x0 + X[N]) / 2;
	py = (y0 + Y[N]) / 2;

	cout << "Mitte: (" << px << "," << py << ")" << endl;
	cout << "(X[N],Y[N]): (" << X[N] << "," << Y[N] << ")" << endl;

	// 2.1 f)
	double Winkel[N+1], Radius[N+1];

	double a_x = x0 - px;
	double a_y = y0 - py;
	for (int i = 0; i <= N; i++) {
		Radius[i] = sqrt((px - X[i])*(px - X[i]) + (py-Y[i])*(py-Y[i]));
		double b_x = X[i] - px;
		double b_y = X[i] - py;
		Winkel[i] = acos((a_x * b_x + a_y * b_y) / (Radius[i]*d/2)) * (180/PI);
	}


	return 0;
}
