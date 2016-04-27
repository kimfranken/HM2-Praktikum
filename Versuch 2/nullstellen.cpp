#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <float.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#define WIDTH 20
#define PREC 15

using namespace std;
const double PI = 4.*atan(1.);
int f_nummer, d, N, fcnt, its;
string f_name = "error";
double x0, xk, bogenlaenge, bogenstuecklaenge;
double tol = 1e-14;

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
  switch (f_nummer) {
    case 1:
      return log(x);
    case 2:
      return (x-2)*(x-2);
    case 3:
      return cosh(x);
    case 4:
      return sqrt(x);
  }
}

double d1f(double x){
  switch (f_nummer) {
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

double newton(double xn, double(*f)(double), double(*d1f)(double)){
  double x;

  for (int i = 1; i < 300; i++)
  {
    x = xn - f(xn)/d1f(xn);
    if (abs(xn - x) <= 10e-10){
			its = i;
    	return x;
    }
		xn = x;

		// evtl Error einbauen, wenn 300 erreicht.
  }
}

double romberg_newton(double x){
	/*
		Gibt (Q-gewünschte Bogenstuecklänge) zurück, um dies in newton zu verwenden.
	*/
	double Q, rc;
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
  cout << "(4) sqrt(x)"<< endl << endl;
  cout << "Funktion durch Eingabe der jeweiligen Nummer wählen: ";
  cin >> f_nummer;
  switch (f_nummer) {
    case 1:
      f_name = "ln(x)";
			break;
    case 2:
      f_name = "(x-2)^2";
			break;
    case 3:
      f_name = "cosh(x)";
			break;
    case 4:
      f_name = "sqrt(x)";
			break;
		case 0:
			f_name = "ln(x)";
			break;
		default:
			f_name = "ERROR";
  }
	if(f_nummer == 0){
		d = 9;
		x0 = 0.09;
		N = 18;
		f_nummer = 1;
	}
	else if (f_nummer > 0 && f_nummer < 5){
		cout << "Abstand d eingeben: ";
	  cin >> d;
	  cout << "Startwert x0 eingeben: ";
	  cin >> x0;
	  cout << "Anzahl N eingeben: ";
	  cin >> N;
	}
	else {
		cout << "Keine Funktion für Eingabe " << f_nummer << " bekannt. Abbruch." << endl;
		return 1;
	}

	// 2.1 b)
  y0 = f(x0);
	// cout << "(x0,y0): (" << x0 << "," << y0 << ")" << endl;

  xn = newton(x0 + d, euklid, euklid_ab);
  yn = f(xn);


	cout << endl << "f(x) = " << f_name << "; N = " << N << "; x0 = " << x0 << "; d = " << d << endl;
	cout << "Endpunkt = (xn,yn) = (" << xn << "," << yn << "); Startwert = " << x0 + d << "; Its = " << its << endl;

  // cout << "d:" << sqrt((x0-xn)*(x0-xn) + (y0-yn)*(y0-yn)) << endl; // Testen ob d richtig ist

  xk = x0;
  int L_out;
  double rc = Romberg(bogenlaenge_integrant, xk, xn, tol, 20, L_out, fcnt, bogenlaenge);
  bogenstuecklaenge = bogenlaenge/N;

	cout << "B(f;x0,xn) = " << bogenlaenge << "; Fehlertoleranz = " << tol << endl;

	// 2.1 c)
	double X[N+1], Y[N+1], X_Start_1[N+1], X_Start_2[N+1], Fcnt[N+1], Its[N+1];
	X[0] = x0;
	Y[0] = y0;
	X_Start_1[0] = x0;
	X_Start_2[0] = x0;
	Fcnt[0], Its[0] = 0;

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

		Fcnt[i+1] = fcnt;
		Its[i+1] = its;
		X_Start_1[i+1] = x_quer;
		X_Start_2[i+1] = x_dach;
	}

	// 2.1 e)
	double px, py;
	px = (x0 + X[N]) / 2;
	py = (y0 + Y[N]) / 2;

	cout << "Mittelpunkt = (px,py) = (" << px << "," << py << ")" << endl;
	cout << "xn = " << xn << endl;
	cout << "x(N) = " << X[N] << endl;
	cout << "xn-x(N) = " << xn - X[N] << endl;

	// 2.1 f)
	double Winkel[N+1], Radius[N+1], ANG[N+1];

	double a_x = x0 - px;
	double a_y = y0 - py;

	fstream file1;
	file1.open("tabelle.txt", ios::out);

	for (int i = 0; i <= N; i++) {
		Radius[i] = sqrt((px - X[i])*(px - X[i]) + (py-Y[i])*(py-Y[i]));
		double b_x = X[i] - px;
		double b_y = X[i] - py;
		Winkel[i] = acos(((a_x * b_x) + (a_y * b_y)) / (Radius[i]*d/2)) * (180.0 / PI);

		file1 << setw(3) << i << fixed << setw(WIDTH) << setprecision(PREC) << X[i] << setw(WIDTH) << Y[i] << setw(WIDTH) << X_Start_1[i] << setw(WIDTH) << X_Start_2[i] << setw(4) << setprecision(0) << Fcnt[i] << setw(2) << Its[i] << setw(WIDTH) << setprecision(PREC) << Winkel[i] << setw(WIDTH) << Radius[i] << endl;
	}

	file1.close();

	return 0;
}
