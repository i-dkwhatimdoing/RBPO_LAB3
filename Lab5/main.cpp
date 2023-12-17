#include <iostream>
#include <complex>
import Math;
using namespace std;
int main()
{
	double re, im;
	cout << "Input Real and Imag num:" << endl;
	cin >> re >> im;
	Complex c_num(re, im);

	double chisl, znam;
	cout << "Input Chisl and Znam num:" << endl;
	cin >> chisl >> znam;
	Rational r_num(chisl, znam);
	double d_num;
	cout << "Input double:" << endl;
	cin >> d_num;

	cout << "f(" << c_num << ") = " << f(c_num) << endl;
	cout << "f(" << r_num << ") = " << f(r_num) << endl;
	cout << "f(" << d_num << ") = " << f(d_num) << endl;

	return 0;
}