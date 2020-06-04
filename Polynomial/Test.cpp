#include<bits/stdc++.h>
#include "Interpolation.h"
using namespace::std;

int main(){
	int len;
	vector< pair<long double, long double> > points;
	cin >> len;
	for(int i=0; i<len; i++){
		long double x, y;
		cin >> x >> y;
		points.emplace_back(make_pair(x, y));
	}
	Polynomial<long double> A = Lagrange(points);
	A.setPretty(false);
	//cout << "Lagrange's interpolation polynomial: " << endl;
	cout << A << endl;/*
	for(auto e : points){
		cout << setprecision(10) << fixed << A.eval(e.first) << " " << e.second << endl;
	}*/
	Polynomial<long double> B = Newton(points);
	B.setPretty(false);
	//cout << "Newton's interpolation polynomial: " << endl;
	cout << B << endl;/*
	for(auto e : points){
		cout << setprecision(10) << fixed << B.eval(e.first) << " " << e.second << endl;
	}*/
	return 0;
}
