#include<bits/stdc++.h>
#include "Polynomial.h"
using namespace::std;

Polynomial<long double> Lagrange(vector< pair<long double, long double> > &points){
	// Lagrange Interpotation O(n^3) where n is the number of points
	int n = points.size();
	Polynomial<long double> ans;
	for(int i = 0; i < n; i++){
		Polynomial<long double> cur({points[i].second});
		for(int j = 0; j < n; j++){
			if(i == j) continue;
			if(fabs(points[i].first - points[j].first) < 1e-10){
				cout << "The given points repeat X coordinate, it's not a function" << endl;
				assert(false);
			}
			Polynomial<long double> linear({-points[j].first, 1});
			linear *= 1.0 / (points[i].first - points[j].first);
			cur *= linear;
		}
		ans += cur;
	}
	return ans;
}

Polynomial<long double> Newton(vector< pair<long double, long double> > &points){
	// Newton Interpotation O(n^2) where n is the number of points
	int n = points.size();
	vector< vector<long double> > memo(n, vector<long double>(n));
	for(int i = 0; i < n; i++) memo[i][0] = points[i].second;
	for(int len = 1; len < n; len++){
		for(int i=0; i + len < n; i++){
			memo[i][len] = (memo[i + 1][len-1] - memo[i][len - 1]) / (points[i + len].first - points[i].first);
		}
	}
	Polynomial<long double> ans;
	Polynomial<long double> prefix({1});
	for(int i = 0; i < n; i++){
		ans += prefix * memo[0][i];
		if(i + 1 < n) prefix *= Polynomial<long double>({-points[i].first, 1});
	}
	return ans;
}
