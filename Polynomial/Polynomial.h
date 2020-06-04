#include<bits/stdc++.h>
using namespace::std;

template<class T>
class Polynomial{
	int n;
	vector<T> p;
	bool pretty;
	const long double PI = acos(-1);
	public:

	Polynomial<T>(){
		n = 1;
		pretty = true;
		p.emplace_back(T(0));
	};

	void setPretty(bool new_pretty){
		pretty = new_pretty;
	}

	static bool isZero(T val){
		if(is_floating_point<T>::value) return fabs(val) < 1e-7;
		return val == T(0);
	}

	Polynomial<T>(vector<T> v){
		while(isZero(v.back()) and v.size() > 1) v.pop_back();
		swap(p, v);
		pretty = true;
		n = p.size();
	}

	~Polynomial<T>(){
		p.clear();
		n = 0;
	}

	Polynomial<T> operator + (Polynomial<T> b) const {
		vector<T> ans(max(n, b.n), 0);
		for(int i = 0; i < ans.size(); i++){
			if(i < n) ans[i] += p[i];
			if(i < b.n) ans[i] += b.p[i];
		}
		return Polynomial<T>(ans);
	}

	void operator += (Polynomial<T> b) {
		while(n < b.n){
			p.emplace_back();
			n += 1;
		}
		for(int i = 0; i < n; i++){
			if(i < b.n){
				p[i] += b.p[i];
			}
		}
		while(isZero(p.back()) and p.size() > 1) p.pop_back();
		n = p.size();
	}

	Polynomial<T> operator - (Polynomial<T> b) const {
		vector<T> ans(max(n, b.n), 0);
		for(int i = 0; i < ans.size(); i++){
			if(i < n) ans[i] += p[i];
			if(i < b.n) ans[i] -= b.p[i];
		}
		return Polynomial<T>(ans);
	}

	void operator -= (Polynomial<T> b) {
		while(n < b.n){
			p.emplace_back();
			n += 1;
		}
		for(int i = 0; i < n; i++){
			if(i < b.n){
				p[i] -= b.p[i];
			}
		}
		while(isZero(p.back()) and p.size() > 1) p.pop_back();
		n = p.size();
	}

	Polynomial<T> operator * (Polynomial<T> b) {
		double totalSlow = 1.0 * n * b.n;
		int len = 1;
		while(len < n + b.n) len <<= 1;
		double totalFFT = 1.0 * len * log(len) / log(2);
		vector<T> res;
		if(totalSlow < totalFFT + 1e-9){
			res = multiplySlow(p, b.p);
		}
		else{
			res = multiplyFFT(p, b.p);
		}
		return Polynomial<T>(res);
	}

	Polynomial<T> operator * (T val) {
		vector<T> res(p.begin(), p.end());
		for(int i=0; i<res.size(); i++){
			res[i] *= val;
		}
		while(isZero(res.back()) and res.size() > 1) res.pop_back();
		return Polynomial<T>(res);
	}

	void operator *= (Polynomial<T> b) {
		double totalSlow = 1.0 * n * b.n;
		int len = 1;
		while(len < n + b.n) len <<= 1;
		double totalFFT = 1.0 * len * log(len) / log(2);
		vector<T> res;
		if(totalSlow < totalFFT + 1e-9){
			res = multiplySlow(p, b.p);
		}
		else{
			res = multiplyFFT(p, b.p);
		}
		p = res;
		n = p.size();
	}

	void operator *= (T val){
		for(int i=0; i<n; i++){
			p[i] *= val;
		}
		while(isZero(p.back()) and p.size() > 1) p.pop_back();
		n = p.size();
	}

	vector<T> multiplySlow(vector<T> &a, vector<T> &b){
		vector<T> res(a.size() + b.size(), T(0));
		for(int i=0; i<a.size(); i++){
			for(int j=0; j<b.size(); j++){
				res[i+j] += a[i] * b[j];
			}
		}
		while(isZero(res.back()) and res.size() > 1) res.pop_back();
		return res;
	}

	// FFT Functions

	vector<T> multiplyFFT(vector<T> &a, vector<T> &b){
		vector< complex<double> > A(a.begin(), a.end()), B(b.begin(), b.end());
		int len = 1;
		while(len < a.size() + b.size()){
			len <<= 1;
		}
		A.resize(len);
		B.resize(len);
		fft(A, false);
		fft(B, false);
		for(int i=0; i<len; i++){
			A[i] *= B[i];
		}
		fft(A, true);
		vector<T> res(len);
		for(int i = 0; i < len; i++){
			res[i] = (T)(A[i].real() + 1e-9);
		}
		while(isZero(res.back()) and res.size() > 1) res.pop_back();
		return res;
	}

	int reverse(int num, int lg_n){
		int res = 0;
		for(int i = 0; i < lg_n; i++){
			if(num & (1<<i)){
				res |= 1<<(lg_n - i - 1);
			}
		}
		return res;
	}

	void fft(vector< complex<double> > &a, bool invert){
		int len = a.size();
		int lg_n = 0;
		while((1<<lg_n) < len) lg_n += 1;

		assert((1<<lg_n) == len); // Len must be a power of 2

		for(int i = 0; i < len; i++){
			if(i < reverse(i, lg_n)){
				swap(a[i], a[reverse(i, lg_n)]);
			}
		}

		for(int l = 2; l <= len; l <<= 1){
			double angle = 2 * PI / l * (invert? -1 : 1);
			complex<double> wlen(cos(angle), sin(angle));
			for(int i = 0; i < len; i += l){
				complex<double> w(1);
				for(int j = 0; j + j < l; j++){
					complex<double> u = a[i + j], v = w * a[i + j + l / 2];
					a[i + j] = u + v;
					a[i + j + l / 2] = u - v;
					w *= wlen;
				}
			}
		}
		if(invert){
			for(complex<double> &x : a){
				x /= len;
			}
		}
	}

	static bool zero(const Polynomial<T> &v){
		return v.n == 1 and isZero(v.p.back());
	}

	// I/O Functions

	friend ostream &operator<<(ostream &output, const Polynomial<T> &v) {
		if(is_floating_point<T>::value) output << setprecision(10) << fixed;
		if(Polynomial<T>::zero(v)) output << T(0);
		else{
			if(v.pretty){
				output << "Polynomial of degree " << v.n-1 << ": \n";
				for(int i=v.n-1; i>=0; i--){
					if(Polynomial<T>::isZero(v.p[i])) continue;
					if(i < v.n - 1){
						if(v.p[i] < -1e-9){
							output << " - ";
							if(fabs(v.p[i] + 1.0) > 1e-9 or i == 0) output << fabs(v.p[i]);
						}
						else{
							output << " + ";
							if(fabs(v.p[i] - 1.0) > 1e-9 or i == 0) output << v.p[i];
						}
					}
					else{
						if(v.p[i] < -1e-9){
							output << "-";
							if(fabs(v.p[i] + 1.0) > 1e-9 or i == 0) output << fabs(v.p[i]);
						}
						else{
							if(fabs(v.p[i] - 1.0) > 1e-9 or i == 0) output << v.p[i];
						}
					}
					if(i > 0){
						output << "x";
						if(i > 1) output << "^" << i;
					}
				}
			}
			else{
				for(int i=v.n-1; i>=0; i--){
					output << v.p[i] << endl;
				}
			}
		}
		return output;
	}

	friend istream &operator>>(istream &input, Polynomial<T> &v){
		input >> v.n;
		v.n++;
		v.p.resize(v.n);
		for(int i=v.n-1; i>=0; i--){
			input >> v.p[i];
		}
		return input;
	}

	// Access Operators
	
	T& lead(){
		assert(not p.empty());
		return p.back();
	}

	T operator [] (int idx) const {
		return idx >= n or idx < 0? T(0) : p[idx];
	}

	bool operator == (const Polynomial<T> &t) const { return p == t.p;}

	bool operator != (const Polynomial<T> &t) const { return p != t.p;}

	// Calculus Functions
	
	Polynomial<T> deriv(){
		int new_n = max(1, n-1);
		vector<T> ans(new_n, T(0));
		for(int i = 1; i < n; i++){
			ans[i-1] = i * p[i];
		}
		return Polynomial<T>(ans);
	}

	Polynomial<long double> integr(long double val){
		int new_n = n + 1;
		vector<long double> ans(new_n, 0);
		for(int i=0; i<n; i++){
			ans[i+1] = 1.0 * p[i] / (i + 1);
		}
		ans[0] = val;
		return Polynomial<long double>(ans);
	}

	long double eval(long double x){
		long double ans = 0.0;
		for(int i = n-1; i >= 0; i--){
			ans *= x;
			ans += p[i];
		}
		return ans;
	}
};
