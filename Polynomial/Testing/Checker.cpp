#include<bits/stdc++.h>
#include "testlib.h"
#define all(v) (v).begin(),(v).end()
#define pb push_back
#define ppb pop_back
#define mp make_pair
#define ri(x) scanf("%d",&(x))
#define ri2(x,y) scanf("%d %d",&(x),&(y))
#define ri3(x,y,z) scanf("%d %d %d",&(x),&(y),&(z))
#define rll(x) scanf("%lld",&(x))
#define rll2(x,y) scanf("%lld %lld",&(x),&(y))
#define rll3(x,y,z) scanf("%lld %lld %lld",&(x),&(y),&(z))
#define gc(x) ((x) = getchar())
using namespace::std;

const long double PI = acos(-1);
const long long MOD = 1000000000 +7;

typedef long long ll;
typedef pair<ll,ll> pll;
typedef pair<ll,pll> tll;
typedef pair<int,int> ii;
typedef pair<int,ii> iii;
typedef vector<int> vi;
typedef vector<ii> vii;
typedef vector<iii> viii;
typedef vector<ll> vll;
typedef vector<pll> vpll;
typedef vector<tll> vtll;
typedef vector<string> vs;
typedef set<int> si;
typedef set<ii> sii;
typedef set<iii> siii;

ll gcd(ll a, ll b){ return b==0?a:gcd(b,a%b);}

ll add(ll a, ll b, ll m = MOD){
	if(a >= m) a %= m;
	if(b >= m) b %= m;
	if(a < 0) a += m;
	if(b < 0) b += m;
	ll res = a+b;
	if(res >= m or res <= -m) res %= m;
	if(res < 0) res += m;
	return res;
}

ll mul(ll a, ll b, ll m = MOD){
	if(a >= m) a %= m;
	if(b >= m) b %= m;
	if(a < 0) a += m;
	if(b < 0) b += m;
	ll res = a*b;
	if(res >= m or res <= -m) res %= m;
	if(res < 0) res += m;
	return res;
}

ll pow_mod(ll a, ll b, ll m = MOD){
	ll res = 1LL;
	a = a%m;
	while(b){
		if(b&1) res = mul(res,a,m);
		b >>= 1;
		a = mul(a,a,m);
	}
	return res;
}

ll fastexp(ll a, ll b){
	ll res = 1LL;
	while(b){
		if(b&1) res = res*a;
		b >>= 1;
		a *= a;
	}
	return res;
}

int gcdExtendido(int a, int b, int *x, int *y){
	if(a == 0){
		*x = 0;
		*y = 1;
		return b;
	}
	int x1, y1;
	int gcd = gcdExtendido(b%a,a,&x1,&y1);
	
	*x = y1-(b/a)*x1;
	*y = x1;
	return gcd;
}

int modInverso(int a, int m){
	int x, y;
	int g = gcdExtendido(a,m,&x,&y);
	if(g!=1) return -1;
	else return (x%m + m)%m;
}

/****************************************
*************P*L*A*N*T*I*L*L*A************
*****************************************/

const int N = 100+5;

int n;
long double x[N];
long double y[N];

void readData(InStream &in){
	n = in.readInt();
	for(int i=1; i<=n; i++){
		x[i] = in.readDouble();
		y[i] = in.readDouble();
	}
}

vector<long double> parse(InStream &in){
	string s = in.readLine();
	vector<long double> ans;
	long double x;
	while(s != ""){
		istringstream is(s);
		is >> x;
		ans.emplace_back(x);
		s = in.readLine();
	}
	reverse(all(ans));
	return ans;
}

const long double EPS = 5e-2;

bool test(vector<long double> &poly){
	for(int i=1; i<=n; i++){
		long double res = 0.0;
		for(int j = poly.size() - 1; j >= 0; j--){
			res *= x[i];
			res += poly[j];
		}
		long double absError = fabs(res - y[i]);
		long double relError = absError / abs(y[i]);
		long double minError = min(absError, relError);
		if(minError > EPS){
			cout << setprecision(10) << fixed << "Incorrect approximation for " << x[i] << ": P(" << x[i] << ") = " << res << " VS y(" << x[i] << ") = " << y[i] << ". Error " << minError << endl; 
			return false;
		}
	}
	return true;
}

int main(int argc, char* argv[]){
	registerTestlibCmd(argc, argv);
	readData(inf);
	vector<long double> jans = parse(ans);
	for(int i=0; i<2; i++){
		vector<long double> pans = parse(ouf);
		if(not test(pans)) quitf(_wa, "Wrong Answer. Incorrect approximation in %d", i);
	}
	quitf(_ok, "Accepted");
	return 0;
}
