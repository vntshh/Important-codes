#include<bits/stdc++.h>
#define rep(i,start,lim) for(lld i=start;i<lim;i++)
#define repd(i,start,lim) for(lld i=start;i>=lim;i--)
#define scan(x) scanf("%lld",&x)
#define print(x) printf("%lld ",x)
#define f first
#define s second
#define pb push_back
#define mp make_pair
#define br printf("\n")
#define sz(a) lld((a).size())
#define YES printf("YES\n")
#define NO printf("NO\n")
#define all(c) (c).begin(),(c).end()
using namespace std;
#define INF         1011111111
#define LLINF       1000111000111000111LL
#define EPS         (double)1e-10
#define MOD         1000000007
#define PI          3.14159265358979323
using namespace std;
typedef long double ldb;
typedef long long lld;
lld powm(lld base,lld exp,lld mod=MOD) {lld ans=1;while(exp){if(exp&1) ans=(ans*base)%mod;exp>>=1,base=(base*base)%mod;}return ans;}
typedef vector<lld> vlld;
typedef pair<lld,lld> plld;
typedef map<lld,lld> mlld;
typedef set<lld> slld;
#define N 100005
#define LOGN 20
lld log_base_2[N],a[N];
struct SparseTable{
	lld n,m1[N][LOGN];
	void pre()
	{
		rep(i,0,n) m1[i][0]=a[i];
		for(lld j=1;(1<<j)<=n;j++) for(lld i=0;i+(1<<j)-1<n;i++) 
			m1[i][j]=m1[i][j-1]|m1[i+(1<<(j-1))][j-1];
	}
	lld get(lld l,lld r)
	{
		lld tmp=log_base_2[r-l+1];
		return m1[l][tmp]|m1[r-(1<<tmp)+1][tmp];
	}
} st;
int main()
{
	lld mid,low=0,high,ans,n,tot=0;
	scan(n),st.n=n,high=n,ans=n;
	rep(i,2,N) log_base_2[i]=log_base_2[i>>1]+1;
	rep(i,0,n) cin>>a[i],tot|=a[i]; 
	st.pre();
	rep(i,0,n) {
		low=i,high=n-1;
		while(high-low>1) {
			mid=(low+high)/2;
			if(st.get(i,mid) >= tot) high=mid;
			else low=mid;
		}
		if(st.get(i,low)>=tot) ans=min(ans,(low-i+1));
		else if(st.get(i,high)>=tot) ans=min(ans,(high-i+1));
	}
	cout<<ans;
	return 0;
}
