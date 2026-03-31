#include <bits/stdc++.h>
using namespace std;

struct Term {
    long long a; // coefficient
    int b, c, d; // exponents for x, sin, cos
    bool operator==(Term const& o) const { return b==o.b && c==o.c && d==o.d; }
};

struct Poly {
    vector<Term> t;
    void simplify(){
        if(t.empty()) return;
        sort(t.begin(), t.end(), [](const Term& x, const Term& y){
            if(x.b!=y.b) return x.b<y.b;
            if(x.c!=y.c) return x.c<y.c;
            return x.d<y.d;
        });
        vector<Term> v; v.reserve(t.size());
        for(auto cur: t){
            if(!v.empty() && v.back()==cur){
                v.back().a += cur.a;
            } else {
                v.push_back(cur);
            }
        }
        t.clear();
        for(auto &x: v) if(x.a!=0) t.push_back(x);
    }
};

static Poly add(const Poly& A, const Poly& B){ Poly R=A; R.t.insert(R.t.end(), B.t.begin(), B.t.end()); R.simplify(); return R; }
static Poly sub(const Poly& A, const Poly& B){ Poly R=A; for(auto x: B.t){ R.t.push_back({-x.a,x.b,x.c,x.d}); } R.simplify(); return R; }
static Poly mul(const Poly& A, const Poly& B){
    Poly R; R.t.reserve((size_t)A.t.size()*B.t.size());
    for(auto &x: A.t) for(auto &y: B.t){
        long long a = x.a * y.a;
        int b = x.b + y.b;
        int c = x.c + y.c;
        int d = x.d + y.d;
        R.t.push_back({a,b,c,d});
    }
    R.simplify();
    return R;
}

static Poly derivate_poly(const Poly& P){
    Poly R;
    for(const auto &x: P.t){
        // x^b term: b*a x^{b-1} sin^c cos^d
        if(x.b>0){ R.t.push_back({x.a * x.b, x.b-1, x.c, x.d}); }
        // sin^c: c*a x^b sin^{c-1} cos^1
        if(x.c>0){ R.t.push_back({x.a * x.c, x.b, x.c-1, x.d+1}); }
        // cos^d: d*a x^b sin^1 cos^{d-1} with negative sign
        if(x.d>0){ R.t.push_back({-x.a * x.d, x.b, x.c+1, x.d-1}); }
    }
    R.simplify();
    return R;
}

struct Frac{ Poly p,q; };

static Frac make_int(long long v){ Frac f; f.p.t.push_back({v,0,0,0}); f.q.t.push_back({1,0,0,0}); return f; }

static Frac add(const Frac& A, const Frac& B){
    Frac R; R.p = add(mul(A.p,B.q), mul(B.p,A.q)); R.q = mul(A.q,B.q); return R;
}
static Frac sub(const Frac& A, const Frac& B){
    Frac R; R.p = sub(mul(A.p,B.q), mul(B.p,A.q)); R.q = mul(A.q,B.q); return R;
}
static Frac mul(const Frac& A, const Frac& B){ Frac R; R.p = mul(A.p,B.p); R.q = mul(A.q,B.q); return R; }
static Frac divi(const Frac& A, const Frac& B){ Frac R; R.p = mul(A.p,B.q); R.q = mul(A.q,B.p); return R; }
static Frac derivate_frac(const Frac& F){ Frac R; Poly p1=derivate_poly(F.p); Poly q1=derivate_poly(F.q); R.p = sub(mul(p1,F.q), mul(q1,F.p)); R.q = mul(F.q,F.q); return R; }

// Output formatting
static void output_poly(const Poly& P){
    if(P.t.empty()){ cout << 0; return; }
    // Sort for output: x-terms by b desc, then sin-terms by c desc, then cos-terms by d desc, then constants
    auto group_of = [](const Term& t){ if(t.b>0) return 0; if(t.c>0) return 1; if(t.d>0) return 2; return 3; };
    vector<Term> out = P.t;
    sort(out.begin(), out.end(), [&](const Term& u, const Term& v){
        int gu = group_of(u), gv = group_of(v);
        if(gu!=gv) return gu<gv;
        if(gu==0){ if(u.b!=v.b) return u.b>v.b; if(u.c!=v.c) return u.c>v.c; if(u.d!=v.d) return u.d>v.d; }
        else if(gu==1){ if(u.c!=v.c) return u.c>v.c; if(u.d!=v.d) return u.d>v.d; if(u.b!=v.b) return u.b>v.b; }
        else if(gu==2){ if(u.d!=v.d) return u.d>v.d; if(u.c!=v.c) return u.c>v.c; if(u.b!=v.b) return u.b>v.b; }
        // constants or tie-breaker by coefficient magnitude
        return llabs(u.a)>llabs(v.a);
    });
    bool first=true;
    for(const auto &x: out){
        long long a=x.a; int b=x.b,c=x.c,d=x.d;
        if(a==0) continue;
        if(first){
            if(a<0){ cout << '-'; a = -a; }
            first=false;
        }else{
            cout << (a<0?'-':'+'); if(a<0) a=-a; }
        bool printed=false;
        if(!(b>0 || c>0 || d>0)){
            cout << a; printed=true;
        } else {
            if(llabs(a)!=1){ cout << llabs(a); printed=true; }
            else if(a==-1){ /* handled sign earlier */ }
        }
        // x^b sin^c x cos^d x in required order: x, sin^k x, cos^k x
        if(b>0){ cout << 'x'; if(b>1) cout << '^' << b; printed=true; }
        if(c>0){ cout << "sin"; if(c>1) cout << '^' << c; cout << 'x'; printed=true; }
        if(d>0){ cout << "cos"; if(d>1) cout << '^' << d; cout << 'x'; printed=true; }
        if(!printed){ cout << 1; }
    }
}

static void output_frac(const Frac& F){
    bool num_single = F.p.t.size()==1;
    bool den_single = F.q.t.size()==1;
    auto is_one = [](const Poly& P){ return P.t.size()==1 && P.t[0].a==1 && P.t[0].b==0 && P.t[0].c==0 && P.t[0].d==0; };
    if(is_one(F.q)){
        output_poly(F.p); return;
    }
    cout << '('; output_poly(F.p); cout << ')' << '/';
    cout << '('; output_poly(F.q); cout << ')';
}

// Parsing utilities
static bool isnum(char c){ return c>='0'&&c<='9'; }

static long long get_number(const string& s, int l, int r){
    // parse leading integer coefficient in s[l:r). If absent, return +/-1 based on leading sign.
    int i=l; bool neg=false; if(i<r && (s[i]=='+'||s[i]=='-')){ neg = (s[i]=='-'); ++i; }
    long long val=0; bool has=false;
    while(i<r && isnum(s[i])){ has=true; val = val*10 + (s[i]-'0'); ++i; }
    if(!has) return neg?-1:1;
    return neg? -val: val;
}

static Term get_term(const string& s, int l, int r){
    // parse a term ax^b sin^c x cos^d x; input segment guaranteed no top-level +/- or * or /
    long long a = get_number(s,l,r);
    int b=0,c=0,d=0;
    // move cursor past coefficient (optional sign+digits)
    int i=l; if(i<r && (s[i]=='+'||s[i]=='-')) ++i; while(i<r && isnum(s[i])) ++i;
    while(i<r){
        if(i+2<r && s[i]=='s' && s[i+1]=='i' && s[i+2]=='n'){
            i+=3; int k=1; if(i<r && s[i]=='^'){ ++i; int j=i; while(j<r && isnum(s[j])) ++j; k = stoi(s.substr(i,j-i)); i=j; }
            if(i<r && s[i]=='x') ++i; c += k; continue;
        }
        if(i+2<r && s[i]=='c' && s[i+1]=='o' && s[i+2]=='s'){
            i+=3; int k=1; if(i<r && s[i]=='^'){ ++i; int j=i; while(j<r && isnum(s[j])) ++j; k = stoi(s.substr(i,j-i)); i=j; }
            if(i<r && s[i]=='x') ++i; d += k; continue;
        }
        if(s[i]=='x'){
            ++i; int k=1; if(i<r && s[i]=='^'){ ++i; int j=i; while(j<r && isnum(s[j])) ++j; k = stoi(s.substr(i,j-i)); i=j; }
            b += k; continue;
        }
        ++i;
    }
    return Term{a,b,c,d};
}

static Frac dfs(const string& s, int l, int r){
    // remove outer parentheses
    while(l<r && s[l]=='(' && s[r-1]==')'){
        int dep=0; bool ok=true; for(int i=l;i<r;i++){ if(s[i]=='(') ++dep; if(s[i]==')'){ --dep; if(dep==0 && i!=r-1){ ok=false; break; } } }
        if(ok){ ++l; --r; } else break;
    }
    // find + or - at top level (excluding the first char if it's unary)
    int dep=0;
    for(int i=r-1;i>=l;--i){
        char ch=s[i];
        if(ch==')') ++dep; else if(ch=='(') --dep;
        else if(dep==0 && (ch=='+'||ch=='-')){
            if(i==l) continue; // leading sign belongs to term
            Frac L = dfs(s,l,i);
            Frac Rf = dfs(s,i+1,r);
            if(ch=='+') return add(L,Rf);
            else return sub(L,Rf);
        }
    }
    // find * or /
    dep=0; for(int i=r-1;i>=l;--i){ char ch=s[i]; if(ch==')') ++dep; else if(ch=='(') --dep; else if(dep==0 && (ch=='*'||ch=='/')){ Frac L=dfs(s,l,i); Frac Rf=dfs(s,i+1,r); if(ch=='*') return mul(L,Rf); else return divi(L,Rf);} }
    // fall to a term segment
    Term t = get_term(s,l,r);
    Frac f; f.p.t.push_back(t); f.q.t.push_back({1,0,0,0}); f.p.simplify(); f.q.simplify(); return f;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string expr; if(!getline(cin, expr)) return 0; // read line
    if(expr.size()==0){ return 0; }
    Frac f = dfs(expr,0,(int)expr.size());
    Frac g = derivate_frac(f);
    output_frac(f); cout << "\n"; output_frac(g); cout << "\n";
    return 0;
}
