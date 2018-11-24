#include <bits/stdc++.h>

using namespace std;

typedef complex<long double> Complex;
typedef long long ll;
const long double PI = acos(-1.0L);
const int mod = 7340033;
const int g = 5; // g^{2^{20}} = 1 mod p
const int gr = 4404020; // g^{p-2} mod p
const int gpw = 1 << 20; // g^{gpw} = 1 mod p

void print(const vector<Complex>& a) {
  for (const Complex& a_k : a) {
    cout << a_k.real() << " ";
  }
  cout << endl;
}

void print(const vector<int>& a) {
  for (const int& a_k : a) {
    cout << a_k << " ";
  }
  cout << endl;
}

vector<Complex> dft(vector<Complex> a, int type) {
  int n = a.size();
  if (n == 1) {
    return a;
  }
  Complex wn(cos(2.0L * type * PI / n), sin(2.0L * type * PI / n));
  Complex w(1.0L, 0.0L);
  vector<Complex> a_even, a_odd;
  for (int k = 0; k <= n-2; k += 2) {
    a_even.push_back(a[k]);
    a_odd.push_back(a[k+1]);
  }
  vector<Complex> y_even = dft(a_even, type);
  vector<Complex> y_odd = dft(a_odd, type);
  vector<Complex> y(n);

  for (int k = 0; k < n/2; k++) {
    y[k] = y_even[k] + w * y_odd[k];
    y[k + n/2] = y_even[k] - w * y_odd[k];
    w = w * wn;
  }

  return y;
}

vector<Complex> dft(vector<Complex> a) {
  return dft(a, 1);
}

vector<Complex> idft(vector<Complex> a) {
  int n = a.size();
  vector<Complex> c = dft(a, -1);
  for(Complex &c_k : c) {
    c_k /= n;
  }
  return c;
}

template <typename T>
vector<T> BitReversalPermutation(const vector<T>& a) {
  int n = a.size();
  int logn = __builtin_ctz(n);
  vector<T> v(n);
  for (int k = 0; k < n; k++) {
    int rev = 0; // rev = reverse(k)
    for (int j = 0; j < logn; j++) 
      if (k&(1<<j)) rev |= (1<<(logn - 1 - j)); 
    v[rev] = a[k];
  }
  return v;
}

vector<Complex> fft(vector<Complex> a, const int r) {
  int n = a.size();
  int logn = __builtin_ctz(n); // n is a power of two
  vector<Complex> v = BitReversalPermutation<Complex>(a);
  for (int s = 1; s <= logn; s++) {
    int m = (1 << s);
    Complex wm(cos(2.0 * r * PI / m), sin(2.0 * r * PI / m));
    for (int k = 0; k < n; k += m) {
      Complex w = 1;
      for (int j = 0; 2*j < m; j++) {
        Complex t = w * v[k + j + m/2];
        Complex u = v[k + j];
        v[k + j] = u + t;
        v[k + j + m/2] = u - t;
        w *= wm;
      }
    }
  }
  if (r == -1)
    for (Complex& c : v) c /= n;
  return v;
}

ll fexp (ll n, ll p) {
    if (p == 0) return 1;
    if (p == 1) return n % mod;
    if (p % 2 == 0) return fexp (n*n % mod, p/2) % mod;
    return (n * fexp (n*n % mod, p/2)) % mod;
}

vector<int> fft(vector<int> a, const int rev) {
  bool r = rev == -1;
  int n = a.size();
  int logn = __builtin_ctz(n);
  vector<int> v = BitReversalPermutation<int>(a);
  for (int s = 1; s <= logn; s++) {
    int m = 1 << s;
    int zm = r ? gr : g;
    int cont = 0;
    for (int i = s; i < 20; i++)
      zm = (int)(1LL * zm * zm % mod);
    for (int k = 0; k < n; k += m) {
      int z = 1;
      for (int j = 0; 2 * j < m; j++) {
        int t = (int)(1LL * v[k + j + m/2] * z % mod);
        int u = v[k + j];
        v[k + j] = u + t < mod ? u + t : u + t - mod;
        v[k + j + m/2] = u - t >= 0 ? u - t : u - t + mod;
        z = (int)(1LL * z * zm % mod);
      }
    }
  }
  if (r == true) {
    int n_1 = (int) fexp(n, mod - 2);
    for (int& x : v) x = (int) (1LL * x * n_1 % mod) ;
  }
  return v;
}

Complex mult(Complex a, Complex b) {
  return a * b;
}

int mult(int a, int b) {
  return (int)(1LL * a * b % mod);
}

template<typename T>
vector<T> fast_multiplication(const vector<T>& a, const vector<T>& b) {
  vector<T> dft_a = fft(a, 1);
  vector<T> dft_b = fft(b, 1);
  vector<T> dft_c(dft_a.size());
  for (int k = 0; k < dft_a.size(); k++) {
    dft_c[k] = mult(dft_a[k], dft_b[k]);
  }
  vector<T> c = fft(dft_c, -1);
  return c;
}

int main() {
  vector<Complex> a = {9, 10, 7, 6, 0, 0, 0, 0};
  vector<Complex> b = {5, 4, 1, 1, 0, 0, 0, 0};
  vector<Complex> c = fast_multiplication<Complex>(a, b);
  print(c);
  vector<int> x = {9, 10, 7, 6, 0, 0, 0, 0};
  vector<int> y = {5, 4, 1, 1, 0, 0, 0, 0};
  vector<int> z = fast_multiplication<int>(x, y);
  print(z);
}


