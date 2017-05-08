/**************************\
|* Autor: Pawel Karpinski *|
|* Indeks: 155085         *|
\**************************/
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cmath>

using namespace std;

const int k = 20;
const int n = k - 1;
double f[k][2];
double S[n][4];

double h(int j) {
	return f[j][0] - f[j-1][0];
}

double mi(int j) {
	if (j == n)
		return 0.0;
	return h(j) / (h(j) + h(j + 1));
}

double lambda(int j) {
	if (j == 0)
		return 0.0;
	return h(j + 1) / h(j) + h(j + 1);
}

double delta(int j) {
	if (j == 0 || j == n)
		return 0.0;
	return 6.0 / (h(j) + h(j + 1))*((f[j + 1][1] - f[j][1]) / h(j + 1) - (f[j][1] - f[j - 1][1]) / h(j));
}

double gauss_seidel_epsilon(double **M, int n, double e, double *wynik) {\
	double max;
	do {
		max = 0;
		for (int i = 0; i < n; i++) {
			double x = 0;
			for (int j = 0; j < n; j++) {
				if (j != i)
					x -= M[i][j] * wynik[j];
			}
			max = std::max(max, fabs(wynik[i] - ((x + M[i][n]) / M[i][i])));
			wynik[i] = (x + M[i][n]) / M[i][i];
		}
	} while (max>e);
	double out = wynik[0];
	return out;
}

double funkcja(double x) {
	double a, b, c, d, xj;
	for (int i = 0; i <= n; i++) {
		if (x < f[i + 1][0] || i == n) {
			a = S[i][0];
			b = S[i][1];
			c = S[i][2];
			d = S[i][3];
			xj = f[i][0];
			break;
		}
	}
	return a + b*(x - xj) + c*(x - xj)*(x - xj) + d*(x - xj)*(x - xj)*(x - xj);
}

int main() {
	double ** mat;
	mat = (double **)malloc((n + 1)*sizeof(double*));
	for (int i = 0; i <= n; i++)
		mat[i] = (double*)calloc(n + 2, sizeof(double));
	double M[n + 1];

	for (int i = 0; i < k; i++) {
		scanf_s("%lf%lf", &(f[i][0]), &(f[i][1]));
		//printf_s("%.20lf\n", 20.0*log10(f[i][1]));
	}

	for (int i = 0; i <= n; i++) {
		mat[i][i] = 2.0;
		if (i > 0)
			mat[i][i - 1] = mi(i);
		if (i < n)
			mat[i][i + 1] = lambda(i);
		mat[i][n + 1] = delta(i);
	}

	gauss_seidel_epsilon(mat, n + 1, 0.0, M);

	for (int j = 0; j <= n - 1; j++) {
		S[j][0] = f[j][1];
		S[j][1] = (f[j + 1][1] - f[j][1]) / h(j + 1) - (2.0*M[j] + M[j + 1]) / 6.0*h(j + 1);
		S[j][2] = M[j] / 2.0;
		S[j][3] = (M[j + 1] - M[j]) / 6.0 / h(j + 1);
	}

	double beg = 2.16;
	double end = 2.66;
	double iterations = 1000;

	for (double i = beg; i <= end; i += (end - beg) / iterations) {
		printf("%.20lf\t%.20lf\n", i, 20.0*log10(funkcja(i)));
	}
	
	return 0;
}