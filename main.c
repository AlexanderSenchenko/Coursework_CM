#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float eps = 1E-2;

float d2f(float x, float y, float d1, float d2)
{
	return powf(d2, 3) - 9 * sin(x) * d2 - sin(x) + 9 * exp(x) * d1 - (y / powf(x + 2, 3));
}

float EquationLine(float a, float b, float x, float y, float d1)
{
	float fa = d2f(x, y, d1, a);
	float fb = d2f(x, y, d1, b);
	return ((a * fb - b * fa) / (fb - fa));
}

float ChordsMethod(float a, float b, float x, float y, float d1)
{
	float fa, fc, c = 0, z = 0;
	// float eps = 1E-5;

	while (fabsf(c - z) > eps) {
		fc = d2f(x, y, d1, c);
		fa = d2f(x, y, d1, a);
		if (fa * fc < 0) {
			b = c;
		} else {
			a = c;
		}

		z = c;
		c = EquationLine(a, b, x, y, d1);
	}
	return c;
}

float Method_half_division(float a, float b, float x, float y, float d1)
{
	double c = 0;
	while(fabs(b - a) > eps) {
		c = (a + b) / 2;
		if(d2f(x, y, d1, a) * d2f(x, y, d1, c) < 0) {
			b = c;
		} else if (d2f(x, y, d1, c) * d2f(x, y, d1, b) < 0) {
			a = c;
		}
	}
	return (a + b) / 2;
}

float f(float x, float y, float d1)
{
	float a = 1, b = 0;
	float fa = 0, fb = 0;
	do {
		fa = d2f(x, y, d1, a--);
		fb = d2f(x, y, d1, b++);
	} while (fa * fb > 0);

	// return ChordsMethod(a, b, x, y, d1);
	return Method_half_division(a, b, x, y, d1);
}

double out_d1 = 0;

float RungeKutt2_time(float x0, float x1, float h, float y, float d1)
{
	float x = x0;
	float yt = 0, d1t = 0;
	for (; x < x1; x += h) {
		yt = y + (h / 2) * d1;
		d1t = d1 + (h / 2) * f(x, y, d1);
		y += h * d1t;
		d1 += h * f(x + h / 2, yt, d1t);
	}

	h = x1 - x;
	yt = y + (h / 2) * d1;
	d1t = d1 + (h / 2) * f(x0, y, d1);
	y += h * d1t;
	// out_d1 = d1 + h * f(x1 + h / 2, yt, d1t);

	return y;
}

float ShootingMethod(float x0, float x1, float y0, float y1, float h)
{
	float a = 1, b = 0;
	float fa = 0, fb = 0;
	// float eps = 1E-5;

	do {
		fa = RungeKutt2_time(x0, x1, h, y0, a) - y1;
		fb = RungeKutt2_time(x0, x1, h, y0, b) - y1;
		a -= h;
		b += h;
	} while (fa * fb > 0);

	float g = 0;
	while(fabs(b - a) > eps) {
		g = (a + b) / 2;
		if((RungeKutt2_time(x0, x1, h, y0, a) - y1) * 
			(RungeKutt2_time(x0, x1, h, y0, g) - y1) < 0) {
			b = g;
		} else if ((RungeKutt2_time(x0, x1, h, y0, g) - y1) * 
			(RungeKutt2_time(x0, x1, h, y0, b) - y1) < 0) {
			a = g;
		}
	}

	return (a + b) / 2;
}

float DoubleCounting(float x0, float x1, float y0, float y1, float h)
{
	float delta = 1;
	// float eps = 1E-5;
	float d1[2];
	int k = 0;
	for (float h0 = h; delta >= eps; h0 /= 2, k ^= 1) {
		d1[k] = ShootingMethod(x0, x1, y0, y1, h);

		if (h > h0)
			delta = fabsf(d1[k] - d1[k ^ 1]);
	}

	return d1[k];
}

float Lagrange(float z, float *x, float *y, int n)
{
	float *Px = calloc(n, sizeof(float));
	float fz = 0.0, mult;

	for (int i = 0; i < n; i++) {
		mult = 1.0;
		for (int j = 0; j < n; j++) {
			if (i != j)
				mult *= (z - x[j]) / (x[i] - x[j]);
		}
		Px[i] = mult * y[i];
	}

	for (int i = 0; i < n; i++) {
		fz += Px[i];
	}

	free(Px);

	return fz;
}

float Integr_f(float x) { return powf(x, 2); }

float NIntegr(float a, float b)
{
	float n0 = 1000000; 
	float n = n0;
	float h;

	float S0 = 0;
	float delta = 1;
	float eps = 1E-2;

	for (;delta >= eps; n *= 2) {
		float S1 = 0;
		h = (b - a) / n;
		for (int i = 0; i < n; i++) {
			S1 += (Integr_f(a + i * h) + Integr_f(a + (i + 1) * h)) / 2;
		}
		S1 *= h;

		delta = fabsf(S1 - S0);

		S0 = S1;
	}

	return S0;
}

int main()
{
	float x0 = 0, x1 = 1, y0 = 3, y1 = 3;
	float h = 0.2;

	float d1 = DoubleCounting(x0, x1, y0, y1, h);

	int n = 6;
	float *X = malloc(sizeof(float) * n);
	float *Y = malloc(sizeof(float) * n);
	for (int i = 0; i < n; i++) {
		X[i] = h * i;
		Y[i] = RungeKutt2_time(x0, X[i], h, y0, d1);
	}

	FILE *res = fopen("res.txt", "w");
	for (int i = 0; i < n ; i++) {
		fprintf(res, "%.3f %.3f\n", X[i], Y[i]);
	}
	fclose(res);

	FILE *out = fopen("out.txt", "w");
	for (float i = x0; i <= x1 ; i += (h / 10)) {
		fprintf(out, "%.3f %.3f\n", i, Lagrange(i, X, Y, n));
	}
	fclose(out);



	system("gnuplot scen.plt");

	free(X);
	free(Y);

	return 0;
}