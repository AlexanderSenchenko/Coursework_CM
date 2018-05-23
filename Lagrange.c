#include "Lagrange.h"

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