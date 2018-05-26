#ifndef MAIN_H
#define MAIN_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Lagrange.h"

float out_d1 = 0;
float eps = 1E-2;

float d2f(float x, float y, float d1, float d2);
float Method_half_division(float a, float b, float x, float y, float d1);
float f(float x, float y, float d1);
float RungeKutt2_time(float x0, float x1, float h, float y, float d1);
float ShootingMethod(float x0, float x1, float y0, float y1, float h);
float DoubleCounting(float x0, float x1, float y0, float y1, float h);
float* DoubleCountingRunge(float *X, int n, float h, float x0, float x1, float y0, float d1);

float NIntegr(float a, float b, float y0, float d1);
float simpsons_rule(float *XIntegr, float *YIntegr, int n, float h);

#endif