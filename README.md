# Numerical Programming Documentation 
**Contents**

* NP myNP.h
* NP matrix.h

<hr>

## NP myNP.h

```c
#include "myNP.h"
```

**List**

* factorial()
* power()
* sinTaylor()
* sindTaylor()
* bisection()
* newtonRaphson()
* secantfzero()
* printVec()
* gradient1D()
* gradientFunc()
* acceleration()
* gradientFunck()
* IntegrateRect()
* trapz()
* simpson13()
* integral()
* odeRK2()
* odeRK3()
* sysRK2()

<hr>

### factorial()

> Multiplies a given natural number _x with every natural number less than or equal

```c
double factorial(int _x);
```
**Parameter**

* _x: given natural number

**Example**

```C
#include "myNP.h"

int main()
{
    int N = 5; // N should be a natural number
    double y;
    y = factorial(N);
    
    printf("\n\n");
	printf("=======================================\n");
	printf("    factorial(%d) Calculation   \n", N);
	printf("=======================================\n");
	printf("   -  My result = %.0f    \n", y);
	printf("=======================================\n");
}
```

**Output**

```c
=======================================
    factorial(5) Calculation
=======================================
   -  My result = 120
=======================================
```

<hr>

### power()

> Multiplies given number _x with N times

```c
double power(double _x, int N);
```

**Parameter**

* _x: base
* N: exponent

**Example**

```C
#include "myNP.h"

int main()
{
    double x = 2.5;
    int N = 5; // N should be a natural number
    double y;
    y = power(x, N);
    
    printf("\n\n");
	printf("=======================================\n");
	printf("    power(%f,%d) Calculation   \n", x, N);
	printf("=======================================\n");
	printf("   -  My     result = %3.12f    \n", y);
	printf("   -  Math.h result = %3.12f    \n", pow(x, N));
	printf("   -  absolute err. = %3.12f    \n", y - pow(x, N));
	printf("=======================================\n");
}
```

**Output**

```c
=======================================
    power(2.500000,5) Calculation
=======================================
   -  My     result = 97.656250000000
   -  Math.h result = 97.656250000000
   -  absolute err. = 0.000000000000
=======================================
```

<hr>

### sinTaylor()

> Taylor series approximation for sin(x) (input unit: [rad])

```c
double sinTaylor(double _x);
```

**Parameter**

* _x: the given angle, the unit should be [rad]

**Example**

```C
#include "myNP.h"

int main()
{
    double x = PI / 6;
    double S_N;
    S_N = sinTaylor(x);
    
    printf("\n\n");
	printf("=======================================\n");
	printf("    sin( %f[rad] ) Calculation   \n", x);
	printf("=======================================\n");
	printf("   -  My     result = %3.12f    \n", S_N);
	printf("   -  Math.h result = %3.12f    \n", sin(x));
	printf("   -  absolute err. = %3.12f    \n", S_N - sin(x));
	printf("=======================================\n");
}
```

**Output**

```c
=======================================
    sin( 0.523599[rad] ) Calculation
=======================================
   -  My     result = 0.500000000000
   -  Math.h result = 0.500000000000
   -  absolute err. = 0.000000000000
=======================================
```

<hr>

### sindTaylor()

> Taylor series approximation for sin(x) (input unit: [deg])

```c
double sindTaylor(double _x)
```

**Parameter**

* _x: the given angle, the unit should be [deg]

**Example**

```C
#include "myNP.h"

int main()
{
    double x = 30;
    double S_N;
    S_N = sindTaylor(x);
    
    printf("\n\n");
	printf("=======================================\n");
	printf("    sin( %f[deg] ) Calculation   \n", x);
	printf("=======================================\n");
	printf("   -  My     result = %3.12f    \n", S_N);
	printf("   -  Math.h result = %3.12f    \n", sin(x * PI /180));
	printf("   -  absolute err. = %3.12f    \n", S_N - sin(x * PI / 180));
	printf("=======================================\n");
}
```

**Output**

```c
=======================================
    sin( 30.000000[deg] ) Calculation
=======================================
   -  My     result = 0.500000000000
   -  Math.h result = 0.500000000000
   -  absolute err. = 0.000000000000
=======================================
```

<hr>

### bisection()

> Non-linear solver using Bisection method

```c
double bisection(double func(const double x), double _a0, double _b0, double _tol);
```

**Parameter**

* func: mathematic function
* _a0: start of interval
* _b0: end of interval
* _tol: tolerance condition

**Example**

```C
#include "myNP.h"

int main()
{
    float tol = 0.00001;
	float a0 = 2;
    float b0 = 3;
    double sol_bm;
    
	printf("Bisection Method:\n");
    sol_bm = bisection(func, a0, b0, tol); // func: F = 8 - 4.5 * (x - sin(x))
	printf("Final Solution: %f \t", sol_bm);
	printf("\n");
}
```

**Output**

```c
Bisection Method:
k:1     Xn(k): 2.500000         Tol: 0.55687535
k:2     Xn(k): 2.250000         Tol: 1.37632942
k:3     Xn(k): 2.375000         Tol: 0.43408266
k:4     Xn(k): 2.437500         Tol: 0.05570867
k:5     Xn(k): 2.406250         Tol: 0.19066088
k:6     Xn(k): 2.421875         Tol: 0.06783819
k:7     Xn(k): 2.429688         Tol: 0.00615447
k:8     Xn(k): 2.433594         Tol: 0.02475477
k:9     Xn(k): 2.431641         Tol: 0.00929455
k:10    Xn(k): 2.430664         Tol: 0.00156864
k:11    Xn(k): 2.430176         Tol: 0.00229327
k:12    Xn(k): 2.430420         Tol: 0.00036240
k:13    Xn(k): 2.430542         Tol: 0.00060310
k:14    Xn(k): 2.430481         Tol: 0.00012034
k:15    Xn(k): 2.430450         Tol: 0.00012103
k:16    Xn(k): 2.430466         Tol: 0.00000034
Final Solution: 2.430466
```

<hr>

### newtonRaphson()

> Non-linear solver using Newton-Raphson method

```c
double newtonRaphson(double func(const double x), double dfunc(const double x), double _x0, double _tol);
```

**Parameter**

* func: mathematic function
* dfunc: derivative of func
* _x0: initial value
* _tol: tolerance condition

**Example**

```C
#include "myNP.h"

int main()
{
    float tol = 0.00001;
	double x0 = 3;
    double sol_nr;
    
	printf("Newton-Raphson Method Result:\n");
    sol_nr = newtonRaphson(func, dfunc, x0, tol);
    // func: F = 8 - 4.5 * (x - sin(x)), dfunc: dF = -4.5 * (1 - cos(x))
	printf("Final Solution: %f \t", sol_nr);
	printf("\n");
}
```

**Output**

```c
Newton-Raphson Method Result:
k:1     x(k): 2.456731  tol: 0.2087399364
k:2     x(k): 2.430590  tol: 0.0009821436
k:3     x(k): 2.430466  tol: 0.0000003439
Final Solution: 2.430466
```

<hr>

### secantfzero()

> Non-linear solver using Secant method

```c
double secantfzero(double func(const double x), double _x0, double _x1, double _tol);
```

**Parameter**

* func: mathematic function
* _x0: 1st initial value
* _x1: 2nd initial value
* _tol: tolerance condition

**Example**

```C
#include "myNP.h"

int main()
{
	float tol = 0.00001;
    double x0 = 0;
	double x1 = 3;
    double sol_sz;
    
    printf("Secant Method Result:\n");
    sol_sz = secantfzero(func, x0, x1, tol); // func: F = 8 - 4.5 * (x - sin(x))
	printf("Final Solution: %f \t", sol_sz);
	printf("\n");
}
```

**Output**

```c
Secant Method Result:
k:1     Xn(k): 1.865532         Tol: 3.91105890
k:2     Xn(k): 2.371111         Tol: 0.46416354
k:3     Xn(k): 2.439193         Tol: 0.06913824
k:4     Xn(k): 2.430367         Tol: 0.00078291
k:5     Xn(k): 2.430466         Tol: 0.00000034
Final Solution: 2.430466
```

<hr>

### printVec()

> printing vector function

```c
void printVec(double* _vec, int _row);
```

**Parameter**

* _vec: the given arrays 
* _row: number of arrays

**Example**

```C
#include "myNP.h"

int main()
{
	int m = 12;
	double t[] = { -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5 };
	printVec(t, m);
}
```

**Output**

```c
Vector[0] = -1.0
Vector[1] = -0.5
Vector[2] = 0.0
Vector[3] = 0.5
Vector[4] = 1.0
Vector[5] = 1.5
Vector[6] = 2.0
Vector[7] = 2.5
Vector[8] = 3.0
Vector[9] = 3.5
Vector[10] = 4.0
Vector[11] = 4.5
```

<hr>

### gradient1D()

> 1st order differentiation using two-point central difference 3-2-3 method

```c
void gradient1D(double _x[], double _y[], double dydx[], int m);
```

**Parameter**

* _x[]: 1-D array for input x
* _y[]: 1-D array for input y
* dydx[]: 1-D array for output
* m: number of data

**Example**

```C
#include "myNP.h"

int main()
{
	int m = 12;
	double t[] = { -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5 };
	double y[] = { -3.632, -0.3935, 1, 0.6487, -1.282, -4.518, -8.611, -12.82, -15.91, -15.88, -9.402, 9.017 };
	double  dydt[12] = { 0 };
    
    gradient1D(t, y, dydt, m);
	printVec(dydt, m);
}
```

**Output**

```c
Vector[0] = 8.3
Vector[1] = 4.6
Vector[2] = 1.0
Vector[3] = -2.3
Vector[4] = -5.2
Vector[5] = -7.3
Vector[6] = -8.3
Vector[7] = -7.3
Vector[8] = -3.1
Vector[9] = 6.5
Vector[10] = 24.9
Vector[11] = 48.8
```

<hr>

### gradientFunc()

> 1st order differentiation using two-point central difference 3-2-3 method form the user-defined mathematic function

```c
void gradientFunc(double func(const double x), double x[], double dydx[], int m);
```

**Parameter**

* func: mathematic function
* x[]: 1-D array for input x
* dydx[]: 1-D array for output
* m: number of data

**Example**

```C
#include "myNP.h"

int main()
{
	int m = 21;
	double x[21] = { 0 };
	for (int i = 0; i < m; i++)
		x[i] = 0.2 * i;
	double dydx[21];
	gradientFunc(myFunc, x, dydx, m); // myFunc: F = x * x * x
	printVec(dydx, m);
}
```

**Output**

```c
Vector[0] = -0.1
Vector[1] = 0.2
Vector[2] = 0.5
Vector[3] = 1.1
Vector[4] = 2.0
Vector[5] = 3.0
Vector[6] = 4.4
Vector[7] = 5.9
Vector[8] = 7.7
Vector[9] = 9.8
Vector[10] = 12.0
Vector[11] = 14.6
Vector[12] = 17.3
Vector[13] = 20.3
Vector[14] = 23.6
Vector[15] = 27.0
Vector[16] = 30.8
Vector[17] = 34.7
Vector[18] = 38.9
Vector[19] = 43.4
Vector[20] = 47.9
```

<hr>

### acceleration()

> 2nd order differentiation using two-point central difference 4-3-4 method

```c
void acceleration(double x[], double y[], double dy2dx2[], int m);
```

**Parameter**

* x[]: 1-D array for input x
* y[]: 1-D array for input y
* dy2dx2[]: 1-D array for output
* m: number of data

**Example**

```C
#include "myNP.h"

int main()
{
	int m = 12;
	double t[] = { -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5 };
	double y[] = { -3.632, -0.3935, 1, 0.6487, -1.282, -4.518, -8.611, -12.82, -15.91, -15.88, -9.402, 9.017 };
	double  dydt[12] = { 0 };
    
    acceleration(t, y, dydt, m);
	printVec(dydt, m);
}
```

**Output**

```c
Vector[0] = -7.8
Vector[1] = -7.4
Vector[2] = -7.0
Vector[3] = -6.3
Vector[4] = -5.2
Vector[5] = -3.4
Vector[6] = -0.5
Vector[7] = 4.5
Vector[8] = 12.5
Vector[9] = 25.8
Vector[10] = 47.8
Vector[11] = 69.7
```

<hr>

### gradientFunck()

> 2nd order differentiation using two-point central difference 4-3-4 method form the user-defined mathematic function

```c
void gradientFunck(double func(const double x), double x[], double dy2dx2[], int m);
```

**Parameter**

* func: mathematic function
* x[]: 1-D array for input x
* dy2dx2[]: 1-D array for output
* m: number of data

**Example**

```C
#include "myNP.h"

int main()
{
	int m = 21;
	double x[21] = { 0 };
	for (int i = 0; i < m; i++)
		x[i] = 0.2 * i;
	double dy2dx2[21];
	gradientFunck(myFunc, x, dy2dx2, m);
	printVec(dy2dx2, m);
}
```

**Output**

```c
Vector[0] = -0.0
Vector[1] = 1.2
Vector[2] = 2.4
Vector[3] = 3.6
Vector[4] = 4.8
Vector[5] = 6.0
Vector[6] = 7.2
Vector[7] = 8.4
Vector[8] = 9.6
Vector[9] = 10.8
Vector[10] = 12.0
Vector[11] = 13.2
Vector[12] = 14.4
Vector[13] = 15.6
Vector[14] = 16.8
Vector[15] = 18.0
Vector[16] = 19.2
Vector[17] = 20.4
Vector[18] = 21.6
Vector[19] = 22.8
Vector[20] = 24.0
```

<hr>

### IntegrateRect()

> Integration rectangular method

```c
double IntegrateRect(double x[], double y[], int m);
```

**Parameter**

* x[]: 1-D array for input x
* y[]: 1-D array for input y
* m: number of data

**Example**

```C
#include "myNP.h"

int main()
{
	double x[] = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60 };
    double y[] = { 0, 3, 8, 20, 33, 42, 40, 48, 60, 12, 8, 4, 3 };
    int m = sizeof(x) / sizeof(x[0]);
    
    double I_rect = IntegrateRect(x, y, m);
	printf("I_rect  = %f\n\n", I_rect);
}
```

**Output** 

```c
I_rect  = 1390.000000
```

<hr>

### trapz()

> Integration trapezoidal method

```c
double trapz(double x[], double y[], int m);
```

**Parameter**

* x[]: 1-D array for input x
* y[]: 1-D array for input y
* m: number of data

**Example**

```C
#include "myNP.h"

int main()
{
	double x[] = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60 };
    double y[] = { 0, 3, 8, 20, 33, 42, 40, 48, 60, 12, 8, 4, 3 };
    int m = sizeof(x) / sizeof(x[0]);
    
    double I_trapz = trapz(x, y, m);
	printf("I_trapz = %f\n\n", I_trapz);
}
```

**Output**

```c
I_trapz = 1397.500000
```

<hr>

### simpson13()

> Integration Simpson's method

```c
double simpson13(double x[], double y[], int m);
```

**Parameter**

* x[]: 1-D array for input x
* y[]: 1-D array for input y
* m: number of data

**Example**

```C
#include "myNP.h"

int main()
{
	double X[] = { -3, -2.25, -1.5, -0.75, 0, 0.75,	1.5, 2.25, 3 };
	double Y[] = { 0, 2.1875, 3.75, 4.6875, 5, 4.6875, 3.75, 2.1875, 0 };
	int M = sizeof(X) / sizeof(X[0]);
    
	double I_simpson13 = simpson13(X, Y, M);
	printf("I_simpson13  = %f\n\n", I_simpson13);
}
```

**Output**

```c
I_simpson13  = 20.000000
```

<hr>

### integral()

> Integration Simpson's method using function

```c
double integral(double func(const double x), double a, double b, int n);
```

**Parameter**

* func: mathematic function
* a: end of the section
* b: start of the section
* n: the number of the interval

**Example**

```C
#include "myNP.h"

int main()
{
	double I_function = integral(myFunc, -1, 1, 12); // myFunc: F = sqrt(1 - x * x)
    printf("I_function  = %f\n\n", I_function);
}
```

**Output**

```c
I_function  = 1.555063
```

<hr>

### odeRK2()

> solve 1st order ODE equation using 2nd order Runge-Kutta method

```c
void odeRK2(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);
```

**Parameter**

* myfunc: mathematic function
* y[]: 1-D array for output 
* t0: start time
* tf: end time
* h: time intervals
* y0: initial value

**Example**

```C
#include "myNP.h"

int main()
{
	double a = 0, b = 0.1, h = 0.001, y0 = 0;
	double n = (b - a) / h;
	double* y;
    y = (double*)malloc(sizeof(double) * (n + 1));
    
    for (int i = 0; i < n + 1; i++)
		y[i] = 0;
   
    odeRK2(myfunc, y, a, b, h, y0);
    // myfunc: F = -T * y + 1 * T * Vm * cos(2 * PI * f * t), tau = 1, T = 1 / tau, f = 10, Vm = 1
	printVec(y, n + 1);
}
```

**Output**

```c
Vector[0] = 0.000000
Vector[1] = 0.000999
Vector[2] = 0.001992
Vector[3] = 0.002977
Vector[4] = 0.003949
....................
Vector[96] = -0.003940
Vector[97] = -0.002962
Vector[98] = -0.001972
Vector[99] = -0.000975
Vector[100] = 0.000024
```

<hr>

### odeRK3()

> solve 1st order ODE equation using 3rd order Runge-Kutta method

```c
void odeRK3(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);
```

**Parameter**

* myfunc: mathematic function
* y[]: 1-D array for output 
* t0: start time
* tf: end time
* h: time intervals
* y0: initial value

**Example**

```C
#include "myNP.h"

int main()
{
	double a = 0, b = 0.1, h = 0.001, y0 = 0;
	double n = (b - a) / h;
	double* y;
    y = (double*)malloc(sizeof(double) * (n + 1));
    
    for (int i = 0; i < n + 1; i++)
		y[i] = 0;
    
    odeRK3(myfunc, y, a, b, h, y0);
    // myfunc: F = -T * y + 1 * T * Vm * cos(2 * PI * f * t), tau = 1, T = 1 / tau, f = 10, Vm = 1
	printVec(y, n + 1);
}
```

**Output**

```c
Vector[0] = 0.000000
Vector[1] = 0.000999
Vector[2] = 0.001993
Vector[3] = 0.002978
Vector[4] = 0.003950
....................
Vector[96] = -0.003942
Vector[97] = -0.002963
Vector[98] = -0.001973
Vector[99] = -0.000976
Vector[100] = 0.000024
```

<hr>

### sysRK2()

> solve 2nd order ODE equation using 2rd order Runge-Kutta method

```c
void sysRK2(void func(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init);
```

**Parameter**

* func: mathematic function
* y1[]: 1-D array for output 
* y2[]: 1-D array for output 
* t0: start time
* tf: end time
* h: time intervals
* y1_init: y1 initial value
* y2_init: y2 initial value

**Example**

```C
#include "myNP.h"

int main()
{
	double a = 0, b = 1, h = 0.01, y_0 = 0, dydt_0 = 0.2;
	double n = (b - a) / h;
	double* y1;
	double* y2;
	y1 = (double*)malloc(sizeof(double) * (n + 1));
	y2 = (double*)malloc(sizeof(double) * (n + 1));
	
    for (int i = 0; i < n + 1; i++)
		y1[i] = 0, y2[i] = 0;
	
	sysRK2(mckfunc, y1, y2, a, b, h, y_0, dydt_0);
    // mckfunc: F0 = Y[1], F1 = (Fin - c * Y[1] - k * Y[0]) / m, m = 1, c = 7, k = 6.9, f = 5, Fin = 2 * cos(2 * PI * f * t)
	printVec(y1, n + 1);
	printVec(y2, n + 1);
}
```

**Output**

```c
Vector[0] = 0.000000
Vector[1] = 0.002010
Vector[2] = 0.004074
Vector[3] = 0.006168
Vector[4] = 0.008251
....................
Vector[96] = 0.012628
Vector[97] = 0.011922
Vector[98] = 0.011348
Vector[99] = 0.010947
Vector[100] = 0.010743

Vector[0] = 0.200000
Vector[1] = 0.205372
Vector[2] = 0.208315
Vector[3] = 0.207275
Vector[4] = 0.201164
....................
Vector[96] = -0.070196
Vector[97] = -0.057136
Vector[98] = -0.039952
Vector[99] = -0.020310
Vector[100] = -0.000119
```

<hr>

## NP myMatrix.h

```c
#include "myMatrix.h"
```
**List**

* createMat()
* freeMat()
* txt2Mat()
* printMat()
* addMat()
* zeros()
* ones()
* initMat()
* multMat()
* eye()
* transpose()
* copyMat()
* copyVal()
* gaussElim()
* gaussElim()
* LUdecompo()
* solveLU()
* fwdsub()
* backsub()
* invMat()
* magV()
* scamultMat()
* diag()
* eig()
* eigvec()
* linearFit()
* polyFit()
* expFit()
* arr2Mat()


<hr>

### createMat()

> Create Matrix with specified size

```c
Matrix createMat(int _rows, int _cols);
```
**Parameter**

* _rows: number of  the row
* _cols: number of the column

**Example**

```C
#include "myMatrix.h"

int main()
{
    int rows = 3;
	int cols = 3;
    Matrix Out = createMat(rows, cols);
    
    for (int i = 0; i < Out.rows; i++)
		for (int j = 0; j < Out.cols; j++)
			Out.at[i][j] = 0;
    
    printMat(Out, "Out");
}
```

**Output**

```c
Out =
       0.000000        0.000000        0.000000
       0.000000        0.000000        0.000000
       0.000000        0.000000        0.000000
```

<hr>

### freeMat()

> Free a memory allocated matrix

```c
void freeMat(Matrix _A);
```
**Parameter**

* _A: defined matrix

**Example**

```C
#include "myMatrix.h"

int main()
{
    int rows = 3;
	int cols = 3;
    Matrix Out = createMat(rows, cols);
    
    for (int i = 0; i < Out.rows; i++)
		for (int j = 0; j < Out.cols; j++)
			Out.at[i][j] = 0;
    
    printMat(Out, "Out");
    freeMat(Out);
}
```

**Output**

```c
Out =
       0.000000        0.000000        0.000000
       0.000000        0.000000        0.000000
       0.000000        0.000000        0.000000
```

<hr>

### txt2Mat()

> Create a matrix from a text file

```c
Matrix txt2Mat(std::string _filePath, std::string _fileName);
```
**Parameter**

* _filePath: the file path in local disk
* _fileName: the name of a file

**Example**

```C
#include "myMatrix.h"

int main()
{
    Matrix Out = txt2Mat(path, "prob1_matA"); // prob1_matA is zero matix
    
    printMat(Out, "Out");
}
```

**Output**

```c
Out =
       0.000000        0.000000        0.000000
       0.000000        0.000000        0.000000
       0.000000        0.000000        0.000000
```

<hr>

### printMat()

> Print matrix

```c
void printMat(Matrix _A, const char* _name);
```
**Parameter**

* _A: name of the matrix to be displayed 
* _name: name to display on screen

**Example**

```C
#include "myMatrix.h"

int main()
{
    Matrix Out = txt2Mat(path, "prob1_matA"); // prob1_matA is 3*3 zero matix
    
    printMat(Out, "Out");
}
```

**Output**

```c
Out =
       0.000000        0.000000        0.000000
       0.000000        0.000000        0.000000
       0.000000        0.000000        0.000000
```

<hr>

### addMat()

> Matrix addition

```c
Matrix addMat(Matrix _A, Matrix _B);
```
**Parameter**

* _A: 1st matrix
* _B: 2nd matrix

**Example**

```C
#include "myMatrix.h"

int main()
{
    Matrix A = txt2Mat(path, "prob1_matA"); // prob1_matA is 3*3 zero matix
    Matrix B = txt2Mat(path, "prob1_matB"); // prob1_matB is 3*3 eye matix
    Matrix Out = addMat(A, B);
    printMat(Out, "Out");
}
```

**Output**

```c
Out =
       1.000000        0.000000        0.000000
       0.000000        1.000000        0.000000
       0.000000        0.000000        1.000000
```

<hr>

### zeros()

> Create matrix of all zeros

```c
Matrix zeros(int _rows, int _cols);
```
**Parameter**

* _rows: number of rows
* _cols: number of columns

**Example**

```C
#include "myMatrix.h"

int main()
{
    Matrix Out = zeros(3, 3);
    
    printMat(Out, "Out");
}
```

**Output**

```c
Out =
       0.000000        0.000000        0.000000
       0.000000        0.000000        0.000000
       0.000000        0.000000        0.000000
```

<hr>

### ones()

> Create matrix of all ones

```c
Matrix ones(int _rows, int _cols);
```
**Parameter**

* _rows: number of rows
* _cols: number of columns

**Example**

```C
#include "myMatrix.h"

int main()
{
    Matrix Out = ones(3, 3);
    
    printMat(Out, "Out");
}
```

**Output**

```c
Out =
       1.000000        1.000000        1.000000
       1.000000        1.000000        1.000000
       1.000000        1.000000        1.000000
```

<hr>

### initMat()

> Initialization of Matrix elements

```c
void initMat(Matrix _A, double _val);
```
**Parameter**

* _A: defined matrix
* _val: value that wants to change

**Example**

```C
#include "myMatrix.h"

int main()
{
    int rows = 3;
	int cols = 3;
    int val = 3;
    Matrix Out = createMat(rows, cols);
    initMat(Out, 3);
    printMat(Out, "Out");
}
```

**Output**

```c
Out =
       3.000000        3.000000        3.000000
       3.000000        3.000000        3.000000
       3.000000        3.000000        3.000000
```

<hr>

### multMat()

> Multiply  matrix A and matrix B

```c
Matrix	multMat(Matrix _A, Matrix _B);
```
**Parameter**

* _A: 1st matrix
* _B: 2nd matrix

**Example**

```C
#include "myMatrix.h"

int main()
{
    Matrix A = txt2Mat(path, "prob1_matA"); // prob1_matA is 3*3 eye matix
    Matrix B = txt2Mat(path, "prob1_matB"); // prob1_matB is 3*3 eye matix
    Matrix Out = multMat(A, B);
    printMat(Out, "Out");
}
```

**Output**

```c
Out =
       1.000000        0.000000        0.000000
       0.000000        1.000000        0.000000
       0.000000        0.000000        1.000000
```

<hr>

### eye()

> Create identity  matrix

```c
Matrix	eye(int _rows, int _cols);
```
**Parameter**

* _rows: number of rows
* _cols: number of columns

**Example**

```C
#include "myMatrix.h"

int main()
{
    int rows = 3;
	int cols = 3;
    Matrix Out = eye(rows, cols);
    printMat(Out, "Out");
}
```

**Output**

```c
Out =
       1.000000        0.000000        0.000000
       0.000000        1.000000        0.000000
       0.000000        0.000000        1.000000
```

<hr>

### transpose()

> Create Transpose matrix

```c
Matrix	transpose(Matrix _A);
```
**Parameter**

* _A: defined matrix

**Example**

```C
#include "myMatrix.h"

int main()
{
    int rows = 3;
	int cols = 3;
    Matrix A = eye(rows, cols);
    Matrix Out = transpose(A);
    printMat(Out, "Out");
}
```

**Output**

```c
Out =
       1.000000        0.000000        0.000000
       0.000000        1.000000        0.000000
       0.000000        0.000000        1.000000
```

<hr>

### copyMat()

> Copy matrix

```c
Matrix copyMat(Matrix _A);
```
**Parameter**

* _A: defined matrix

**Example**

```C
#include "myMatrix.h"

int main()
{
    int rows = 3;
	int cols = 3;
    Matrix A = eye(rows, cols);
    Matrix Out = copyMat(A);
    printMat(Out, "Out");
}
```

**Output**

```c
Out =
       1.000000        0.000000        0.000000
       0.000000        1.000000        0.000000
       0.000000        0.000000        1.000000
```

<hr>

### copyVal()

> Copy matrix Elements from A to B

```c
void copyVal(Matrix _A, Matrix _B);
```
**Parameter**

* _A: defined matrix
* _B: matrix that wants to change the elements same as _A

**Example**

```C
#include "myMatrix.h"

int main()
{
    int rows = 3;
	int cols = 3;
    Matrix A = eye(rows, cols);
    Matrix Out = createMat(rows, cols);
    copyVal(A, Out);
    printMat(Out, "Out");
}
```

**Output**

```c
Out =
       1.000000        0.000000        0.000000
       0.000000        1.000000        0.000000
       0.000000        0.000000        1.000000
```

<hr>

### gaussElim()

> Gauss Elimination method without partial pivoting

```c
void gaussElim(Matrix A, Matrix b, Matrix U, Matrix d);
```
**Parameter**

* A: the input matrix A
* b: the input matrix b
* U: the output matrix A that passed elimination process
* d: the output matrix b that passed elimination process

**Example**

```C
#include "myMatrix.h"

int main()
{
	double arrA[] = { 0.5, -0.5, 0, 0, 0.87, 0.87,0, 0, 0,0,1,0,0,0,0,1 };	
	double arrB[] = { 0,-10,-5.75 * cos(PI / 3), 5.75 * sin(PI / 3) };
	Matrix matA = arr2Mat(arrA, 4, 4);
	Matrix vecb = arr2Mat(arrB, 4, 1);
	Matrix U = zeros(4, 4);
	Matrix d = zeros(4, 1);
    
    gaussElim(matA, vecb, U, d);
    printMat(U, "U");
    printMat(d, "d");
}
```

**Output**

```c
U =
       0.500000       -0.500000        0.000000        0.000000
       0.000000        1.740000        0.000000        0.000000
       0.000000        0.000000        1.000000        0.000000
       0.000000        0.000000        0.000000        1.000000
d =
       0.000000
     -10.000000
      -2.875000
       4.979646
```

<hr>

### gaussElim()

> Gauss elimination method with partial pivoting

```c
void gaussElim(Matrix A, Matrix b, Matrix U, Matrix d, Matrix P);
```
**Parameter**

* A: the input matrix A
* b: the input matrix b
* U: the output matrix A that passed elimination process
* d: the output matrix b that passed elimination process
* P: permutation matrix

**Example**

```C
#include "myMatrix.h"

int main()
{
	double arrA[] = { 0.5, -0.5, 0, 0, 0.87, 0.87,0, 0, 0,0,1,0,0,0,0,1 };	
	double arrB[] = { 0,-10,-5.75 * cos(PI / 3), 5.75 * sin(PI / 3) };
	Matrix matA = arr2Mat(arrA, 4, 4);
	Matrix vecb = arr2Mat(arrB, 4, 1);
    Matrix P = eye(4, 4);
	Matrix U = zeros(4, 4);
	Matrix d = zeros(4, 1);
    
    gaussElim(matA, vecb, U, d, P);
    printMat(U, "U");
    printMat(d, "d");
    printMat(P, "P");
}
```

**Output**

```c
U =
       0.870000        0.870000        0.000000        0.000000
       0.000000       -1.000000        0.000000        0.000000
       0.000000        0.000000        1.000000        0.000000
       0.000000        0.000000        0.000000        1.000000
d =
       0.000000
     -10.000000
      -2.875000
       4.979646
P =
       0.000000        1.000000        0.000000        0.000000
       1.000000        0.000000        0.000000        0.000000
       0.000000        0.000000        1.000000        0.000000
       0.000000        0.000000        0.000000        1.000000
```

<hr>

### LUdecompo()

> LU docomposition

```c
void LUdecompo(Matrix A, Matrix L, Matrix U, Matrix P);
```
**Parameter**

* A: the input matrix A
* L: the output matrix forms A = LU
* U: the output matrix forms A = LU
* P: permutation matrix

**Example**

```C
#include "myMatrix.h"

int main()
{
	double arrA[] = { 0.5, -0.5, 0, 0, 0.87, 0.87,0, 0, 0,0,1,0,0,0,0,1 };
	Matrix matA = arr2Mat(arrA, 4, 4);
    Matrix P = eye(4, 4);
    Matrix L = zeros(4, 4);
	Matrix U = zeros(4, 4);
    
    LUdecompo(matA, L, U, P);
    printMat(L, "L");
    printMat(U, "U");
    printMat(P, "P");
}
```

**Output**

```c
L =
       1.000000        0.000000        0.000000        0.000000
       0.574713        1.000000        0.000000        0.000000
       0.000000       -0.000000        1.000000        0.000000
       0.000000       -0.000000        0.000000        1.000000

U =
       0.870000        0.870000        0.000000        0.000000
       0.000000       -1.000000        0.000000        0.000000
       0.000000        0.000000        1.000000        0.000000
       0.000000        0.000000        0.000000        1.000000

P =
       0.000000        1.000000        0.000000        0.000000
       1.000000        0.000000        0.000000        0.000000
       0.000000        0.000000        1.000000        0.000000
       0.000000        0.000000        0.000000        1.000000
```

<hr>

### solveLU()

> Solving LU decomposition

```c
void solveLU(Matrix L, Matrix U, Matrix P, Matrix b, Matrix x);
```
**Parameter**

* L: the input matrix through LUdecompo()
* U: the input matrix through LUdecompo()
* P: permutation matrix through LUdecompo()
* b: the input matrix Ax = b
* x: the output matrix

**Example**

```C
#include "myMatrix.h"

int main()
{
	double arrA[] = { 0.5, -0.5, 0, 0, 0.87, 0.87,0, 0, 0,0,1,0,0,0,0,1 };
    double arrB[] = { 0,-10,-5.75 * cos(PI / 3), 5.75 * sin(PI / 3) };
	Matrix matA = arr2Mat(arrA, 4, 4);
    Matrix vecb = arr2Mat(arrB, 4, 1);
    Matrix vecx = zeros(4, 1);
    Matrix P = eye(4, 4);
    Matrix L = zeros(4, 4);
	Matrix U = zeros(4, 4);
    
    LUdecompo(matA, L, U, P);
    solveLU(L, U, P, vecb, vecx);
    printMat(vecx, "vecx");

}
```

**Output**

```c
vecx =
      -5.747126
      -5.747126
      -2.875000
       4.979646
```

<hr>

### fwdsub()

> Forward substitution

```c
void fwdsub(Matrix L, Matrix b, Matrix y);
```
**Parameter**

* L: the lower triangular matrix
* b: the input matrix Ly = b
* y: the output matrix

**Example**

```C
#include "myMatrix.h"

int main()
{
	double arrA[] = { 0.5, -0.5, 0, 0, 0.87, 0.87,0, 0, 0,0,1,0,0,0,0,1 };
    double arrB[] = { 0,-10,-5.75 * cos(PI / 3), 5.75 * sin(PI / 3) };
	Matrix matA = arr2Mat(arrA, 4, 4);
    Matrix vecb = arr2Mat(arrB, 4, 1);
    Matrix vecy = zeros(4, 1);
    Matrix P = eye(4, 4);
    Matrix L = zeros(4, 4);
	Matrix U = zeros(4, 4);
    
    LUdecompo(matA, L, U, P);
    Matrix _b = multMat(P, vecb);
    fwdsub(L, _b, vecy);
    printMat(vecy, "vecy");

}
```

**Output**

```c
vecy =
     -10.000000
       5.747126
      -2.875000
       4.979646
```

<hr>

### backsub()

> Back substitution

```c
void backsub(Matrix U, Matrix d, Matrix x);
```
**Parameter**

* U: the upper triangular matrix
* d: the input matrix Ux = d
* x: the output matrix

**Example**

```C
#include "myMatrix.h"

int main()
{
	double arrA[] = { 0.5, -0.5, 0, 0, 0.87, 0.87,0, 0, 0,0,1,0,0,0,0,1 };
    double arrB[] = { 0,-10,-5.75 * cos(PI / 3), 5.75 * sin(PI / 3) };
	Matrix matA = arr2Mat(arrA, 4, 4);
    Matrix vecb = arr2Mat(arrB, 4, 1);
    Matrix vecx = zeros(4, 1);
    Matrix vecy = zeros(4, 1);
    Matrix P = eye(4, 4);
    Matrix L = zeros(4, 4);
	Matrix U = zeros(4, 4);
    
    LUdecompo(matA, L, U, P);
    Matrix _b = multMat(P, vecb);
    fwdsub(L, _b, vecy);
    backsub(U, vecy, vecx);
    printMat(vecx, "vecx");

}
```

**Output**

```c
vecx =
      -5.747126
      -5.747126
      -2.875000
       4.979646
```

<hr>

### invMat()

> Inverse Matrix

```c
double invMat(Matrix A, Matrix Ainv);
```
**Parameter**

* A: the defined matrix
* Ainv: the output matrix (initially Ainv should be eye matrix)

**Example**

```C
#include "myMatrix.h"

int main()
{
	Matrix matA = txt2Mat(path, "prob1_matK"); // (75, -20, 0; -20, 35, -15; 0, -15, 15)
    int rows = matA.rows;
	int cols = matA.cols;
    Matrix matF = eye(rows, cols);
    
    invMat(matA, matF);
    printMat(matF, "matF");
}
```

**Output**

```c
matF =
       0.018182        0.018182        0.018182
       0.018182        0.068182        0.068182
       0.018182        0.068182        0.134848
```

<hr>

### magV()

> Magunitude of Vector

```c
double magV(Matrix V);
```
**Parameter**

* V: the input matrix

**Example**

```C
#include "myMatrix.h"

int main()
{
	double arrA[] = { 0.5, -0.5, 0, 0, 0.87, 0.87,0, 0, 0,0,1,0,0,0,0,1 };
	Matrix matA = arr2Mat(arrA, 4, 4);
	double val = magV(matA);
	printf("val = %f\n", val);
}
```

**Output**

```c
val = 1.003444
```

<hr>

### scamultMat()

> Multiply  matrix A and scalar val

```c
Matrix scamultMat(Matrix _A, double _val);
```
**Parameter**

* V: the input matrix

**Example**

```C
#include "myMatrix.h"

int main()
{
	double arrA[] = { 0.5, -0.5, 0, 0, 0.87, 0.87,0, 0, 0,0,1,0,0,0,0,1 };
	Matrix matA = arr2Mat(arrA, 4, 4);
	Matrix Out = scamultMat(matA, 10);
    
    printMat(Out,"Out");
}
```

**Output**

```c
Out =
       5.000000       -5.000000        0.000000        0.000000
       8.700000        8.700000        0.000000        0.000000
       0.000000        0.000000       10.000000        0.000000
       0.000000        0.000000        0.000000       10.000000
```

<hr>

### diag()

> vector, diagonal term of matrix A

```c
Matrix diag(Matrix _A);
```
**Parameter**

* _A: the input matrix

**Example**

```C
#include "myMatrix.h"

int main()
{
	double arrA[] = { 0.5, -0.5, 0, 0, 0.87, 0.87,0, 0, 0,0,1,0,0,0,0,1 };
	Matrix matA = arr2Mat(arrA, 4, 4);
	Matrix Out = diag(matA);
    
    printMat(Out,"Out");
}
```

**Output**

```c
Out =
       0.500000
       0.870000
       1.000000
       1.000000
```

<hr>

### eig()

> finding Eigenvalues

```c
Matrix eig(Matrix A);
```
**Parameter**

* A: the input matrix

**Example**

```C
#include "myMatrix.h"

int main()
{
	Matrix matA = txt2Mat(path, "prob1_matA"); (30, 15, 20; 15, 22, 26; 20, 26, 40)
    Matrix eigVals = eig(matA);
    printMat(eigVals, "eigVals");
}
```

**Output**

```c
eigVals =
      73.015425
      15.522467
       3.462108
```

<hr>

### eigvec()

> finding Eigenvetors

```c
Matrix eigvec(Matrix A);
```
**Parameter**

* A: the input matrix

**Example**

```C
#include "myMatrix.h"

int main()
{
	Matrix matA = txt2Mat(path, "prob1_matA"); (30, 15, 20; 15, 22, 26; 20, 26, 40)
    Matrix eigVecs = eigvec(matA);
    printMat(eigVecs, "eigVecs");
}
```

**Output**

```c

```

<hr>

### linearFit()

> Linear least square regression

```c
Matrix	linearFit(Matrix _X, Matrix _Y);
```
**Parameter**

* _X: the input matrix x data
* _Y: the input matrix y data

**Example**

```C
#include "myMatrix.h"

int main()
{
	int M = 6;
	double T_array[] = { 30, 40, 50, 60, 70, 80 };
    double P_array[] = { 1.05, 1.07, 1.09, 1.14, 1.17, 1.21 };
    Matrix T = arr2Mat(T_array, M, 1);
	Matrix P = arr2Mat(P_array, M, 1);
    Matrix z = linearFit(T, P);
    
	printMat(z, "z");
}
```

**Output**

```c
z =
       0.940952
       0.003286
```

<hr>

### polyFit()

> Curve fit with a high-order polynomial

```c
Matrix	polyFit(Matrix _X, Matrix _Y, int n);
```
**Parameter**

* _X: the input matrix x data
* _Y: the input matrix y data
* n: the number of the order

**Example**

```C
#include "myMatrix.h"

int main()
{
	int M = 16;
	double strain[] = { 0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.4, 4.8, 5.2, 5.6, 6.0 };
	double stress[] = { 0, 3, 4.5, 5.8, 5.9, 5.8, 6.2, 7.4, 9.6, 15.6, 20.7, 26.7, 31.1, 35.6, 39.3, 41.5 };
	Matrix Stra = arr2Mat(strain, M, 1);
	Matrix Stre = arr2Mat(stress, M, 1);
	Matrix z = polyFit(Stra, Stre, 4);
    
	printMat(z, "z");
}
```

**Output**

```c
z =
      -0.274607
      12.877980
     -10.192668
       3.118549
      -0.264389
```

<hr>

### expFit()

> Curve fit with an exponential function

```c
Matrix	expFit(Matrix _X, Matrix _Y);
```
**Parameter**

* _X: the input matrix x data
* _Y: the input matrix y data

**Example**

```C
#include "myMatrix.h"

int main()
{
	int M = 15;
	double time[] = { 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30 };
	double Voltage[] = { 9.7, 8.1, 6.6, 5.1, 4.4, 3.7, 2.8, 2.4, 2.0, 1.6, 1.4, 1.1, 0.85, 0.69, 0.6 };
	Matrix t = arr2Mat(time, M, 1);
	Matrix V= arr2Mat(Voltage, M, 1);
	Matrix z = expFit(t, V);
    
	printMat(z, "z");
}
```

**Output**

```c
z =
      11.913118
      -0.100161
```

<hr>

### arr2Mat()

> Create a matrix from 1D-array

```c
Matrix	arr2Mat(double* _1Darray, int _rows, int _cols);
```
**Parameter**

* _1Darray: the input array
* _rows: number of rows 
* _cols: number of column

**Example**

```C
#include "myMatrix.h"

int main()
{
	int M = 6;
    double T_array[] = { 30, 40, 50, 60, 70, 80 };
    Matrix T = arr2Mat(T_array, M, 1);
	printMat(T, "T");
}
```

**Output**

```c
T =
      30.000000
      40.000000
      50.000000
      60.000000
      70.000000
      80.000000
```