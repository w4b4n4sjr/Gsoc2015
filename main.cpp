#include <iostream>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <cstdlib>
#define COLS 3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
#define ROWS 3
using namespace std;
void gramSchmidt(double matb[ROWS][COLS], double matr[ROWS][COLS]);
void step2(double array1[ROWS][COLS], double d[COLS], double array3[ROWS][COLS], double array2[ROWS][COLS]);
void step4(double array[ROWS][COLS], double array1[ROWS][COLS], int k);
double scalar_product(double a[ROWS][COLS], double b[ROWS][COLS], int p, int q);
void swap(double matrix[ROWS][COLS], int r, int s);
int main(int argc, char** argv) {
    double u[ROWS][COLS] = {1,2,3,3,2,1,2,1,3};
    double v[ROWS][COLS] = {0,0,0,0,0,0,0,0,0};
    gramSchmidt(u,v);
    return 0;
}
void gramSchmidt(double matb[ROWS][COLS], double matr[ROWS][COLS])
{
    /* some iterators*/
    int i = 0, j = 0, k = 1;
    /* data structures for some storage*/
    double den[COLS];
    double coeffs[ROWS][COLS];
    /* intializations */
    for (i = 0; i < COLS; i++)
    {
        den[i] = 0;
    }
    for(i = 0; i < ROWS; i++ )
    {
        for(j = 0; j < COLS; j++)
        {
            coeffs[i][j] = 0;
        }
    }
    /* step1 of procedure, b1* = b1*/
    for(i = 0; i < ROWS; i++ )
    {
        matr[i][0] = matb[i][0];
    }
    /* calculation of B1*/
    den[0] = scalar_product(matr, matr, 0,0);
    /* call the routine step2*/
    step2(matb, den, coeffs,matr);
    /* step 3*/
    k = 2;
    /*call the routine step 4*/
    step4(coeffs, matb, k);
    /* step 5, calculation of u(i,j) */
    step2(matb, den, coeffs,matr);
    /* step 6*/
    while(k <= COLS)
    {
        step4(coeffs, matb, k);
        if(den[k-1] < (0.75 - coeffs[k-1][k-2]*coeffs[k-1][k-2])*den[k-2])
        {
            swap(matb, k-1, k-2);
            k = max(2, k-1);
        }
        else
        {
            k = k+1;
        }
    }
    /* output basis*/
    for( i = 0; i < ROWS; i++)
    {
        for(j = 0; j < COLS; j++)
        {
            cout << matb[i][j] << " ";
        }
        cout << endl;
    }
}

/* step 2 routine*/
void step2(double array1[ROWS][COLS], double d[COLS], double array3[ROWS][COLS], double array2[ROWS][COLS])
{
    /*some iterators*/
    int i = 0, j = 0, k = 0, p = 0;
    for(i = 2; i < COLS; i++)
    {
        for(j = 1; j < i-1; j++)
        {
            array3[i][j] = scalar_product(array1, array2, i, j)/d[j];
            /* b(i)* = b(i)* - u(i,j)*b(j)* */
            for(p = 0; p < ROWS; p++)
            {
                array2[p][i] = array2[p][i] - array3[i][j]*array2[p][j];
            }
        }
        d[i] = scalar_product(array2, array2, i, i);
    }
    
}
void step4(double array[ROWS][COLS], double array1[ROWS][COLS], int k)
{
    /*some iterators*/
    int i = 0, j = 0, l = 0, r;
    for(i = 2; i < ROWS; i++)
    {
        for(j = 1; j < i-1; j++)
        {
            if(sqrt(array1[i][j]*array1[i][j]) > 0.5)
            {
                for(l = 1; l < k - 1; l++ )
                {
                    r =  (int)(0.5 + array1[k][l]);
                    for(int s = 0; s < ROWS; s++)
                    {
                        array[s][k] = array[s][k] - r*array[s][l];
                    }
                    for(int t = 1; t < l - 1; t++)
                    {
                        array1[k][t] = array1[k][t] - r*array1[l][t];
                        array1[k][l] = array1[k][l] - r;
                    }
                }
            }
        }
    }
}
double scalar_product(double a[ROWS][COLS], double b[ROWS][COLS], int p, int q)
{
    double result = 0;
    for(int i = 0; i < ROWS; i++)
    {
        result += a[i][p]*b[i][q];
    }
    return result;
}
void swap(double matrix[ROWS][COLS], int r, int s)
{
    /* some temporary variables necessary for swapping*/
    double tmp = 0;
    for(int i = 0; i < ROWS; i ++)
    {
        tmp = matrix[i][r];
        matrix[i][r] = matrix[i][s];
        matrix[i][s] = tmp;
    }
}
int max(int x, int y)
{
    int maximum;
    if( x < y)
    {
        maximum = y;
    }
    else
    {
        maximum = x;
    }
}