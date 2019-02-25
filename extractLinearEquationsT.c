/*
 * This is a c code for compiling a mex function, that given rotations for the problem, 
 * solves linearly the "d" of the polynomial system for multiple solutions
 **/
// Copyright (C) Yoni Kasten, Weizmann Institute, 2019
#include "mex.h"



void getMatrixIndex(size_t * matrixInd,size_t row,size_t column,size_t numRows,size_t numCols)
{
    *matrixInd =(column-1)*(numRows)+(row-1);
}


/**
 * Extract a linear system for multiple solutions of d to be represented by as Ad=b where
 * each solution for d has 4 entries, and A,b depends on the known solutions for rotation (represented as quaternions) and on the coefficients of the problem
 * inputs - coefficients of the original problem
 * outputA - the A matrix 
 * outputB - the b vector
 * q_2s - solutions for q_2 parameter
 * q_3s - solutions for q_3 parameter
 * q_4s - solutions for q_4 parameter
 */
void getAB(double *inputs, double *outputA,double *outputB,double *q_2s,double *q_3s,double *q_4s, size_t numSols)
{
    size_t i,j,k,count=0;
    size_t matrixInd;
    double q_2,q_3,q_4;
    size_t m,n;
    double	 t1_1, t1_2, t1_3, t2_1, t2_2, t2_3, t3_1, t3_2, t3_3, t4_1, t4_2, t4_3, t5_1, t5_2, t5_3, t6_1, t6_2, t6_3, p1_x, p1_y, p1_z, p1_x_t, p1_y_t, p2_x, p2_y, p2_z, p2_x_t, p2_y_t, p3_x, p3_y, p3_z, p3_x_t, p3_y_t, p4_x, p4_y, p4_z, p4_x_t, p4_y_t, p5_x, p5_y, p5_z, p5_x_t, p5_y_t, p6_x, p6_y, p6_z, p6_x_t, p6_y_t;
    m=7*numSols;
    n=4*numSols;
    
    t1_1=inputs[0];
    t1_2 = inputs[1];
    t1_3 =inputs[2];
    t2_1 = inputs[3];
    t2_2 = inputs[4];
    t2_3 = inputs[5];
    t3_1 = inputs[6];
    t3_2 = inputs[7];
    t3_3 = inputs[8];
    t4_1 = inputs[9];
    t4_2 = inputs[10];
    t4_3 = inputs[11];
    t5_1 = inputs[12];
    t5_2 = inputs[13];
    t5_3 = inputs[14];
    t6_1 =inputs[15];
    t6_2 = inputs[16];
    t6_3 = inputs[17];
    p1_x = inputs[18];
    p1_y =inputs[19];
    p1_z = inputs[20];
    p1_x_t = inputs[21];
    p1_y_t = inputs[22];
    p2_x = inputs[23];
    p2_y = inputs[24];
    p2_z = inputs[25];
    p2_x_t = inputs[26];
    p2_y_t =inputs[27];
    p3_x = inputs[28];
    p3_y = inputs[29];
    p3_z = inputs[30];
    p3_x_t = inputs[31];
    p3_y_t = inputs[32];
    p4_x = inputs[33];
    p4_y = inputs[34];
    p4_z = inputs[35];
    p4_x_t = inputs[36];
    p4_y_t = inputs[37];
    p5_x = inputs[38];
    p5_y = inputs[39];
    p5_z = inputs[40];
    p5_x_t = inputs[41];
    p5_y_t = inputs[42];
    p6_x = inputs[43];
    p6_y = inputs[44];
    p6_z =inputs[45];
    p6_x_t = inputs[46];
    p6_y_t = inputs[47];
    

    k=0;
    for (count=0;count<numSols;count++)
    {
        q_2=q_2s[count];
        q_3=q_3s[count];
        q_4=q_4s[count];
 
        i=(count)*7+1;
        j=(count)*4+1;
        outputB[k]=-(p1_x*t1_2 - p1_y*t1_1 - p1_x*q_2*q_2*t1_2 + p1_y*q_2*q_2*t1_1 - p1_x*q_3*q_3*t1_2 + p1_y*q_3*q_3*t1_1 + p1_x*q_4*q_4*t1_2 - p1_y*q_4*q_4*t1_1 - p1_x*p1_y_t*t1_3 + p1_y*p1_x_t*t1_3 - p1_z*p1_x_t*t1_2 + p1_z*p1_y_t*t1_1 + 2*p1_x*q_2*t1_3 - 2*p1_z*q_2*t1_1 + 2*p1_y*q_3*t1_3 - 2*p1_z*q_3*t1_2 - 2*p1_x*p1_x_t*q_3*t1_2 + 2*p1_x*p1_y_t*q_2*t1_2 + 2*p1_y*p1_x_t*q_3*t1_1 - 2*p1_y*p1_y_t*q_2*t1_1 - 2*p1_x*p1_x_t*q_4*t1_3 + 2*p1_z*p1_x_t*q_4*t1_1 - 2*p1_y*p1_y_t*q_4*t1_3 + 2*p1_z*p1_y_t*q_4*t1_2 - 2*p1_x*q_3*q_4*t1_3 + 2*p1_y*q_2*q_4*t1_3 - 2*p1_z*q_2*q_4*t1_2 + 2*p1_z*q_3*q_4*t1_1 + p1_x*p1_y_t*q_2*q_2*t1_3 + p1_y*p1_x_t*q_2*q_2*t1_3 - p1_z*p1_x_t*q_2*q_2*t1_2 - p1_z*p1_y_t*q_2*q_2*t1_1 - p1_x*p1_y_t*q_3*q_3*t1_3 - p1_y*p1_x_t*q_3*q_3*t1_3 + p1_z*p1_x_t*q_3*q_3*t1_2 + p1_z*p1_y_t*q_3*q_3*t1_1 + p1_x*p1_y_t*q_4*q_4*t1_3 - p1_y*p1_x_t*q_4*q_4*t1_3 + p1_z*p1_x_t*q_4*q_4*t1_2 - p1_z*p1_y_t*q_4*q_4*t1_1 - 2*p1_x*p1_x_t*q_2*q_3*t1_3 + 2*p1_x*p1_x_t*q_2*q_4*t1_2 - 2*p1_y*p1_x_t*q_2*q_4*t1_1 + 2*p1_z*p1_x_t*q_2*q_3*t1_1 + 2*p1_x*p1_y_t*q_3*q_4*t1_2 + 2*p1_y*p1_y_t*q_2*q_3*t1_3 - 2*p1_y*p1_y_t*q_3*q_4*t1_1 - 2*p1_z*p1_y_t*q_2*q_3*t1_2);
        k++;
        getMatrixIndex(&matrixInd, i, j, m, n);
        outputA[matrixInd]=p1_y*q_2 - p1_x*p1_x_t - p1_y*p1_y_t - p1_x*q_3 - p1_z + p1_x*p1_y_t*q_4 - p1_y*p1_x_t*q_4 + p1_z*p1_x_t*q_3 - p1_z*p1_y_t*q_2;
        j++;
        getMatrixIndex(&matrixInd, i, j, m, n);
        outputA[matrixInd]=p1_y - p1_z*p1_y_t - p1_x*q_4 + p1_z*q_2 - p1_x*p1_x_t*q_2 - p1_x*p1_y_t*q_3 - p1_y*p1_x_t*q_3 + p1_y*p1_y_t*q_2 - p1_z*p1_x_t*q_4;
        j++;
        getMatrixIndex(&matrixInd, i, j, m, n);
        outputA[matrixInd]=p1_z*p1_x_t - p1_x - p1_y*q_4 + p1_z*q_3 + p1_x*p1_x_t*q_3 - p1_x*p1_y_t*q_2 - p1_y*p1_x_t*q_2 - p1_y*p1_y_t*q_3 - p1_z*p1_y_t*q_4;
        j++;
        getMatrixIndex(&matrixInd, i, j, m, n);
        outputA[matrixInd]=p1_x*p1_y_t - p1_y*p1_x_t - p1_x*q_2 - p1_y*q_3 - p1_z*q_4 + p1_x*p1_x_t*q_4 - p1_z*p1_x_t*q_2 + p1_y*p1_y_t*q_4 - p1_z*p1_y_t*q_3;
        j++;
        i++;
        j=j-4;
        outputB[k]=p2_x*t2_2 - p2_y*t2_1 - p2_x*q_2*q_2*t2_2 + p2_y*q_2*q_2*t2_1 - p2_x*q_3*q_3*t2_2 + p2_y*q_3*q_3*t2_1 + p2_x*q_4*q_4*t2_2 - p2_y*q_4*q_4*t2_1 - p2_x*p2_y_t*t2_3 + p2_y*p2_x_t*t2_3 - p2_z*p2_x_t*t2_2 + p2_z*p2_y_t*t2_1 + 2*p2_x*q_2*t2_3 - 2*p2_z*q_2*t2_1 + 2*p2_y*q_3*t2_3 - 2*p2_z*q_3*t2_2 - 2*p2_x*p2_x_t*q_3*t2_2 + 2*p2_x*p2_y_t*q_2*t2_2 + 2*p2_y*p2_x_t*q_3*t2_1 - 2*p2_y*p2_y_t*q_2*t2_1 - 2*p2_x*p2_x_t*q_4*t2_3 + 2*p2_z*p2_x_t*q_4*t2_1 - 2*p2_y*p2_y_t*q_4*t2_3 + 2*p2_z*p2_y_t*q_4*t2_2 - 2*p2_x*q_3*q_4*t2_3 + 2*p2_y*q_2*q_4*t2_3 - 2*p2_z*q_2*q_4*t2_2 + 2*p2_z*q_3*q_4*t2_1 + p2_x*p2_y_t*q_2*q_2*t2_3 + p2_y*p2_x_t*q_2*q_2*t2_3 - p2_z*p2_x_t*q_2*q_2*t2_2 - p2_z*p2_y_t*q_2*q_2*t2_1 - p2_x*p2_y_t*q_3*q_3*t2_3 - p2_y*p2_x_t*q_3*q_3*t2_3 + p2_z*p2_x_t*q_3*q_3*t2_2 + p2_z*p2_y_t*q_3*q_3*t2_1 + p2_x*p2_y_t*q_4*q_4*t2_3 - p2_y*p2_x_t*q_4*q_4*t2_3 + p2_z*p2_x_t*q_4*q_4*t2_2 - p2_z*p2_y_t*q_4*q_4*t2_1 - 2*p2_x*p2_x_t*q_2*q_3*t2_3 + 2*p2_x*p2_x_t*q_2*q_4*t2_2 - 2*p2_y*p2_x_t*q_2*q_4*t2_1 + 2*p2_z*p2_x_t*q_2*q_3*t2_1 + 2*p2_x*p2_y_t*q_3*q_4*t2_2 + 2*p2_y*p2_y_t*q_2*q_3*t2_3 - 2*p2_y*p2_y_t*q_3*q_4*t2_1 - 2*p2_z*p2_y_t*q_2*q_3*t2_2;
        outputB[k]=-outputB[k];
        k++;
        getMatrixIndex(&matrixInd, i, j, m, n);
        outputA[matrixInd]=p2_y*q_2 - p2_x*p2_x_t - p2_y*p2_y_t - p2_x*q_3 - p2_z + p2_x*p2_y_t*q_4 - p2_y*p2_x_t*q_4 + p2_z*p2_x_t*q_3 - p2_z*p2_y_t*q_2;
        j++;
        getMatrixIndex(&matrixInd, i, j, m, n);
        outputA[matrixInd]=p2_y - p2_z*p2_y_t - p2_x*q_4 + p2_z*q_2 - p2_x*p2_x_t*q_2 - p2_x*p2_y_t*q_3 - p2_y*p2_x_t*q_3 + p2_y*p2_y_t*q_2 - p2_z*p2_x_t*q_4;
        j++;
        getMatrixIndex(&matrixInd, i, j, m, n);
        outputA[matrixInd]=p2_z*p2_x_t - p2_x - p2_y*q_4 + p2_z*q_3 + p2_x*p2_x_t*q_3 - p2_x*p2_y_t*q_2 - p2_y*p2_x_t*q_2 - p2_y*p2_y_t*q_3 - p2_z*p2_y_t*q_4;
        j++;
        getMatrixIndex(&matrixInd, i, j, m, n);
        outputA[matrixInd]=p2_x*p2_y_t - p2_y*p2_x_t - p2_x*q_2 - p2_y*q_3 - p2_z*q_4 + p2_x*p2_x_t*q_4 - p2_z*p2_x_t*q_2 + p2_y*p2_y_t*q_4 - p2_z*p2_y_t*q_3;
        j++;
        i++;
        j=j-4;
        outputB[k]=p3_x*t3_2 - p3_y*t3_1 - p3_x*q_2*q_2*t3_2 + p3_y*q_2*q_2*t3_1 - p3_x*q_3*q_3*t3_2 + p3_y*q_3*q_3*t3_1 + p3_x*q_4*q_4*t3_2 - p3_y*q_4*q_4*t3_1 - p3_x*p3_y_t*t3_3 + p3_y*p3_x_t*t3_3 - p3_z*p3_x_t*t3_2 + p3_z*p3_y_t*t3_1 + 2*p3_x*q_2*t3_3 - 2*p3_z*q_2*t3_1 + 2*p3_y*q_3*t3_3 - 2*p3_z*q_3*t3_2 - 2*p3_x*p3_x_t*q_3*t3_2 + 2*p3_x*p3_y_t*q_2*t3_2 + 2*p3_y*p3_x_t*q_3*t3_1 - 2*p3_y*p3_y_t*q_2*t3_1 - 2*p3_x*p3_x_t*q_4*t3_3 + 2*p3_z*p3_x_t*q_4*t3_1 - 2*p3_y*p3_y_t*q_4*t3_3 + 2*p3_z*p3_y_t*q_4*t3_2 - 2*p3_x*q_3*q_4*t3_3 + 2*p3_y*q_2*q_4*t3_3 - 2*p3_z*q_2*q_4*t3_2 + 2*p3_z*q_3*q_4*t3_1 + p3_x*p3_y_t*q_2*q_2*t3_3 + p3_y*p3_x_t*q_2*q_2*t3_3 - p3_z*p3_x_t*q_2*q_2*t3_2 - p3_z*p3_y_t*q_2*q_2*t3_1 - p3_x*p3_y_t*q_3*q_3*t3_3 - p3_y*p3_x_t*q_3*q_3*t3_3 + p3_z*p3_x_t*q_3*q_3*t3_2 + p3_z*p3_y_t*q_3*q_3*t3_1 + p3_x*p3_y_t*q_4*q_4*t3_3 - p3_y*p3_x_t*q_4*q_4*t3_3 + p3_z*p3_x_t*q_4*q_4*t3_2 - p3_z*p3_y_t*q_4*q_4*t3_1 - 2*p3_x*p3_x_t*q_2*q_3*t3_3 + 2*p3_x*p3_x_t*q_2*q_4*t3_2 - 2*p3_y*p3_x_t*q_2*q_4*t3_1 + 2*p3_z*p3_x_t*q_2*q_3*t3_1 + 2*p3_x*p3_y_t*q_3*q_4*t3_2 + 2*p3_y*p3_y_t*q_2*q_3*t3_3 - 2*p3_y*p3_y_t*q_3*q_4*t3_1 - 2*p3_z*p3_y_t*q_2*q_3*t3_2;
        outputB[k]=-outputB[k];
        k++;
        getMatrixIndex(&matrixInd, i, j, m, n);
        outputA[matrixInd]=p3_y*q_2 - p3_x*p3_x_t - p3_y*p3_y_t - p3_x*q_3 - p3_z + p3_x*p3_y_t*q_4 - p3_y*p3_x_t*q_4 + p3_z*p3_x_t*q_3 - p3_z*p3_y_t*q_2;
        j++;
        getMatrixIndex(&matrixInd, i, j, m, n);
        outputA[matrixInd]=p3_y - p3_z*p3_y_t - p3_x*q_4 + p3_z*q_2 - p3_x*p3_x_t*q_2 - p3_x*p3_y_t*q_3 - p3_y*p3_x_t*q_3 + p3_y*p3_y_t*q_2 - p3_z*p3_x_t*q_4;
        j++;
        getMatrixIndex(&matrixInd, i, j, m, n);
        outputA[matrixInd]=p3_z*p3_x_t - p3_x - p3_y*q_4 + p3_z*q_3 + p3_x*p3_x_t*q_3 - p3_x*p3_y_t*q_2 - p3_y*p3_x_t*q_2 - p3_y*p3_y_t*q_3 - p3_z*p3_y_t*q_4;
        j++;
        getMatrixIndex(&matrixInd, i, j, m, n);
        outputA[matrixInd]=p3_x*p3_y_t - p3_y*p3_x_t - p3_x*q_2 - p3_y*q_3 - p3_z*q_4 + p3_x*p3_x_t*q_4 - p3_z*p3_x_t*q_2 + p3_y*p3_y_t*q_4 - p3_z*p3_y_t*q_3;
        j++;
        i++;j=j-4;
        outputB[k]=p4_x*t4_2 - p4_y*t4_1 - p4_x*q_2*q_2*t4_2 + p4_y*q_2*q_2*t4_1 - p4_x*q_3*q_3*t4_2 + p4_y*q_3*q_3*t4_1 + p4_x*q_4*q_4*t4_2 - p4_y*q_4*q_4*t4_1 - p4_x*p4_y_t*t4_3 + p4_y*p4_x_t*t4_3 - p4_z*p4_x_t*t4_2 + p4_z*p4_y_t*t4_1 + 2*p4_x*q_2*t4_3 - 2*p4_z*q_2*t4_1 + 2*p4_y*q_3*t4_3 - 2*p4_z*q_3*t4_2 - 2*p4_x*p4_x_t*q_3*t4_2 + 2*p4_x*p4_y_t*q_2*t4_2 + 2*p4_y*p4_x_t*q_3*t4_1 - 2*p4_y*p4_y_t*q_2*t4_1 - 2*p4_x*p4_x_t*q_4*t4_3 + 2*p4_z*p4_x_t*q_4*t4_1 - 2*p4_y*p4_y_t*q_4*t4_3 + 2*p4_z*p4_y_t*q_4*t4_2 - 2*p4_x*q_3*q_4*t4_3 + 2*p4_y*q_2*q_4*t4_3 - 2*p4_z*q_2*q_4*t4_2 + 2*p4_z*q_3*q_4*t4_1 + p4_x*p4_y_t*q_2*q_2*t4_3 + p4_y*p4_x_t*q_2*q_2*t4_3 - p4_z*p4_x_t*q_2*q_2*t4_2 - p4_z*p4_y_t*q_2*q_2*t4_1 - p4_x*p4_y_t*q_3*q_3*t4_3 - p4_y*p4_x_t*q_3*q_3*t4_3 + p4_z*p4_x_t*q_3*q_3*t4_2 + p4_z*p4_y_t*q_3*q_3*t4_1 + p4_x*p4_y_t*q_4*q_4*t4_3 - p4_y*p4_x_t*q_4*q_4*t4_3 + p4_z*p4_x_t*q_4*q_4*t4_2 - p4_z*p4_y_t*q_4*q_4*t4_1 - 2*p4_x*p4_x_t*q_2*q_3*t4_3 + 2*p4_x*p4_x_t*q_2*q_4*t4_2 - 2*p4_y*p4_x_t*q_2*q_4*t4_1 + 2*p4_z*p4_x_t*q_2*q_3*t4_1 + 2*p4_x*p4_y_t*q_3*q_4*t4_2 + 2*p4_y*p4_y_t*q_2*q_3*t4_3 - 2*p4_y*p4_y_t*q_3*q_4*t4_1 - 2*p4_z*p4_y_t*q_2*q_3*t4_2;
        outputB[k]=-outputB[k];
        k++;
        getMatrixIndex(&matrixInd, i, j, m, n);
        outputA[matrixInd]=p4_y*q_2 - p4_x*p4_x_t - p4_y*p4_y_t - p4_x*q_3 - p4_z + p4_x*p4_y_t*q_4 - p4_y*p4_x_t*q_4 + p4_z*p4_x_t*q_3 - p4_z*p4_y_t*q_2;
        j++;
        getMatrixIndex(&matrixInd, i, j, m, n);
        outputA[matrixInd]=p4_y - p4_z*p4_y_t - p4_x*q_4 + p4_z*q_2 - p4_x*p4_x_t*q_2 - p4_x*p4_y_t*q_3 - p4_y*p4_x_t*q_3 + p4_y*p4_y_t*q_2 - p4_z*p4_x_t*q_4;
        j++;
        getMatrixIndex(&matrixInd, i, j, m, n);
        outputA[matrixInd]=p4_z*p4_x_t - p4_x - p4_y*q_4 + p4_z*q_3 + p4_x*p4_x_t*q_3 - p4_x*p4_y_t*q_2 - p4_y*p4_x_t*q_2 - p4_y*p4_y_t*q_3 - p4_z*p4_y_t*q_4;
        j++;
        getMatrixIndex(&matrixInd, i, j, m, n);
        outputA[matrixInd]=p4_x*p4_y_t - p4_y*p4_x_t - p4_x*q_2 - p4_y*q_3 - p4_z*q_4 + p4_x*p4_x_t*q_4 - p4_z*p4_x_t*q_2 + p4_y*p4_y_t*q_4 - p4_z*p4_y_t*q_3;
        j++;
        i++;j=j-4;
        outputB[k]=p5_x*t5_2 - p5_y*t5_1 - p5_x*q_2*q_2*t5_2 + p5_y*q_2*q_2*t5_1 - p5_x*q_3*q_3*t5_2 + p5_y*q_3*q_3*t5_1 + p5_x*q_4*q_4*t5_2 - p5_y*q_4*q_4*t5_1 - p5_x*p5_y_t*t5_3 + p5_y*p5_x_t*t5_3 - p5_z*p5_x_t*t5_2 + p5_z*p5_y_t*t5_1 + 2*p5_x*q_2*t5_3 - 2*p5_z*q_2*t5_1 + 2*p5_y*q_3*t5_3 - 2*p5_z*q_3*t5_2 - 2*p5_x*p5_x_t*q_3*t5_2 + 2*p5_x*p5_y_t*q_2*t5_2 + 2*p5_y*p5_x_t*q_3*t5_1 - 2*p5_y*p5_y_t*q_2*t5_1 - 2*p5_x*p5_x_t*q_4*t5_3 + 2*p5_z*p5_x_t*q_4*t5_1 - 2*p5_y*p5_y_t*q_4*t5_3 + 2*p5_z*p5_y_t*q_4*t5_2 - 2*p5_x*q_3*q_4*t5_3 + 2*p5_y*q_2*q_4*t5_3 - 2*p5_z*q_2*q_4*t5_2 + 2*p5_z*q_3*q_4*t5_1 + p5_x*p5_y_t*q_2*q_2*t5_3 + p5_y*p5_x_t*q_2*q_2*t5_3 - p5_z*p5_x_t*q_2*q_2*t5_2 - p5_z*p5_y_t*q_2*q_2*t5_1 - p5_x*p5_y_t*q_3*q_3*t5_3 - p5_y*p5_x_t*q_3*q_3*t5_3 + p5_z*p5_x_t*q_3*q_3*t5_2 + p5_z*p5_y_t*q_3*q_3*t5_1 + p5_x*p5_y_t*q_4*q_4*t5_3 - p5_y*p5_x_t*q_4*q_4*t5_3 + p5_z*p5_x_t*q_4*q_4*t5_2 - p5_z*p5_y_t*q_4*q_4*t5_1 - 2*p5_x*p5_x_t*q_2*q_3*t5_3 + 2*p5_x*p5_x_t*q_2*q_4*t5_2 - 2*p5_y*p5_x_t*q_2*q_4*t5_1 + 2*p5_z*p5_x_t*q_2*q_3*t5_1 + 2*p5_x*p5_y_t*q_3*q_4*t5_2 + 2*p5_y*p5_y_t*q_2*q_3*t5_3 - 2*p5_y*p5_y_t*q_3*q_4*t5_1 - 2*p5_z*p5_y_t*q_2*q_3*t5_2;
        outputB[k]=-outputB[k];
        k++;
        getMatrixIndex(&matrixInd, i, j, m, n);
        outputA[matrixInd]=p5_y*q_2 - p5_x*p5_x_t - p5_y*p5_y_t - p5_x*q_3 - p5_z + p5_x*p5_y_t*q_4 - p5_y*p5_x_t*q_4 + p5_z*p5_x_t*q_3 - p5_z*p5_y_t*q_2;
        j++;
        getMatrixIndex(&matrixInd, i, j, m, n);
        outputA[matrixInd]=p5_y - p5_z*p5_y_t - p5_x*q_4 + p5_z*q_2 - p5_x*p5_x_t*q_2 - p5_x*p5_y_t*q_3 - p5_y*p5_x_t*q_3 + p5_y*p5_y_t*q_2 - p5_z*p5_x_t*q_4;
        j++;
        getMatrixIndex(&matrixInd, i, j, m, n);
        outputA[matrixInd]=p5_z*p5_x_t - p5_x - p5_y*q_4 + p5_z*q_3 + p5_x*p5_x_t*q_3 - p5_x*p5_y_t*q_2 - p5_y*p5_x_t*q_2 - p5_y*p5_y_t*q_3 - p5_z*p5_y_t*q_4;
        j++;
        getMatrixIndex(&matrixInd, i, j, m, n);
        outputA[matrixInd]=p5_x*p5_y_t - p5_y*p5_x_t - p5_x*q_2 - p5_y*q_3 - p5_z*q_4 + p5_x*p5_x_t*q_4 - p5_z*p5_x_t*q_2 + p5_y*p5_y_t*q_4 - p5_z*p5_y_t*q_3;
        j++;
        i++;j=j-4;
        outputB[k]=p6_x*t6_2 - p6_y*t6_1 - p6_x*q_2*q_2*t6_2 + p6_y*q_2*q_2*t6_1 - p6_x*q_3*q_3*t6_2 + p6_y*q_3*q_3*t6_1 + p6_x*q_4*q_4*t6_2 - p6_y*q_4*q_4*t6_1 - p6_x*p6_y_t*t6_3 + p6_y*p6_x_t*t6_3 - p6_z*p6_x_t*t6_2 + p6_z*p6_y_t*t6_1 + 2*p6_x*q_2*t6_3 - 2*p6_z*q_2*t6_1 + 2*p6_y*q_3*t6_3 - 2*p6_z*q_3*t6_2 - 2*p6_x*p6_x_t*q_3*t6_2 + 2*p6_x*p6_y_t*q_2*t6_2 + 2*p6_y*p6_x_t*q_3*t6_1 - 2*p6_y*p6_y_t*q_2*t6_1 - 2*p6_x*p6_x_t*q_4*t6_3 + 2*p6_z*p6_x_t*q_4*t6_1 - 2*p6_y*p6_y_t*q_4*t6_3 + 2*p6_z*p6_y_t*q_4*t6_2 - 2*p6_x*q_3*q_4*t6_3 + 2*p6_y*q_2*q_4*t6_3 - 2*p6_z*q_2*q_4*t6_2 + 2*p6_z*q_3*q_4*t6_1 + p6_x*p6_y_t*q_2*q_2*t6_3 + p6_y*p6_x_t*q_2*q_2*t6_3 - p6_z*p6_x_t*q_2*q_2*t6_2 - p6_z*p6_y_t*q_2*q_2*t6_1 - p6_x*p6_y_t*q_3*q_3*t6_3 - p6_y*p6_x_t*q_3*q_3*t6_3 + p6_z*p6_x_t*q_3*q_3*t6_2 + p6_z*p6_y_t*q_3*q_3*t6_1 + p6_x*p6_y_t*q_4*q_4*t6_3 - p6_y*p6_x_t*q_4*q_4*t6_3 + p6_z*p6_x_t*q_4*q_4*t6_2 - p6_z*p6_y_t*q_4*q_4*t6_1 - 2*p6_x*p6_x_t*q_2*q_3*t6_3 + 2*p6_x*p6_x_t*q_2*q_4*t6_2 - 2*p6_y*p6_x_t*q_2*q_4*t6_1 + 2*p6_z*p6_x_t*q_2*q_3*t6_1 + 2*p6_x*p6_y_t*q_3*q_4*t6_2 + 2*p6_y*p6_y_t*q_2*q_3*t6_3 - 2*p6_y*p6_y_t*q_3*q_4*t6_1 - 2*p6_z*p6_y_t*q_2*q_3*t6_2;
        outputB[k]=-outputB[k];
        k++;
        getMatrixIndex(&matrixInd, i, j, m, n);
        outputA[matrixInd]=p6_y*q_2 - p6_x*p6_x_t - p6_y*p6_y_t - p6_x*q_3 - p6_z + p6_x*p6_y_t*q_4 - p6_y*p6_x_t*q_4 + p6_z*p6_x_t*q_3 - p6_z*p6_y_t*q_2;
        j++;
        getMatrixIndex(&matrixInd, i, j, m, n);
        outputA[matrixInd]=p6_y - p6_z*p6_y_t - p6_x*q_4 + p6_z*q_2 - p6_x*p6_x_t*q_2 - p6_x*p6_y_t*q_3 - p6_y*p6_x_t*q_3 + p6_y*p6_y_t*q_2 - p6_z*p6_x_t*q_4;
        j++;
        getMatrixIndex(&matrixInd, i, j, m, n);
        outputA[matrixInd]=p6_z*p6_x_t - p6_x - p6_y*q_4 + p6_z*q_3 + p6_x*p6_x_t*q_3 - p6_x*p6_y_t*q_2 - p6_y*p6_x_t*q_2 - p6_y*p6_y_t*q_3 - p6_z*p6_y_t*q_4;
        j++;
        getMatrixIndex(&matrixInd, i, j, m, n);
        outputA[matrixInd]=p6_x*p6_y_t - p6_y*p6_x_t - p6_x*q_2 - p6_y*q_3 - p6_z*q_4 + p6_x*p6_x_t*q_4 - p6_z*p6_x_t*q_2 + p6_y*p6_y_t*q_4 - p6_z*p6_y_t*q_3;
        j++;
        i++;j=j-4;
        outputB[k]=0;
        k++;
        getMatrixIndex(&matrixInd, i, j, m, n);
        outputA[matrixInd]=1;
        j++;
        getMatrixIndex(&matrixInd, i, j, m, n);
        outputA[matrixInd]=q_2;
        j++;
        getMatrixIndex(&matrixInd, i, j, m, n);
        outputA[matrixInd]=q_3;
        j++;
        getMatrixIndex(&matrixInd, i, j, m, n);
        outputA[matrixInd]=q_4;
        j++;
        i++;
    }
    
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    double *inputs; double *outputA;double *outputB;double *q_2s;double *q_3s;double *q_4s;
    size_t numSols,numRows,numCols;
    
    
    
    /*  check for proper number of arguments */
    /* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks you out of
     the MEX-file) */
    if(nrhs!=4)
        mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumInputs",
                "4 inputs required.");
    if(nlhs!=2)
        mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumOutputs",
                "Two outputs required.");
    
    
    /*  get the scalar input x */
    // x = mxGetScalar(prhs[0]);
    
    /*  create a pointer to the input matrix y */
    inputs = mxGetPr(prhs[0]);
    q_2s= mxGetPr(prhs[1]);
    q_3s= mxGetPr(prhs[2]);
    q_4s= mxGetPr(prhs[3]);
    /*  get the dimensions of the matrix input y */
    numSols = mxGetN(prhs[1]);
    numRows=numSols*7;
    numCols=numSols*4;
    
    /*  set the output pointer to the output matrix */
    plhs[0] = mxCreateDoubleMatrix( (mwSize)numRows, (mwSize)numCols, mxREAL);
    
    /*  create a C pointer to a copy of the output matrix */
    outputA = mxGetPr(plhs[0]);
    
    
    plhs[1] = mxCreateDoubleMatrix( (mwSize)numRows, (mwSize)1, mxREAL);
    
    outputB = mxGetPr(plhs[1]);
    /*  call the C subroutine */
    getAB(inputs,outputA,outputB,q_2s,q_3s,q_4s,numSols);
  
}
