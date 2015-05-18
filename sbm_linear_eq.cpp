//
//  sbm_linear_eq.cpp
//  computational physics
//
//  Created by Haoyan Huo on 5/3/15.
//  Copyright (c) 2015 Haoyan Huo. All rights reserved.
//

#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>

const bool IN_DEBUG = false;

/*
 * A n*n(2m+1) symmetric band matrix is represented by a n*(m+1)
 * array, for example, matrix:
 * 
 *      a00 ...
 *      a10 a11 ...
 *      ... ... ...
 *      am  ... ... amm
 *      ... ... ... ... ...
 *                      an-1,m  an-1n-1
 *
 * will be represented by:
 * 
 *      a00 a11 a22 a33 ... an-1n-1
 *      a10 a21 a32 ... an-1,n-2 0
 *      a21 a32 ... an-1,n-3 0 0
 *      ...
 *      am0  ...
 *
 * so, anm = 
 *          c[n-m][m] n>=m
 *          c[m-n][n] n<=m
 *          0         abs(n-m) >= M
 */

class SBMMatrix{
    double* A;
    unsigned fN, fM;
    double fZero;
    
    // the zero elements in SBM matrix share one common data pointer,
    // in case this data pointer are accessed and changed by outer
    // functions, this class checks the data and prints error message
    // if it's changed.
    void checkZeroElements(){
        if(fZero != 0.0){
            printf("ERROR: zero representation changed in SBM!");
        }
    }
public:
    ~SBMMatrix(){
        delete A;
    }
    
    SBMMatrix(double* a, unsigned N, unsigned M){
        A = new double[N*(M+1)];
        memset(A, 0, sizeof(double)*(M+1)*N);
        fZero = 0.0;
        fN = N;
        fM = M;
        memcpy(A, a, sizeof(double)*(M+1)*N);
    }
    
    SBMMatrix(unsigned N, unsigned M){
        A = new double[N*(M+1)];
        memset(A, 0, sizeof(double)*(M+1)*N);
        fZero = 0.0;
        fN = N;
        fM = M;
    }
    
    double& operator()(unsigned i, unsigned j){
        checkZeroElements();
        
        i-=1;
        j-=1;
        
        if(i>=j && i-j<=fM){
            return A[fN*(i-j)+j];
        }else if(i<=j && j-i<=fM){
            return A[fN*(j-i)+i];
        }else{
            return fZero;
        }
    }
    
    unsigned N(){
        return fN;
    }
    
    unsigned M(){
        return fM;
    }
    
    void print(){
        for (int i = 1; i<=N(); ++i) {
            for (int j = 1; j<=N(); ++j) {
                printf("%f ", operator()(i,j));
            }
            printf("\n");
        }
        printf("\n");
    }
};

void solve_sbm_linear_eq(SBMMatrix& mat, std::vector<double>& B, std::vector<double>& X){
    SBMMatrix L(mat.N(), mat.M());
    
    for (int i = 1; i<=mat.N(); ++i) {
        for (int j = std::max<int>(1, i-mat.M()); j<=i; ++j) {
            L(i,j) = mat(i,j);
            for (int k = (i<=mat.M()+1 ? 1 : i-mat.M()); k<=j-1; ++k) {
                L(i,j) -= L(i,k) * L(j,k) / L(k,k);
            }
        }
    }
    
    std::vector<double> B1;
    B1.insert(B1.begin(), mat.N(), 0);
    
    for (int i = 1; i<=mat.N(); ++i) {
        B1[i-1] = B[i-1];
        for (int j = (i<=mat.M()+1 ? 1: i-mat.M()); j<=i-1; ++j) {
            B1[i-1] -= L(i,j)*B1[j-1]/L(j,j);
        }
    }
    
    X.clear();
    X.insert(X.begin(), mat.N(), 0);
    for (int i = mat.N(); i>=1; --i) {
        X[i-1] = B1[i-1];
        for (int j=i+1; j<= (i>mat.N()-mat.M()-1 ? mat.N() : i+mat.M()); ++j) {
            // eq.65 of lecture 2 is not right...
            X[i-1] -= L(j,i) * X[j-1];
        }
        X[i-1] /= L(i,i);
    }
}

int solve_special(int N){
    double* c = new double [N * 3];
    for (int i = 0; i<N; ++i) {
        if(i == 0 || i == N-1)
            c[i] = 5.0;
        else
            c[i] = 6.0;
        
        if(i != N-1){
            c[N+i] = 4.0;
        }
        
        if(i != N-1 && i != N-2){
            c[2 * N + i] = 1.0;
        }
    }
    
    SBMMatrix mat(c, N, 2);
    
    if(IN_DEBUG){
        printf("A:\n");
        mat.print();
    }
    
    std::vector<double> B;
    for (int i = 0; i<N; ++i) {
        if(i == 0 || i == N-1)
            B.push_back(60);
        else
            B.push_back(120);
    }
    
    std::vector<double> X;
    
    solve_sbm_linear_eq(mat, B, X);
    
    if(IN_DEBUG){
        printf("\nB:\n");
        for (int i = 0; i<N; ++i) {
            printf("%f ", B[i]);
        }
    }
    
    printf("\nX:\n");
    for (int i = 0; i<N; ++i) {
        printf("%f ", X[i]);
    }
    printf("\n");
    
    return 0;
}

int sbm_linear_eq_main(){
    printf("Solving N=100...\n");
    solve_special(100);
    printf("\n\nSolving N=10000...\n");
    solve_special(10000);
    
    return 0;
}