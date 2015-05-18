//
//  solve_linear_eq.cpp
//  computational physics
//
//  Created by Haoyan Huo on 5/15/15.
//  Copyright (c) 2015 Haoyan Huo. All rights reserved.
//

#include <iostream>

class matrix{
    double* _data;
    unsigned _Nx, _Ny;
public:
    ~matrix(){
        delete[] _data;
    }
    
    matrix(double* init, unsigned Ny, unsigned Nx){
        _data = new double[Nx * Ny];
        memcpy(_data, init, sizeof(double) * Nx * Ny);
        _Nx = Nx;
        _Ny = Ny;
    }
    
    double& operator()(unsigned y, unsigned x){
        return _data[y*_Nx+x];
    }
    
    unsigned width(){
        return _Nx;
    }
    
    unsigned height(){
        return _Ny;
    }
    
    void print(){
        for (int i = 0; i<height(); ++i) {
            for (int j = 0; j<width(); ++j) {
                printf("%.5f ", operator()(i,j));
            }
            std::cout<<std::endl;
        }
    }
};

void solve_linear(matrix& AB){
    unsigned N = AB.height();
    
    for (unsigned i = 0; i<N; ++i) {
        // find the primary element
        {
            unsigned k = i;
            double ak = AB(k,i);
            for (unsigned j; j<N; ++j) {
                if (AB(j,i) > ak) {
                    ak = AB(j,i);
                    k = j;
                }
            }
            
            if(k!=i){
                for (unsigned j=0; j<=N; ++j) {
                    std::swap(AB(i,j), AB(k,j));
                }
            }
        }
        
        {
            for (unsigned j=i+1; j<N; ++j) {
                double l = AB(j,i) / AB(i,i);
                for (unsigned k = i; k<=N; ++k) {
                    AB(j,k) -= AB(i,k) * l;
                }
            }
            
        }
    }
    
    AB(N-1,N) /= AB(N-1,N-1);
    for (int j = N-2; j>=0; --j) {
        double sigma = 0;
        for (unsigned k = j+1; k<N; ++k) {
            sigma += AB(j,k) * AB(k,N);
        }
        
        AB(j,N) = (AB(j,N) - sigma) / AB(j,j);
    }
}

int solve_linear_eq_main(){
    double AB[] = {
        2.0/3.0, 2.0/5, 2.0/7, -2.0/9,
        2.0/5, 2.0/7, 2.0/9, -2.0/11,
        2.0/21, 2.0/15, 10.0/77, -14.0/117
    };
    matrix ab(AB, 3,4);
    
    solve_linear(ab);
    
    ab.print();
    
    return 0;
}