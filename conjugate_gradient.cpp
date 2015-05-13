//
//  conjugate_gradient.cpp
//  computational physics
//
//  Created by Haoyan Huo on 5/1/15.
//  Copyright (c) 2015 Haoyan Huo. All rights reserved.
//

// codes here are fully optimized for solving the problem,
// however code convention is compromised.
// still a lot of stuff must be done to it...

#include <cstdio>
#include <cstring>
#include <cmath>

const bool IN_DEBUG = false;

void vector_times_A(double* v, double* ret, unsigned N){
    static double* result = NULL;
    if(result == NULL){
        // FIXME: here we leak N doubles, but it should be okay.
        result = new double[N];
    }
    
    for(int i = 0;i<N;++i){
        result[i] = v[i] * 3.0;
        if(i-1>=0)
            result[i] += v[i-1] * -1.0;
        if(i+1<=N-1)
            result[i] += v[i+1] * -1.0;
        if(N-1-i != i && N-1-i != i-1 && N-1-i != i+1)
            result[i] += v[N-1-i] * 0.5;
    }
    
    memcpy(ret, result, sizeof(double) * N);
}

void A_times_vector(double* v, double* ret, unsigned N){
    // v * A = A * v
    return vector_times_A(v, ret, N);
}

void vector_subtract(double* v1, double* v2, double* ret, unsigned N){
    for(int i = 0;i<N;++i){
        ret[i] = v1[i] - v2[i];
    }
}

double vector_abs(double* v, unsigned N){
    double result = 0.0;
    for(int i = 0;i<N;++i){
        result += v[i] * v[i];
    }
    
    return result;
}

double vector_times_vector(double* v1, double* v2, unsigned N){
    double result = 0.0;
    for(int i = 0;i<N;++i){
        result += v1[i] * v2[i];
    }
    
    return result;
}

int solve(double* x, double* b, unsigned N, double epsilon, bool init=true){
    double* r = new double[N];
    double* tmp = new double[N];
    double* r1 = new double[N];
    double* p = new double[N];
    
    int iter_step = 0;
    bool should_break = false;
    
    A_times_vector(x, tmp, N);
    vector_subtract(b, tmp, r, N);
    if(init)
        memset(x, 0, sizeof(double) * N);
    memcpy(p, r, sizeof(double) * N);
    
    if(IN_DEBUG){
        printf("step %.3d: ", iter_step);
        for(int i = 0;i<N;++i){
            printf("%.4f ", x[i]);
        }printf("\n");
    }
    
    while (!should_break) {
        double alpha = 0;
        
        // r(T) * r
        vector_times_A(p, tmp, N);
        alpha = vector_abs(r, N);
        if(alpha == 0.0)
            // ax = b
            break;
        alpha /= vector_times_vector(tmp, p, N);
        
        // update X
        should_break = true;
        for (int i = 0; i<N; ++i) {
            double update = alpha * p[i];
            if(std::abs(update) > epsilon)
                should_break = false;
            
            x[i] += update;
        }
        
        // r1 = r
        memcpy(r1, r, sizeof(double) * N);
        
        // r = r - alpha * A * p
        A_times_vector(p, tmp, N);
        for (int i = 0; i<N; ++i) {
            r[i] -= alpha * tmp[i];
        }
        
        // beta = - (r(T) * r) / (r1(T) * r1)
        double beta = - vector_abs(r, N) / vector_abs(r1, N);
        
        // p = r + beta * p
        for(int i = 0;i<N;++i){
            p[i] = r[i] + beta * p[i];
        }
        
        iter_step += 1;
        
        if(IN_DEBUG){
            printf("step %.3d: ", iter_step);
            for(int i = 0;i<N;++i){
                printf("%.4f ", x[i]);
            }printf("\n");
        }
    }
    
    delete[] r;
    delete[] tmp;
    delete[] r1;
    delete[] p;
    
    return iter_step;
}

int conjugate_gradient_main(){
    unsigned N = 10000;
    
    double* x = new double[N];
    double* b = new double[N];
    
    for(int i = 0;i<N;++i){
        if(i == 0 || i == N-1)
            b[i] = 2.5;
        else if(i == N/2)
            b[i] = 1.0;
        else if(i == N/2 - 1 && N % 2 == 0)
            b[i] = 1.0;
        else
            b[i] = 1.5;
    }
    
    printf("converged in %d steps\n", solve(x, b, N, 1e-6));
    printf("\nsolution in 3 digits precision:\n");
    for (int i = 0; i<N; ++i) {
        printf("%.3f ", x[i]);
    }printf("\n");
    
    return 0;
}