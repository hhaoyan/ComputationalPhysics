//
//  problem7.cpp
//  computational physics
//
//  Created by Haoyan Huo on 5/21/15.
//  Copyright (c) 2015 Haoyan Huo. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <vector>
#include <iostream>


const int N = 3;
const double diff_step = 1e-6;

int factorial(int N){
    return N >= 2 ? N * factorial(N-1) : 1;
}

double eval_function(double x){
    const double A = sqrt(pow(2.0/N, 3) * (factorial(N-1-1)) / 2.0 / N / factorial(N+1));
    double rho = 2.0 * x / N;
    
    return A * exp(-rho / 2.0) * rho * (-rho + (2*1)+1 + 1);
}

double D(double x){
    double ddxp = eval_function(x + diff_step), ddxn = eval_function(x - diff_step);
    return (ddxp - ddxn) / diff_step / 2.0;
}

double laplacian(double x){
    double r_2_D_foo_p = (x+diff_step) * (x+diff_step) * D(x+diff_step);
    double r_2_D_foo_n = (x-diff_step) * (x-diff_step) * D(x-diff_step);
    
    double partials = (r_2_D_foo_p - r_2_D_foo_n) / diff_step;
    return partials / 2.0 / x / x - 2.0 / x / x * eval_function(x);
}

double hamilton(double x){
    return -0.5 * laplacian(x) - 1.0 / x * eval_function(x);
}

double integral_body(double x){
    double f = eval_function(x);
    double hami = hamilton(x);
    return f * hami * x * x;
}

double integrator(double start, double end){
    int N=2;
    double result = 0.0;
    double h = (end - start) / N;
    int step = 1;
    
    // bootstrap
    for (int i = 0; i<N; ++i) {
        double f0 = integral_body(start + h * i);
        double f1 = integral_body(start + h * (i+1));
        result += h / 2 * (f0 + f1);
    }
    
    while (true) {
        std::cout<<"step "<< step<<": "<<result<<std::endl;
        
        double last_result = result;
        result /= 2.0;
        
        for (int i = 0; i<N; ++i) {
            result += h / 2.0 * integral_body(start + h * i + h / 2.0);
        }
        
        if(fabs(last_result - result)< 1e-6)
            break;
        
        h /= 2.0;
        N *= 2;
        step += 1;
    }
    
    return result;
}

int problem7_main(){
    std::cout<<"integrate..."<<std::endl;
    double result = integrator(1e-9, 60);
    std::cout<<"result: "<<result<<std::endl;
    return 0;
}