//
//  main.cpp
//  computational physics
//
//  Created by Haoyan Huo on 5/1/15.
//  Copyright (c) 2015 Haoyan Huo. All rights reserved.
//

#include <stdio.h>

int hermite_main();
int conjugate_gradient_main();
int spline_interpolate_main();
int fdm_radiation_main();
int sbm_linear_eq_main();
int problem3_main();
int solve_linear_eq_main();

int main(int argc, char** argv){
    return sbm_linear_eq_main();
}