//
//  problem3.cpp
//  computational physics
//
//  Created by Haoyan Huo on 5/12/15.
//  Copyright (c) 2015 Haoyan Huo. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <cmath>

const double G = 9.81;
const double DT = 0.01;
double B2_M = 4e-5;

void iter_func(double& x, double& y, double& vx, double& vy){
    x += vx * DT;
    y += vy * DT;
    
    double dvxi, dvyi, sqrted = sqrt(vx*vx+vy*vy);
    dvxi = - B2_M * vx * sqrted * DT;
    dvyi = (- B2_M * vy * sqrted - G) * DT;
    vx += dvxi;
    vy += dvyi;
}

void solve(double theta, std::ofstream& fout){
    double v0 = 700.0;
    double x, y, vx, vy, t;
    
    // initial condition
    x = 0;
    y = 0;
    vx = v0 * cos(theta);
    vy = v0 * sin(theta);
    t = 0;
    
    int steps = 0;
    std::cout<<"begin solving problem"<<std::endl;
    while (y >= 0.0) {
        fout<<t<<' '<<x<<' '<<y<<std::endl;
        iter_func(x, y, vx, vy);
        t += DT;
        steps += 1;
    }
    fout<<t<<' '<<x<<' '<<y<<std::endl;
    
    std::cout<<"problem solved after "<<steps<<" steps"<<std::endl;
    std::cout<<"t: "<<t<<" s"<<std::endl
        <<"xf: "<<x<<" m"<<std::endl
        <<"v: "<<sqrt(vx*vx+vy*vy)<<" m/s"<<std::endl
        <<"vx: "<<vx<<" m/s"<<std::endl<<"vy: "<<vy<<" m/s"<<std::endl<<std::endl;
}

int problem3_main(){
    std::ofstream fout1("problem3.30.txt");
    solve(30.0 / 180.0 * acos(-1), fout1);
    fout1.close();
    
    std::ofstream fout2("problem3.40.txt");
    solve(40.0 / 180.0 * acos(-1), fout2);
    fout2.close();
    
    std::ofstream fout3("problem3.50.txt");
    solve(50.0 / 180.0 * acos(-1), fout3);
    fout3.close();
    
    B2_M = 0.0;
    
    std::ofstream fout4("problem3-no-air.30.txt");
    solve(30.0 / 180.0 * acos(-1), fout4);
    fout4.close();
    
    std::ofstream fout5("problem3-no-air.40.txt");
    solve(40.0 / 180.0 * acos(-1), fout5);
    fout5.close();
    
    std::ofstream fout6("problem3-no-air.50.txt");
    solve(50.0 / 180.0 * acos(-1), fout6);
    fout6.close();
    return 0;
}