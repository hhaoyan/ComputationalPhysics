//
//  fdm_radiation.cpp
//  computational physics
//
//  Created by Haoyan Huo on 5/2/15.
//  Copyright (c) 2015 Haoyan Huo. All rights reserved.
//

#include <iostream>
#include <fstream>

void solve_fdm(double dt, std::ofstream & fout){
    // dN / dt = - N
    double N = 100.0;
    
    for (double t=0.0; t<=5.0; t+= dt) {
        fout<<t<<' '<<N<<std::endl;
        N -= dt * N;
    }
}

// this function solves the problem by using FDM
int fdm_radiation_main(){
    double dts[] = {0.4, 0.2, 0.1, 0.05};
    
    for (int i = 0; i<sizeof(dts)/sizeof(dts[0]); ++i) {
        char buf[100];
        sprintf(buf, "/Users/huohaoyan/Desktop/output%d.txt", i);
        std::ofstream fout(buf);
        solve_fdm(dts[i], fout);
        fout.close();
    }
    
    return 0;
}