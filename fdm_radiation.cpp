//
//  fdm_radiation.cpp
//  computational physics
//
//  Created by Haoyan Huo on 5/2/15.
//  Copyright (c) 2015 Haoyan Huo. All rights reserved.
//

#include <iostream>

void solve_fdm(double dt){
    // dN / dt = - N
    double N = 100.0;
    
    printf("ListLinePlot[{");
    for (double t=0.0; t<=5.0; t+= dt) {
        printf("%s{%.4f, %.4f}", t==0.0 ? "" : ",", t, N);
        
        N -= dt * N;
    }
    printf("}]");
}

// this function solves the problem by using FDM, and writes mathematica
// commands directly to stdout.
int fdm_radiation_main(){
    double dts[] = {0.4, 0.2, 0.1, 0.05};
    printf("Show[");
    for (int i = 0; i<sizeof(dts)/sizeof(dts[0]); ++i) {
        if(i!=0)
            printf(",");
        solve_fdm(dts[i]);
    }
    printf("]");
    
    return 0;
}