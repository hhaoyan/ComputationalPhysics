//
//  spline_interpolate.cpp
//  computational physics
//
//  Created by Haoyan Huo on 5/2/15.
//  Copyright (c) 2015 Haoyan Huo. All rights reserved.
//

#include <vector>
#include <iostream>
#include <fstream>

template<typename E>
class SplineInterpolate{
    std::vector<E> fXn;
    std::vector<E> fHn; /* the distances between x[i+1] and x[i], N=fLength-1 */
    std::vector<E> fYn;
    unsigned fN; /* NOTE: this is *NOT* the length of fXn!!! fN = len(fXn) - 1*/
    
    std::vector<E> fMn;
    std::vector<E> fC1n, fC2n;
    
    void solveDiagLinearEq(std::vector<E> b, // b1 -> bN
                         std::vector<E>& a, // c1 -> cN-1
                         std::vector<E>& c, // a2 -> aN
                         std::vector<E>& x, // empty: x1 -> xN
                         std::vector<E> f, // f1 -> fN
                         unsigned N){
        x.clear();
        x.reserve(N);
        x.insert(x.begin(), N, 0);
        
        for (int i = 1; i<N; ++i) {
            E m = a[i] / b[i-1];
            b[i] -= m * c[i-1];
            f[i] -= m * f[i-1];
        }
        
        x[N-1] = f[N-1] / b[N-1];
        for (int i = N-2; i>=0; --i) {
            x[i] = (f[i] - c[i] * x[i+1]) / b[i];
        }
    }
    
    // computes the second derivative Mi of spline function at Xi
    void computeM(){
        std::vector<E> B;
        B.reserve(fN-1);
        
        {
            std::vector<E> delta1;
            delta1.reserve(fN);
            
            for (int i = 0; i<fN; ++i) {
                delta1.push_back((fYn[i+1]-fYn[i])/fHn[i]);
                if(i!=0){
                    B.push_back((delta1[i]-delta1[i-1])*6.0/(fHn[i-1]+fHn[i]));
                }
            }
        }
        
        std::vector<E> lambda, miu;
        lambda.reserve(fN-2);
        miu.reserve(fN-1);
        
        {
            lambda.push_back(fHn[1] / (fHn[0] + fHn[1]));
            for (int i = 2; i<fN-1; ++i) {
                E tmp = fHn[i-1] + fHn[i];
                lambda.push_back(fHn[i]/tmp);
                miu.push_back(fHn[i-1]/tmp);
            }
            miu.push_back(fHn[fN-2]/(fHn[fN-2]+fHn[fN-1]));
        }
        
        {
            fMn.reserve(fN+1);
            fMn.push_back(0);
            std::vector<E> diag(fN-1, 2);
            std::vector<E> result;
            
            solveDiagLinearEq(diag, miu, lambda, result, B, fN-1);
            
            fMn.insert(fMn.end(), result.begin(), result.end());
            fMn.push_back(0);
        }
    }
    
    void computeC1C2(){
        for (int i = 0; i<fN; ++i) {
            fC1n.push_back((fYn[i+1] - fYn[i]) / fHn[i] - fHn[i] * (fMn[i+1] - fMn[i]) / 6.0);
            fC2n.push_back((fYn[i] * fXn[i+1] - fYn[i+1] * fXn[i]) / fHn[i] - fHn[i] * (fXn[i+1] * fMn[i] - fXn[i] * fMn[i+1]) / 6.0);
        }
    }
    
    E pow3(E x){
        return x * x * x;
    }
public:
    SplineInterpolate(E* xi, E* yi, unsigned N){
        fXn.reserve(N);
        fYn.reserve(N);
        for (int i = 0; i<N; ++i) {
            fXn.push_back(xi[i]);
            fYn.push_back(yi[i]);
            fN = N - 1;
        }
        
        fHn.reserve(fN);
        for(int i = 0;i<fN;++i){
            fHn.push_back(fXn[i+1]-fXn[i]);
        }
        
        computeM();
        computeC1C2();
    }
    
    E evaluate(E x){
        unsigned idx = 0;
        for (; idx<fN; ++idx) {
            if(x >= fXn[idx] && x <= fXn[idx + 1])
                break;
        }
        
        return pow3(fXn[idx+1] - x) * fMn[idx] / fHn[idx] / 6 +
            pow3(x - fXn[idx]) * fMn[idx+1] / fHn[idx] / 6 +
            fC1n[idx] * x + fC2n[idx];
    }
};

// this function solves the problem 6. and prints the standard mathematica
// commands to stdout, which draws both the sample data point and the 101
// interpolated data point.
int spline_interpolate_main(){
    // set up x and y
    double a[] = {0, 3, 5, 7, 9, 11, 12, 13, 14, 15};
    double b[] = {0, 1.2, 1.7, 2.0, 2.1, 2.0, 1.8, 1.2, 1.0, 1.6};
    
    // build our interpolator, using double as the value type.
    SplineInterpolate<double> S(a, b, sizeof(a)/sizeof(a[0]));
    std::ofstream fout("output.txt");
    
    for (double x = a[0]; x<a[sizeof(a)/sizeof(a[0])-1]; x+=0.05) {
        fout<<x<<' '<<S.evaluate(x)<<std::endl;
    }
    fout.close();
    
    return 0;
}

