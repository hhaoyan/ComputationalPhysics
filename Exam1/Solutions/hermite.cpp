//
//  main.cpp
//  computational physics
//
//  Created by Haoyan Huo on 5/1/15.
//  Copyright (c) 2015 Haoyan Huo. All rights reserved.
//

#include <vector>
#include <iostream>


template<typename E>
class HermiteInterpolate{
    std::vector<E> mXi, mYi, mDi;
    unsigned length;
    
    E l_func(E x, unsigned j){
        E result = 1.0;
        for(int i = 0;i<length;++i){
            if(i == j)
                continue;
            result *= (x - mXi[i]) / (mXi[j] - mXi[i]);
        }
        return result;
    }
    
    E dl_func(E x, unsigned j){
        E result = 0.0;
        
        for(int k = 0;k<length;++k){
            if(k == j)
                continue;
            
            result += 1 / (mXi[j] - mXi[k]);
        }
        return result;
    }
    
    E alpha(E x, unsigned j){
        E result;
        
        E l_func_v = l_func(x, j);
        
        result = (1 - 2 * (x - mXi[j]) * dl_func(x, j))
        * l_func_v * l_func_v;
        return result;
    }
    
    E beta(E x, unsigned j){
        E result;
        
        E l_func_v = l_func(x, j);
        result = (x - mXi[j]) * l_func_v * l_func_v;
        
        return result;
    }
public:
    HermiteInterpolate(E* xi, E* yi, E* di, unsigned length){
        for(int i = 0;i<length;++i){
            mXi.push_back(xi[i]);
            mYi.push_back(yi[i]);
            mDi.push_back(di[i]);
        }
        
        this->length = length;
    }
    
    E evaluate(E x){
        E result = 0.0;
        
        for(int j = 0;j<length;++j){
            result += mYi[j] * alpha(x, j) + mDi[j] * beta(x, j);
        }
        
        return result;
    }
};

int hermite_main(){
    double xn[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5};
    double xyn[] = {0.0, 0.3, 0.4, 0.4, 0.3, 0.0};
    double xdn[] = {2.5, 0.6, 0.2,-0.2,-0.6,-2.5};
    HermiteInterpolate<double> hermite(xn, xyn, xdn, 6);
    
    double left = -1.5, right = 1.5;
    
    for(int i = 0;i<100;++i){
        double x = left + (right - left) * i / 100;
        std::cout<<x<<' '<<hermite.evaluate(x)<<std::endl;
    }
    
    return 0;
}

#ifdef ROOT_VERSION

int hermite(){
    TCanvas * c1 = new TCanvas();
    TGraph * g1 = new TGraph();
    
    double xn[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5};
    double xyn[] = {0.0, 0.3, 0.4, 0.4, 0.3, 0.0};
    double xdn[] = {2.5, 0.6, 0.2,-0.2,-0.6,-2.5};
    HermiteInterpolate<double> hermite(xn, xyn, xdn, 6);
    
    double left = 0, right = 0.5;
    
    for(int i = 0;i<100;++i){
        double x = left + (right - left) * i / 100;
        g1->SetPoint(i, x, hermite.evaluate(x));
    }
    g1->Draw();
    
    return 0;
}

#endif