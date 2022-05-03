#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

void printvec(vector<float> v){
    for (auto i: v)
        cout << i << ' ';
} 

vector<float> ele_op(vector<float>v1, char op,vector<float>v2){
    vector<float> v;
    vector<float>::iterator it1 = v1.begin();
    vector<float>::iterator it2 = v2.begin();
    float solu;
    for (int i = 0; i < v1.size();++i){
        if (op == '*')
            solu = *it1 * *it2;
        if (op == '/')
            solu = *it1 / *it2;
        if (op == '+')
            solu = *it1 + *it2;
        if (op == '-')
            solu = *it1 - *it2;
        v.push_back(solu);
        it1++;
        it2++;
    }
    return v;
}

vector<float> all_op(vector<float>v1,char op,float num,bool flip = false){
    vector<float> v;
    vector<float>::iterator it1 = v1.begin();
    float solu;
    for (int i = 0; i < v1.size();++i){
        if (!flip){
            if (op == '^')
                solu = pow(*it1,num);
            if (op == '+')
                solu = *it1 + num;
            if (op == '-')
                solu = *it1 - num;
            if (op == '*')
                solu = *it1 * num;
            if (op == '/')
                solu = *it1 / num;
        }
        else{
            if (op == '/')
                solu = num / *it1;
            if (op == '-')
                solu = num - *it1;
            if (op == '^')
                solu = pow(num,*it1);
        }

        if (op == 'g')
            solu = log(tgamma(*it1));
        if (op == 'e')
            solu = exp(*it1);

        v.push_back(solu);
        it1++;
    }
    return v;
}

vector<float> gammaln(vector<float>v1){
    vector<float> v;
    vector<float>::iterator it1 = v1.begin();
    float solu;
    for (int i = 0; i < v1.size();++i){
        solu = log(tgamma(*it1));
        v.push_back(solu);
        it1++;
    }
    return v;
}

vector<float> expvec(vector<float>v1){
    vector<float> v;
    vector<float>::iterator it1 = v1.begin();
    float solu;
    for (int i = 0; i < v1.size();++i){
        solu = exp(*it1);
        v.push_back(solu);
        it1++;
    }
    return v;
}

vector<float> neg(vector<float>v1){
    vector<float> v;
    vector<float>::iterator it1 = v1.begin();
    float solu;
    for (int i = 0; i < v1.size();++i){
        solu = *it1 * (-1);
        v.push_back(solu);
        it1++;
    }
    return v;
}

vector<float> hazard(vector<float> r, float lambda){
    if (lambda == 0){
        cout << "hazard failed division by zero";
        return NULL;
    }
    
    vector<float> ones = all_op(r,'^',0);
    vector<float> v = all_op(ones,'*',1/lambda);
    return v;
}
// function p = constant_hazard(r, lambda)
//   p = 1/lambda * ones(size(r));
vector<float> studentpdf(vector<float> x,vector<float> mu,vector<float> var,vector<float> nu){
    vector<float> c,p,x1,x2,x3,x4,x5,x6,x7,x8;
    v1 = gammaln(all_op(all_op(nu,'/',2),'+',0.5));
    v2 = gammaln(all_op(nu,'/',2));
    v3 = all_op(ele_op(all_op(nu,'*',M_PI),'*',var),'^',(-0.5));
    c  = ele_op(expvec(ele_op(v1,'-',v2)),'*',v3);
    
    v4 = all_op(ele_op(nu,'*',var),'/',1,true);
    v5 = all_op(ele_op(x,'-',mu),'^',2);
    v6 = neg(all_op(all_op(nu,'+',1),'/',2));
    
    v7 = all_op(ele_op(v4,'*',v5),'+',1);
    v8 = ele_op(v7,'^',v6);

    p = ele_op(c,'*',v8);
    return p;
    //             x1                  x2                    x3
    // exp(gammaln(nu/2 + 0.5) - gammaln(nu/2)) .* (nu.*pi.*var).^(-0.5);
    //                    x4              x5          x6                            
    //   p = c .* (1 + (1./(nu.*var)).*(x-mu).^2).^(-(nu+1)/2);
    //                             x7
    //                                    x8
}
// vector<int> nums(begin(numsArr), end(numsArr));

int main() {    
    
    float muT    = 0;
    float kappaT = 1;
    float alphaT = 1;
    float betaT  = 1;
    float a[] = {1,2,3,4,5,6};
    int count = sizeof(a)/sizeof(a[0]);
    cout << count << endl;
    int number;
    
    vector<float> tt(begin(a), end(a));
    vector<float> t2 = neg(tt);
    printvec(t2);
    return 0;
}
