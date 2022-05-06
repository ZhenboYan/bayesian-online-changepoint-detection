#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

void printvec(vector<float> v){
    for (auto i: v)
        cout << i << ' ';
} 

vector<float> ele_mul(const vector<float> &v1,const vector<float> &v2){
    vector<float> v;
    vector<float>::const_iterator it1 = v1.begin();
    vector<float>::const_iterator it2 = v2.begin();
    for (int i = 0; i < v1.size();++i){
        v.push_back(*it1 * *it2);
        it1++;
        it2++;
    }
    return v;
}

vector<float> ele_mul3(const vector<float> &v1,const vector<float> &v2,const vector<float> &v3){
    vector<float> v;
    vector<float>::const_iterator it1 = v1.begin();
    vector<float>::const_iterator it2 = v2.begin();
    vector<float>::const_iterator it3 = v3.begin();
    for (int i = 0; i < v1.size();++i){
        v.push_back(*it1 * *it2 * *it3);
        it1++;
        it2++;
        it3++;
    }
    return v;
}

float sum_ele_mul3(const vector<float> &v1,const vector<float> &v2,const vector<float> &v3){
    float sum = 0;
    vector<float>::const_iterator it1 = v1.begin();
    vector<float>::const_iterator it2 = v2.begin();
    vector<float>::const_iterator it3 = v3.begin();
    for (int i = 0; i < v1.size();++i){
        sum += *it1 * *it2 * *it3;
        it1++;
        it2++;
        it3++;
    }
    return sum;
}

vector<float> ele_pow(const vector<float> &v1,const vector<float> &v2){
    vector<float> v;
    vector<float>::const_iterator it1 = v1.begin();
    vector<float>::const_iterator it2 = v2.begin();
    for (int i = 0; i < v1.size();++i){
        v.push_back(pow(*it1, *it2));
        it1++;
        it2++;
    }
    return v;
}

vector<float> ele_div(const vector<float> &v1, const vector<float> &v2){
    vector<float> v;
    vector<float>::const_iterator it1 = v1.begin();
    vector<float>::const_iterator it2 = v2.begin();
    for (int i = 0; i < v1.size();++i){
        v.push_back(*it1 / *it2);
        it1++;
        it2++;
    }
    return v;
}

vector<float> ele_add(const vector<float> &v1, const vector<float> &v2){
    vector<float> v;
    vector<float>::const_iterator it1 = v1.begin();
    vector<float>::const_iterator it2 = v2.begin();
    for (int i = 0; i < v1.size();++i){
        v.push_back(*it1 + *it2);
        it1++;
        it2++;
    }
    return v;
}

vector<float> ele_sub(const vector<float> &v1, const vector<float> &v2){
    vector<float> v;
    vector<float>::const_iterator it1 = v1.begin();
    vector<float>::const_iterator it2 = v2.begin();
    for (int i = 0; i < v1.size();++i){
        v.push_back(*it1 - *it2);
        it1++;
        it2++;
    }
    return v;
}

vector<float> all_add(const vector<float> &v1, const float &num){
    vector<float> v;
    for (auto it = v1.begin(); it != v1.end(); ++it) 
        v.push_back(*it + num);
    return v;
}

vector<float> all_sub(const vector<float> &v1, const float &num){
    vector<float> v;
    for (auto it = v1.begin(); it != v1.end(); ++it) 
        v.push_back(*it - num);
    return v;
}

vector<float> all_mul(const vector<float> &v1, const float &num){
    vector<float> v;
    for (auto it = v1.begin(); it != v1.end(); ++it) 
        v.push_back(*it * num);
    return v;
}

vector<float> all_div(const vector<float> &v1, const float &num){
    vector<float> v;
    for (auto it = v1.begin(); it != v1.end(); ++it) 
        v.push_back(*it / num);
    return v;
}

vector<float> all_pow(const vector<float> &v1, const float &num){
    vector<float> v;
    for (auto it = v1.begin(); it != v1.end(); ++it) 
        v.push_back(pow(*it,num));
    return v;
}

vector<float> div_all(const vector<float> &v1, const float &num){
    vector<float> v;
    for (auto it = v1.begin(); it != v1.end(); ++it) 
        v.push_back(num / *it);
    return v;
}

vector<float> sub_all(const vector<float> &v1, const float &num){
    vector<float> v;
    for (auto it = v1.begin(); it != v1.end(); ++it) 
        v.push_back(num - *it);
    return v;
}

vector<float> gammaln(const vector<float> &v1){
    vector<float> v;
    vector<float>::const_iterator it1 = v1.begin();
    float solu;
    for (int i = 0; i < v1.size();++i){
        solu = log(tgamma(*it1));
        v.push_back(solu);
        it1++;
    }
    return v;
}

vector<float> expvec(const vector<float> &v1){
    vector<float> v;
    vector<float>::const_iterator it1 = v1.begin();
    float solu;
    for (int i = 0; i < v1.size();++i){
        solu = exp(*it1);
        v.push_back(solu);
        it1++;
    }
    return v;
}

vector<float> neg(const vector<float> &v1){
    vector<float> v;
    vector<float>::const_iterator it1 = v1.begin();
    float solu;
    for (int i = 0; i < v1.size();++i){
        solu = *it1 * (-1);
        v.push_back(solu);
        it1++;
    }
    return v;
}

vector<float> hazard(const vector<float> &r, const float &lambda){
    if (lambda == 0)
        cout << "hazard failed division by zero";
    vector<float> ones(r.size(), 1);
    vector<float> v = all_mul(ones,1/lambda);
    return v;
}

vector<float> cal_c(const vector<float> &v2,const vector<float> &v3){
    vector<float> v;
    vector<float>::const_iterator var = v2.begin();
    vector<float>::const_iterator nu = v3.begin();
    float div2,add05,gl1,gl2,part3,mul_pi,mul_var,part4;
    for (int i = 0; i < v1.size();++i){
        div2 = *nu / 2;
        add05 = div + 0.5;
        gl1 = log(tgamma(add05));
        gl2 = log(tgamma(div2));
        part3 = exp(gl1 - gl2);
        mul_pi = *nu * M_PI;
        mul_var = mul_pi * *var;
        part4 = pow(mul_var,-0.5);
        v.push_back(part3*part4);
        var++;
        nu++;
    }
    return v;
}

vector<float> studentpdf(const float &x, const vector<float> &mu, const vector<float> &var, const vector<float> &nu){
    vector<float> v;
    cc = cal_c(var,nu);
    vector<float> v;
    vector<float>::const_iterator v = var.begin();
    vector<float>::const_iterator n = nu.begin();
    vector<float>::const_iterator m = mu.begin();
    vector<float>::const_iterator c = cc.begin();

    float part5,part6,part7,part8;
    for (int i = 0; i < v1.size();++i){
        part5 = 1 / (*n * *v);
        part6 = pow((*m - x),2);
        part7 = (part5 * part6)+1;
        part8 = pow(part7,(-(*n + 1) / 2));
        v.push_back(*c * part8);
        v++;
        n++;
        m++;
        c++;
    }
    return v;
}

// vector<float> studentpdf(const float &x, const vector<float> &mu, const vector<float> &var, const vector<float> &nu){
//     vector<float> c,p,v1,v2,v3,v4,v5,v6,v7,v8;
//     v1 = gammaln(all_add(all_div(nu,2),0.5));
//     v2 = gammaln(all_div(nu,2));
//     v3 = expvec(ele_sub(v1,v2));
//     v4 = all_pow(ele_mul(all_mul(nu,M_PI),var),(-0.5));
//     c = ele_mul (v3,v4);

//     v5 = div_all(ele_mul(nu,var),1);
//     v6 = all_pow(sub_all(mu,x),2);
//     v7 = all_add(ele_mul(v5,v6),1);
//     v8 = ele_pow(v7,neg(all_div(all_add(nu,1),2)));
//     p = ele_mul(c,v8);

//     return p;
// }
//             x1                  x2                    x3
// exp(gammaln(nu/2 + 0.5) - gammaln(nu/2)) .* (nu.*pi.*var).^(-0.5);
//                    x4              x5          x6                            
//   p = c .* (1 + (1./(nu.*var)).*(x-mu).^2).^(-(nu+1)/2);
//                             x7
//                                    x8

void print_matrix(const vector<vector<float>> &matrix){
    for(int i = 0; i<matrix.size();i++){
        for(int j = 0; j<matrix[i].size();j++){
            cout << matrix[i][j] << " ";
        }    
        cout << endl;
    }
}

vector<float> select(const vector<vector<float>> &matrix, int ls, int le, int r){
    vector<float> v;
    for (int i=ls;i<le;++i)
        v.push_back(matrix[i][r]);
    return v;
}

void matrix_assign(vector<vector<float>> &matrix, const vector<float> &v, const int &ls, const int &le, const int &r){
    int j = 0;
    for (int i=ls;i<le;++i){
        if (j >= v.size()){
            cout << "error length";
        }
        matrix[i][r] = v[j];
        j++;
    }
}

float sumvec(const vector<float> &v1){
    float sum = 0;
    for (int i = 0; i < v1.size();++i){
        sum += *it1;
        it1++;
    }
    return sum;
}

// vector<float> var1 = ele_mul(all_add(kappaT,1),betaT);
// vector<float> var2 = ele_mul(kappaT,alphaT);
// vector<float> var = ele_div(var1,var2);

vector<float> cal_var(const vector<float> &v1,const vector<float> &v2,const vector<float> &v3){
    vector<float> v;
    vector<float>::const_iterator a = v1.begin();
    vector<float>::const_iterator b = v2.begin();
    vector<float>::const_iterator k = v3.begin();
    float add_one,var1,var2;
    for (int i = 0; i < v1.size();++i){
        add_one = *k + 1;
        var1 = add_one * *b;
        var2 = *k * *a;
        v.push_back(var1/var2);
        a++;
        b++;
        k++;
    }
    return v;
}

vector<float> cal_mu(float x, const vector<float> &v1,const vector<float> &v2){
    vector<float> v;
    vector<float>::const_iterator mu = v1.begin();
    vector<float>::const_iterator kappa = v2.begin();
    v.push_back(0);
    for (int i = 0; i < v1.size();++i){
        v.push_back((*kappa * *mu + x) / (*kappa + 1));
        mu++;
        kappa++;
    }
    return v;
}

vector<float> cal_kappa(const vector<float> &v1){
    vector<float> v;
    vector<float>::const_iterator kappa = v1.begin();
    v.push_back(1);
    for (int i = 0; i < v1.size();++i){
        v.push_back(*kappa + 1);
        kappa++;
    }
    return v;
}

vector<float> cal_alpha(const vector<float> &v1){
    vector<float> v;
    vector<float>::const_iterator alpha = v1.begin();
    v.push_back(1);
    for (int i = 0; i < v1.size();++i){
        v.push_back(*alpha + 0.5);
        alpha++;
    }
    return v;
}

vector<float> cal_beta(float x, const vector<float> &v1,const vector<float> &v2,const vector<float> &v3){
    vector<float> v;
    vector<float>::const_iterator beta = v1.begin();
    vector<float>::const_iterator kappa = v2.begin();
    vector<float>::const_iterator mu = v3.begin();
    v.push_back(1);
    for (int i = 0; i < v1.size();++i){
        v.push_back(*beta + (*kappa * pow((x - *mu),2)) / (2*(*kappa+1)));
        beta++;
        kappa++;
        mu++;
    }
    return v;
}

//   betaT0  = [ beta0  ; betaT + (kappaT .*(X(t)-muT).^2)./(2*(kappaT+1)) ];
int max_index(const vector<float> &v){
    int index = 0;
    float max = 0;
    for (int i=0;i<v.size();++i){
        if (v[i] > max){
            max = v[i];
            index = i;
        }
    }
    return index;
}

int main() {
    vector<float> muT = {0}; //PASSED INTO FUNCTION
    vector<float> kappaT = {1}; //PASSED INTO FUNCTION
    vector<float> alphaT = {1}; //PASSED INTO FUNCTION
    vector<float> betaT  = {1}; //PASSED INTO FUNCTION
    
    // test examples
    float a[] = {1,2,3,4,5,6};
    vector<float> tt(begin(a), end(a));
    vector<float> t2 = all_add(tt,1);
    
    int n = 500;
    // 2D vector
    vector<vector<float>> matrix; //PASSED INTO FUNCTION
    vector<float> maxes(n+1, 0); //PASSED INTO FUNCTION

    for (float i = 0; i < n+1; i++){
        vector<float> row(n, 0);
        matrix.push_back(row);
    }
    
    vector<float> example; 
    
    // trial run
    for (int t = 0; t < n; ++t){
        // MATLAB
        // predprobs = studentpdf(X(t), muT, ...
        //                  betaT.*(kappaT+1)./(alphaT.*kappaT), ...
        //                  2 * alphaT);
        // C++ 
        vector<float> var = cal_var(alphaT,betaT,kappaT);
        vector<float> predprobs = studentpdf(example[t],muT,var,all_mul(alphaT,2));

        // MATLAB
        //   H = hazard_func([1:t]');
        // C++
        vector<float> ht(t+1,0);
        vector<float> h = hazard(ht,200);

        // MATLAB
        // R(2:t+1,t+1) = R(1:t,t) .* predprobs .* (1-H);
        // R(1,t+1) = sum( R(1:t,t) .* predprobs .* H );
        // R(:,t+1) = R(:,t+1) ./ sum(R(:,t+1));
        // C++
        vector<float> one_h = sub_all(h,1);
        vector<float> r1 = select(matrix,0,t,t);
        vector<float> r2 = select(matrix,0,matrix.size(),t);
        vector<float> assign1 = ele_mul3(r1,predprobs,one_h);
        matrix_assign(matrix,assign1,1,t+1,t+1);
        matrix[0][t+1] = sum_ele_mul3(r1,predprobs,h);
        vector<float> assign2 = all_div(r2,sum(r2));
        matrix_assign(matrix,assign2,0,matrix.size(),t+1);
        
        // MATLAB
        // muT0    = [ mu0    ; (kappaT.*muT + X(t)) ./ (kappaT+1) ];
        // kappaT0 = [ kappa0 ; kappaT + 1 ];
        // alphaT0 = [ alpha0 ; alphaT + 0.5 ];
        // betaT0  = [ beta0  ; betaT + (kappaT .*(X(t)-muT).^2)./(2*(kappaT+1)) ];
        // muT     = muT0;
        // kappaT  = kappaT0;
        // alphaT  = alphaT0;
        // betaT   = betaT0;
        // C++
        muT = cal_mu(example(t),muT,kappaT);
        kappaT = cal_kappa(kappaT);
        alphaT = cal_alpha(alphaT);
        betaT = cal_beta(example(t),betaT,kappaT,muT);
        
        // MATLAB
        // maxes(t) = find(R(:,t)==max(R(:,t)));
        // C++
        vector<float> r3 = select(matrix,0,matrix.size(),t);
        maxes(t) = max_index(r3);
    }
    
    return 0;
}
