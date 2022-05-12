#include <iostream>
#include <vector>
#include <cmath>
// #include <chrono> // comment out in the h file in device

void printvec(std::vector<float> v){
    for (auto i: v)
        std::cout << i << ' ';
} 

std::vector<float> ele_mul3(const std::vector<float> &v1,const std::vector<float> &v2,const std::vector<float> &v3){
    std::vector<float> v;
    std::vector<float>::const_iterator it1 = v1.begin();
    std::vector<float>::const_iterator it2 = v2.begin();
    std::vector<float>::const_iterator it3 = v3.begin();
    for (int i = 0; i < v1.size();++i){
        v.push_back(*it1 * *it2 * *it3);
        it1++;
        it2++;
        it3++;
    }
    return v;
}

float sum_ele_mul3(const std::vector<float> &v1,const std::vector<float> &v2,const std::vector<float> &v3){
    float sum = 0;
    std::vector<float>::const_iterator it1 = v1.begin();
    std::vector<float>::const_iterator it2 = v2.begin();
    std::vector<float>::const_iterator it3 = v3.begin();
    for (int i = 0; i < v1.size();++i){
        sum += *it1 * *it2 * *it3;
        it1++;
        it2++;
        it3++;
    }
    return sum;
}

std::vector<float> all_mul(const std::vector<float> &v1, const float &num){
    std::vector<float> v;
    for (auto it = v1.begin(); it != v1.end(); ++it) 
        v.push_back(*it * num);
    return v;
}

std::vector<float> all_div(const std::vector<float> &v1, const float &num){
    std::vector<float> v;
    for (auto it = v1.begin(); it != v1.end(); ++it) 
        v.push_back(*it / num);
    return v;
}

std::vector<float> sub_all(const std::vector<float> &v1, const float &num){
    std::vector<float> v;
    for (auto it = v1.begin(); it != v1.end(); ++it) 
        v.push_back(num - *it);
    return v;
}

std::vector<float> hazard(const std::vector<float> &r, const float &lambda){
    if (lambda == 0)
        std::cout << "hazard failed division by zero";
    std::vector<float> ones(r.size(), 1);
    std::vector<float> v = all_mul(ones,1/lambda);
    return v;
}

std::vector<float> studentpdf(const float &x, const std::vector<float> &mu, const std::vector<float> &var, const std::vector<float> &alpha){
    std::vector<float> v;
    std::vector<float>::const_iterator va = var.begin();
    std::vector<float>::const_iterator al = alpha.begin();
    std::vector<float>::const_iterator m = mu.begin();
    float part5,part6,part7,part8,n,add05,gl1,gl2,partc1,mul_pi,mul_var,partc2,c;
    for (int i = 0; i < var.size();++i){
        n = *al * 2;
        part5 = 1 / (n * *va);
        part6 = pow((*m - x),2);
        part7 = (part5 * part6)+1;
        part8 = pow(part7,(-(n + 1) / 2));

        add05 = *al + 0.5;
        gl1 = lgamma(add05);
        gl2 = lgamma(*al);
        partc1 = exp(gl1 - gl2);
        mul_pi = n * M_PI;
        mul_var = mul_pi * *va;
        partc2 = pow(mul_var,-0.5);
        c = partc1 * partc2;

        v.push_back(c * part8);
        va++;
        al++;
        m++;
    }
    return v;
}

void print_matrix(const std::vector<std::vector<float>> &matrix){
    for(int i = 0; i<matrix.size();i++){
        for(int j = 0; j<matrix[i].size();j++){
            std::cout << matrix[i][j] << " ";
        }    
        std::cout << std::endl;
    }
}

std::vector<float> select(const std::vector<std::vector<float>> &matrix, int ls, int le, int r){
    std::vector<float> v;
    for (int i=ls;i<=le;++i)
        v.push_back(matrix[i][r]);
    return v;
}

float sumvec(const std::vector<float> &v1){
    float sum = 0;
    std::vector<float>::const_iterator it1 = v1.begin();
    for (int i = 0; i < v1.size();++i){
        sum += *it1;
        it1++;
    }
    return sum;
}

float t_kappa(const std::vector<float> &kappa,const int &t){
    return(kappa[t] + 1);
}

float t_alpha(const std::vector<float> &alpha,const int &t){
    return (alpha[t] + 0.5);
}

std::vector<float> cal_mu(float x, const std::vector<float> &v1,const std::vector<float> &v2){
    std::vector<float> v;
    std::vector<float>::const_iterator mu = v1.begin();
    std::vector<float>::const_iterator kappa = v2.begin();
    v.push_back(0);
    for (int i = 0; i < v1.size();++i){
        v.push_back((*kappa * *mu + x) / (*kappa + 1));
        mu++;
        kappa++;
    }
    return v;
}

std::vector<float> cal_beta(float x, const std::vector<float> &v1,const std::vector<float> &v2,const std::vector<float> &v3){
    std::vector<float> v;
    std::vector<float>::const_iterator beta = v1.begin();
    std::vector<float>::const_iterator kappa = v2.begin();
    std::vector<float>::const_iterator mu = v3.begin();
    v.push_back(1);
    for (int i = 0; i < v1.size();++i){
        v.push_back(*beta + (*kappa * pow((x - *mu),2)) / (2*(*kappa+1)));
        beta++;
        kappa++;
        mu++;
    }
    return v;
}

std::vector<float> cal_var(const std::vector<float> &v1,const std::vector<float> &v2,const std::vector<float> &v3){
    std::vector<float> v;
    std::vector<float>::const_iterator a = v1.begin();
    std::vector<float>::const_iterator b = v2.begin();
    std::vector<float>::const_iterator k = v3.begin();
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

int max_index(const std::vector<float> &v){
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

// self update variables
static int end = 350;
static std::vector<float> data;
static std::vector<float> kappaT = {1}; // the past doesn't change
static std::vector<float> alphaT = {1}; // the past doesn't change
static std::vector<float> maxes(end+1, 1); 
static std::vector<std::vector<float>> matrix;
static std::vector<float> ht = {0};

bool detection(bool initialization, float datapoint,int t) {
    if (t == end){
        std::cout << "exceeds time limit" << std::endl;
        exit(EXIT_SUCCESS);    
    }

    // anything depends on these two will be updated because past chagnes
    std::vector<float> muT = {0}; // the past changes
    std::vector<float> betaT  = {1};// the past changes

    if (initialization){
        for (float i = 0; i < end+1; i++){
            std::vector<float> row(end+1, 0);
            matrix.push_back(row);
        }
        matrix[0][0] = 1;
    }

    std::vector<float> var,predprobs,h,one_h,r1,r2,assign1,assign2;
    data.push_back(datapoint);
    
    // algorithm itself
    var = cal_var(alphaT,betaT,kappaT);
    predprobs = studentpdf(data[t],muT,var,alphaT);
    
    ht.push_back(0);
    h = hazard(ht,200);
    
    r1 = select(matrix,0,t,t);
    one_h = sub_all(h,1);        
    assign1 = ele_mul3(r1,predprobs,one_h);

    for (int i=1;i<=t+1;++i)
        matrix[i][t+1] = assign1[i-1];
    
    matrix[0][t+1] = sum_ele_mul3(r1,predprobs,h);
    
    r2 = select(matrix,0,end,t+1);
    assign2 = all_div(r2,sumvec(r2));        
    
    for (int i=0;i<=end;++i)
        matrix[i][t+1] = assign2[i];
        
    betaT = cal_beta(data[t],betaT,kappaT,muT);
    muT = cal_mu(data[t],muT,kappaT);
    kappaT.push_back(t_kappa(kappaT,t));
    alphaT.push_back(t_alpha(alphaT,t));
        
    maxes[t] = max_index(select(matrix,0,end,t));
    
    bool change = false;

    if (t > 0 and t != end){
        if (maxes[t-1] > maxes[t]){
            change = true;
        }
    }

    return change;
}