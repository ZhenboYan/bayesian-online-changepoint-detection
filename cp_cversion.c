// #include <iostream>
// #include <>
#include <math.h>
#include <stdbool.h>
// void printvec(<float> v){
//     for (auto i: v)
//         cout << i << ' ';
// } 
const int end = 350;

void reset( int *array, int size) {
   memset(array,0,size * sizeof(*array));
}

float* ele_mul3(float *v1,float *v2,float *v3,int t){
    static float v[end];
    for (int i = 0; i < t;++i){
        v[t] = (v1[t] * v2[t] * v3[t]);
    }
    return v;
}

float sum_ele_mul3(const <float> &v1,const <float> &v2,const <float> &v3){
    float sum = 0;
    <float>::const_iterator it1 = v1.begin();
    <float>::const_iterator it2 = v2.begin();
    <float>::const_iterator it3 = v3.begin();
    for (int i = 0; i < v1.size();++i){
        sum += *it1 * *it2 * *it3;
        it1++;
        it2++;
        it3++;
    }
    return sum;
}

<float> all_mul(const <float> &v1, const float &num){
    <float> v;
    for (auto it = v1.begin(); it != v1.end(); ++it) 
        v.push_back(*it * num);
    return v;
}

<float> all_div(const <float> &v1, const float &num){
    <float> v;
    for (auto it = v1.begin(); it != v1.end(); ++it) 
        v.push_back(*it / num);
    return v;
}

<float> sub_all(const <float> &v1, const float &num){
    <float> v;
    for (auto it = v1.begin(); it != v1.end(); ++it) 
        v.push_back(num - *it);
    return v;
}

<float> hazard(const <float> &r, const float &lambda){
    if (lambda == 0)
        cout << "hazard failed division by zero";
    <float> ones(r.size(), 1);
    <float> v = all_mul(ones,1/lambda);
    return v;
}

<float> studentpdf(const float &x, const <float> &mu, const <float> &var, const <float> &alpha){
    <float> v;
    <float>::const_iterator va = var.begin();
    <float>::const_iterator al = alpha.begin();
    <float>::const_iterator m = mu.begin();
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

void print_matrix(const <<float>> &matrix){
    for(int i = 0; i<matrix.size();i++){
        for(int j = 0; j<matrix[i].size();j++){
            cout << matrix[i][j] << " ";
        }    
        cout << endl;
    }
}

<float> select_row_col(const <<float>> &matrix, int ls, int le, int r){
    <float> v;
    for (int i=ls;i<=le;++i)
        v.push_back(matrix[i][r]);
    return v;
}

float sumvec(const <float> &v1){
    float sum = 0;
    <float>::const_iterator it1 = v1.begin();
    for (int i = 0; i < v1.size();++i){
        sum += *it1;
        it1++;
    }
    return sum;
}

<float> cal_mu(float x, const <float> &v1,const <float> &v2){
    <float> v;
    <float>::const_iterator mu = v1.begin();
    <float>::const_iterator kappa = v2.begin();
    v.push_back(0);
    for (int i = 0; i < v1.size();++i){
        v.push_back((*kappa * *mu + x) / (*kappa + 1));
        mu++;
        kappa++;
    }
    return v;
}

<float> cal_beta(float x, const <float> &v1,const <float> &v2,const <float> &v3){
    <float> v;
    <float>::const_iterator beta = v1.begin();
    <float>::const_iterator kappa = v2.begin();
    <float>::const_iterator mu = v3.begin();
    v.push_back(1);
    for (int i = 0; i < v1.size();++i){
        v.push_back(*beta + (*kappa * pow((x - *mu),2)) / (2*(*kappa+1)));
        beta++;
        kappa++;
        mu++;
    }
    return v;
}

<float> cal_var(const <float> &v1,const <float> &v2,const <float> &v3){
    <float> v;
    <float>::const_iterator a = v1.begin();
    <float>::const_iterator b = v2.begin();
    <float>::const_iterator k = v3.begin();
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

int max_index(const <float> &v){
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


bool detection(float datapoint,int t) {
    static float data[350];
    static float kappaT[350]; // the past doesn't change 1
    static float alphaT[350]; // the past doesn't change 1
    static float maxes[351]; 
    static float matrix[351][351];
    static float ht[350]; // 0

    static float muT[350]; // the past changes 0 
    static float betaT[350];// the past changes 1

    // anything depends on these two will be updated because past chagnes
    if (t==0){
        for (int i = 0; i < end+1; i++){
            for(int j =0; j< end+1;++j){
                matrix[i][j] = 0;
            }
        }
        matrix[0][0] = 1;
    }
    float var,predprobs,h,one_h,r1,r2,assign1,assign2;
    
    data[t] = datapoint;
    // algorithm itself
    // var = cal_var(alphaT,betaT,kappaT);
    // predprobs = studentpdf(data[t],muT,var,alphaT);
    
    float ht[t] = {0};
    ht.push_back(0);
    h = hazard(ht,200);
    
    r1 = select_row_col(matrix,0,t,t);
    one_h = sub_all(h,1);        
    assign1 = ele_mul3(r1,predprobs,one_h);

    for (int i=1;i<=t+1;++i)
        matrix[i][t+1] = assign1[i-1];
    
    matrix[0][t+1] = sum_ele_mul3(r1,predprobs,h);
    
    r2 = select_row_col(matrix,0,end,t+1);
    assign2 = all_div(r2,sumvec(r2));        
    
    for (int i=0;i<=end;++i)
        matrix[i][t+1] = assign2[i];
        
    betaT = cal_beta(data[t],betaT,kappaT,muT);
    muT = cal_mu(data[t],muT,kappaT);
    kappaT.push_back(kappaT[t]+1);
    alphaT.push_back(alphaT[t]+0.5);
    
    maxes[t] = max_index(select_row_col(matrix,0,end,t));
    
    bool change = false;

    if (t > 0 and t != end){
        if (maxes[t-1] > maxes[t]){
            change = true;
        }
    }
    // if (t==349)
    //     printvec(predprobs);
    return change;
}