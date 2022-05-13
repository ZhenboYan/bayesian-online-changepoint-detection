#include <iostream>
#include <vector>
#include <cmath>
// #include <chrono> // comment out in the h file in device

// void printvec(std::vector<float> v){
//     for (auto i: v)
//         std::cout << i << ' ';
// }

// void print_matrix(const std::vector<std::vector<float>> &matrix){
//     for(int i = 0; i<matrix.size();i++){
//         for(int j = 0; j<matrix[i].size();j++){
//             std::cout << matrix[i][j] << " ";
//         }    
//         std::cout << std::endl;
//     }
// }

// allocate at once no more dynamic allocation
static int end = 350;
static std::vector<float> data,maxes,selection;
static std::vector<float> alphaT,betaT,muT,kappaT;
static std::vector<std::vector<float> > matrix;
static std::vector<float> var,predprobs,h,one_h,r1,r2,assign1,assign2,betaT0,muT0;

void initialize_vectors(int len,float lambda =200){
    data = std::vector<float> (len,0);
    maxes = std::vector<float> (len+1,1);
    var = std::vector<float> (len,0);
    predprobs = std::vector<float> (len,0);
    h = std::vector<float> (len+1,1.0/lambda);
    one_h = std::vector<float> (len+1,1.0-(1.0/lambda));
    assign1 = std::vector<float> (len+1,0);
    assign2 = std::vector<float> (len+1,0);

    alphaT = std::vector<float> (len+1,0);
    betaT = std::vector<float> (len+1,0);
    muT = std::vector<float> (len+1,0);
    betaT0 = std::vector<float> (len+1,0);
    muT0 = std::vector<float> (len+1,0);
    kappaT = std::vector<float> (len+1,0);
    selection = std::vector<float> (len+1,0);

    for (float i = 0; i < len+1; i++){
        std::vector<float> row(len+1, 0);
        matrix.push_back(row);
    }
    matrix[0][0] = 1;
    alphaT[0] = 1; 
    betaT[0] = 1;
    muT[0] = 0;
    betaT0[0] = 1;
    muT0[0] = 0;
    kappaT[0] =1;
    
}

bool detection(float datapoint,int t) {
    if (t == end)
        exit(EXIT_SUCCESS);    
    
    data[t] = datapoint;
    // algorithm itself

    // var = cal_var(alphaT,betaT,kappaT);
    for (int i = 0; i < t+1;i++)
        var[i] = ((kappaT[i] + 1) * betaT[i])/(kappaT[i] * alphaT[i]);
    
    float part5,part6,part7,part8,n,add05,gl1,gl2,partc1,mul_pi,mul_var,partc2,c;
        
    // predprobs = studentpdf(data[t],muT,var,alphaT);
    for (int i=0;i<t+1;i++){
        n = alphaT[i] * 2;
        part5 = 1 / (n * var[i]);
        part6 = pow((muT[i] - data[t]),2);
        part7 = (part5 * part6)+1;
        part8 = pow(part7,(-(n + 1) / 2));

        add05 = alphaT[i] + 0.5;
        gl1 = lgamma(add05);
        gl2 = lgamma(alphaT[i]);
        partc1 = exp(gl1 - gl2);
        mul_pi = n * M_PI;
        mul_var = mul_pi * var[i];
        partc2 = pow(mul_var,-0.5);
        c = partc1 * partc2;
        predprobs[i] = c * part8;
    }

    // r1 = select(matrix,0,t,t);
    int s1 = 0;
    for (int i=0;i<=t;i++){
        selection[s1] = matrix[i][t];
        s1++;
    }

    for (int i=0;i<=t;i++){
        assign1[i] = selection[i] * predprobs[i] * one_h[i];
    }
    // assign1 = ele_mul3(r1,predprobs,one_h);
    
    for (int i=1;i<=t+1;i++)
        matrix[i][t+1] = assign1[i-1];
    
    float sum = 0;
    for (int i = 0; i < t+2;i++)
        sum += selection[i] * predprobs[i] * h[i];
    
    matrix[0][t+1] = sum; // sum_ele_mul3(r1,predprobs,h)
    
    // r2 = select(matrix,0,end,t+1);
    int s2 = 0;
    for (int i=0;i<=end;i++){
        selection[s2] = matrix[i][t+1];
        s2++;
    }

    float sum1 = 0;
    for (int i = 0; i <=end;i++)
        sum1 += selection[i];
            
    // assign2 = all_div(r2,sumvec(r2));        
    
    for (int i=0;i<=end;i++)
        matrix[i][t+1] = selection[i] / sum1;
        // matrix[i][t+1] = assign2[i];
    // std::vector<float> v;
    // std::vector<float>::const_iterator beta = v1.begin();
    // std::vector<float>::const_iterator kappa = v2.begin();
    // std::vector<float>::const_iterator mu = v3.begin();
    // v.push_back(1);
    // for (int i = 0; i < v1.size();++i){
    //     v.push_back(*beta + (*kappa * pow((x - *mu),2)) / (2*(*kappa+1)));
    //     beta++;
    //     kappa++;
    //     mu++;
    // }
    // return v;

    for (int i=0; i<=t;i++){
        betaT0[i+1] = betaT[i] + (kappaT[i] * pow((data[t] - muT[i]),2)) / (2*(kappaT[i]+1));
        muT0[i+1] = ((kappaT[i] * muT[i] + data[t]) / (kappaT[i] + 1));
    }

    for (int i=0; i<=t+1;i++){
        betaT[i] = betaT0[i];
        muT[i] = muT0[i];
    }

    // betaT = cal_beta(data[t],betaT,kappaT,muT);
    // muT = cal_mu(data[t],muT,kappaT);
    kappaT[t+1] = kappaT[t] + 1;
    alphaT[t+1] = alphaT[t] + 0.5;
        
    int s3 = 0;
    for (int i=0;i<=end;i++){
        selection[s3] = matrix[i][t];
        s3++;
    }

    int index = 0;
    float max = 0;
    for (int i=0;i<t+1;i++){
        if (selection[i] > max){
            max = selection[i];
            index = i;
        }
    }
    maxes[t] = index;
    // maxes[t] = max_index(select(matrix,0,end,t));
    
    bool change = false;

    if (t > 0 and t != end){
        if (maxes[t-1] > maxes[t]){
            change = true;
        }
    }
    
    return change;
}