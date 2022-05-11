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
    for (int i = 0; i < v2.size();++i){
        div2 = *nu / 2;
        add05 = div2 + 0.5;
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
    vector<float> cc = cal_c(var,nu);
    vector<float>::const_iterator va = var.begin();
    vector<float>::const_iterator n = nu.begin();
    vector<float>::const_iterator m = mu.begin();
    vector<float>::const_iterator c = cc.begin();

    float part5,part6,part7,part8;
    for (int i = 0; i < var.size();++i){
        part5 = 1 / (*n * *va);
        part6 = pow((*m - x),2);
        part7 = (part5 * part6)+1;
        part8 = pow(part7,(-(*n + 1) / 2));
        v.push_back(*c * part8);
        va++;
        n++;
        m++;
        c++;
    }
    return v;
}

vector<float> printstudentpdf(const float &x, const vector<float> &mu, const vector<float> &var, const vector<float> &nu){
    vector<float> v;
    vector<float> cc = cal_c(var,nu);
    vector<float>::const_iterator va = var.begin();
    vector<float>::const_iterator n = nu.begin();
    vector<float>::const_iterator m = mu.begin();
    vector<float>::const_iterator c = cc.begin();

    float part5,part6,part7,part8;
    for (int i = 0; i < var.size();++i){
        part5 = 1 / (*n * *va);
        part6 = pow((*m - x),2);
        part7 = (part5 * part6)+1;
        part8 = pow(part7,(-(*n + 1) / 2));
        v.push_back(*c * part8);
        cout << *c << " " << part8 << endl;
        
        va++;
        n++;
        m++;
        c++;
    }
    return v;
}


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
    for (int i=ls;i<=le;++i)
        v.push_back(matrix[i][r]);
    return v;
}

float sumvec(const vector<float> &v1){
    float sum = 0;
    vector<float>::const_iterator it1 = v1.begin();
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
    
    int n = 350;
    // 2D vector
    vector<vector<float>> matrix; //PASSED INTO FUNCTION
    vector<float> maxes(n+1, 1); //PASSED INTO FUNCTION

    for (float i = 0; i < n+1; i++){
        vector<float> row(n+1, 0);
        matrix.push_back(row);
    }

    vector<float> x = {-2.89413585733584,-3.13757954265236,-3.89881375312536,-8.98925077892571,-2.62153756230441,-7.99188227755037,-4.28585422225016,-8.03078925444828,-0.172459089806954,-2.84331004725389,-2.98647890054335,-5.08194794999613,-6.72247984422000,-5.27199784309344,-0.372848281456204,-2.37505662334214,-6.40443539840911,-7.43055394448026,-4.65611002007197,-4.84900085123620,-3.68224323448457,-3.26802422199164,-4.94940491550879,-1.67841960736000,-7.21653858212221,-3.10166793730401,-4.86615787460141,-5.53451036441451,-4.47263576725563,-7.20642421133086,-6.17088140621246,-4.75948724401598,-0.744424549672364,-4.90369098692676,-5.18527501301778,-1.84297548686373,-5.12576114021551,-2.27393238765045,-0.382952202971286,-1.71766193088686,-3.77950789294723,-1.78203543645149,-3.83557584032453,-3.02263927059802,-7.37278642832170,-2.58287527688366,-5.00672185420727,-5.48845018646094,-4.27512727833949,-2.24664760594248,-5.78889975403171,-5.35547246509533,-2.14561441442919,-8.45760380594335,-6.17788576969700,-3.79035914154741,-5.84710102885098,-3.75923483318685,-4.70281512261768,-9.89363941366485,-1.35040152765250,-2.50195587656568,-3.49504269536162,-1.11321986090822,-0.547573753032270,-4.06352201757900,-0.370056301856919,-8.02534345556231,0.720428225083553,-0.601745589054245,-3.45352515299219,-5.08722908615257,-1.96169052732168,1.07263794478203,-5.67188780866544,-3.34178063869127,-9.04767286600750,-8.40988212099603,-2.04420939844551,-2.07365844174133,-2.01861326423309,-4.01254836458308,-1.25467239537884,-8.32108966780085,-7.89311805649111,-2.82291044330034,-4.12727676853064,-7.31407971763147,-5.75401674494276,-4.99427489041277,-5.16139981749749,-3.08451902679629,-1.52965882231012,-0.323653144493811,-2.62008632371280,-8.84414873160441,-9.26333731217383,-0.600561212975651,2.64882901116051,-1.45870252812687,-4.60312303864546,-3.52427213513064,-6.00066795377230,-3.50584577266910,1.36714449556860,-3.06879096131358,-6.56756516998363,-5.07417294546920,-2.03559441289650,-4.29233063879900,-4.85709737866494,-8.19472896706759,-7.23027219734924,-5.02132544014981,-7.84733298959818,-4.39108123415873,-4.38689724404473,-0.663622841066089,-5.81897202870696,-4.10342681348347,-4.95830814365757,-3.36330514858325,-5.92020816425332,-4.17135393328728,-5.86415929403152,-2.64029160000639,-1.89443592327620,-4.20118110998920,-5.24089358039068,-7.07564526964819,-4.44686396613911,-7.38479319309219,-3.81962653520402,-1.27450686794820,-0.138629944569164,-12.6152026068577,-3.15442608761371,-3.22404750696673,-5.25945227985581,-2.66158024214514,-6.10815898059910,-5.68361472587095,-0.0657106251193667,-0.633101373670186,-5.88219716760777,0.310943482949662,-1.64190800588205,-2.51929426250565,-2.56272158051646,-5.48462114855544,-6.22929460729310,-3.16627006043735,-4.70985306703890,-0.447900529505465,-8.60623428850987,-8.08421755616614,-0.520039187699381,-6.14355334233532,-5.16931798189571,-7.00923390164420,-3.24264409607031,-1.98665202999172,-4.43190644413263,-5.81791138164921,-4.96707234110365,-8.40475003388859,-3.63755026958384,-3.64385281964593,-6.49140455327287,-7.39460175455404,-6.87274466593925,-0.432371194712541,-5.77015693056854,-5.53284097200706,0.575482887146192,-5.00383490581456,-1.37273298542833,-4.48907176246101,-5.83889233477264,-4.30198522585837,-4.47662385451683,-3.83521826868708,-6.63918652091615,-8.01814186386279,-3.20596103186956,-4.22263990149217,-6.90121446142971,-6.41491737575562,-6.24676150836705,-3.25700801696393,-5.97351873104954,-3.12693581433078,-6.83139096645900,-1.78011763535227,-3.80373468317779,-3.69477435805760,-4.47273124937999,-0.612138658926004,-4.42266031587525,-2.88672485934955,-1.97521887218878,-5.12494467615146,-2.69713330384374,-4.14594048194887,-2.21341924428212,-4.57573685122901,-5.23937862819972,-2.32410091736290,-1.60990417876420,-4.56752129574926,-0.689374169366328,-0.984454181347688,-0.813945131957150,-2.51961501863286,-1.99637513868993,2.08601548070483,-1.10503308079890,-0.618146531198050,-3.74270527717942,-2.07976494468715,0.0427933136757805,-2.16206213475408,-3.30053646905225,-1.00821571159960,-0.0138612244597844,-1.50467426660837,-0.635172345012147,-0.781625959397180,1.31323732221030,-0.260319063672296,-1.73623401235730,-0.702531964695428,-1.15825036160384,2.13090829746537,-1.49750162909367,-0.801956768423833,-2.26956147983030,-0.601940725300159,-1.83896725158010,-1.81250493832757,-1.89024212488391,-4.33423470799806,-1.11573432027804,-2.49149070648842,0.622518571679779,-3.35828819696919,-2.35747322759157,-2.10003395378007,1.54055214103733,-1.36822356621177,-3.59640649421619,-0.490528330707827,-2.79076663281157,-1.24736779140742,-0.999100975367136,-3.73612192039640,-3.36051157032502,-0.860790736708343,-3.10846987466421,-2.86050896888865,-0.302467170064138,-4.04574277545687,-1.38194498476233,-1.94646624204525,-1.51799940952857,-1.91290032110976,-0.459834288900089,-0.641331211702238,-2.75882513902424,-1.58565260706514,0.0802915649331681,-1.83342953452294,-0.853386253645416,-2.38439002319934,-1.45708862120137,-1.06902546225508,-1.51881907220493,0.633468638886240,-1.94838047737155,-0.587843575201906,-4.69700163468728,1.80239288005130,-2.06182201854234,-3.70746020335692,-0.728094612364433,-0.795742352950786,-1.99003675171507,-0.355663675309980,1.60618644436810,-0.765277311486608,-2.17122823793244,0.388208344545117,0.687155323340297,-0.676290766310433,-2.22143227642808,-0.0810575295411029,-0.448531854767670,-1.07004990346418,-2.38379326026827,0.401581477666425,-1.02426488516619,-1.41143670715780,-2.14395305592448,1.25810059165313,-1.01609012163358,-0.126420557189876,-3.68010764789636,-2.38539656800316,-2.10986324495678,0.398542575486797,-2.95296252604373,-2.11933251043006,-1.47944854856645,-0.334735090598214,-1.99953461077805,-0.0335179484472321,-3.34336364596985,-1.63517415705056,-0.652780328407879,0.569092527577346,-0.950624764476918,-0.928845605969305,-1.17588848008877,0.557090246823577,-3.54175456786769,-1.23417383172420,0.653901097685375,-6.40313700242553,-2.32108126619026,-0.936787740528060,-2.13082150967901,-2.26476857101059,-2.09452386380253,-1.07620407209794,-1.59356029207629,-0.00546149138165242,-0.579590475156016,-0.517789782729861,-1.59904479405530,-0.295836721969105,-2.21061778828779,-1.74427923164077,-0.669726100283744,-1.77739718297530,1.77626106252061,-1.20013041543134,-1.05490968852755,-2.62602711203832,0.877476477402013,-2.20280380931230}; 
    // vector<float> x = {-2.76370360381049,0.0573163170589316,-0.706511493773318,-2.01534788966410,-0.407413074203505,-0.0836110224233372,-0.458329822538812,-1.25231490529184,-0.768216948938031,-0.213274650375708,-1.67322683109909,-2.92167216759239,0.283115035986327,-2.29526909498525,-2.17408024291826,-0.811688007665232,-0.00486050587537068,-0.669374961460806,-0.0428962374940307,1.00390573302637,-0.711984902194184,0.708635941581847,-0.462635148285531,0.327365703970926,-1.07209702483171,0.592907800405136,0.199497851067636,-1.23354129763229,-1.01661833257235,-1.42049814812400,-2.44637001152788,-0.346964431251829,0.536232709609980,-0.813033754841561,-0.258750555691758,1.17676683684095,-0.0702560653950652,1.43084411062057,0.203520465579950,-1.95354066335159,0.650823285836168,-2.14730164560560,-1.92582622870725,1.92837305480922,0.961273317196426,1.06965464638227,1.15535717472873,1.43754998539321,0.946659665823033,2.16535169719395,1.24992138532915,0.0121951634313825,1.46216496576516,2.36936047177808,-0.0994168196765182,0.0336673232559529,0.708441991165384,1.34394520485098,0.375698295085537,-0.0368962406737367,1.47737148592373,0.776520229206616,-0.268081671193695,0.0548807436368809,1.39704005675358,1.52406752121771,-0.154231362781193,0.0757285503874314,1.86171396125556,0.515471577198627,-0.202856354204293,0.819572073228988,-0.118774506813265,1.23663820989621,0.0993311383559840,0.163147172994370,1.97032855010959,0.881585570635071,0.0673744287914302,0.884725153379912,1.54601544281183,1.25481156867620,1.68773050068664,0.221786678106811,0.854047250626853,0.262645190529329,0.524470489024543,1.46734874855561,0.136140510128175,2.19420258648142,0.774399090288669,1.57428975290576,-1.03787462564839,-0.268346834828350,1.01573401034109,1.82861259499774,1.20861655441001,2.02357730206710,0.151522618287368,0.524909991537515,0.693419118737008,0.739382705605774,0.724909146457398,0.241754829992018,1.27524341530032,1.01037628341288,0.979783709254315,-0.668155199198708,0.935108190206481,0.0424434436118272,0.883196820752743,0.887216927343987,0.686529834447601,0.604303055918822,0.721088746232720,0.980959397450372,-0.208709667006959,0.386159558712483,0.355627029311358,0.230649461682905,1.90002113187085,2.31505395378497,0.544118367496148,1.04651905909166,0.759198239941882,2.01686431686251,1.12856252336368,1.53069140257490,0.309874078168937,1.39441892296008,-0.391792478849823,0.562310983940428,-0.696162322777339,-0.136707576193802,0.926635093928321,1.45152649074573,0.0154838540588851,1.84037564347030,-0.396516105220739,2.60624927022471,0.0304046516218504,0.208660963314265,1.03336397411501,-0.136993272097423,1.81955232392767,0.638859045689281,0.651169400186277,1.21760362338248,1.84904236400273,1.87366266140761,0.778473564719651,0.337909240777062,-0.0944212116181119,1.64920162532335,0.361906800039614,1.43124324639019,1.39286481494686,0.799707075331078,0.569407249188919,1.37604180420732,1.65215254186396,-0.235395384781440,0.383020574860179,0.505699817947355,-0.190869017362302,0.478833955219883,0.196853396518090,1.55682676304324,0.686443980876115,0.0270802673341642,0.970642228328129,1.13789733163280,0.274990112117435,1.39844342177188,1.85431768319935,-0.273178966989566,0.483200257319732,0.845948814382700,0.00603163591592526,1.52141906676214,0.745459653425324,0.0577962897163991,0.447592517664731,1.29569202744024,-0.348611458563669,0.607943889143789,0.493449250952099,0.782880258542840,0.652081542238760,-0.404572386885279,1.82706926311258,-0.359354485150894,0.314485894068977,-0.174235056914031,-0.262211070652827,-0.274201975025466,0.0970354456384213,0.679704916451667,0.150610795058247,2.03815215876880,0.531926904094175,-0.00743599693182084,0.659959728639740,0.0549452396647523,0.387940704931653,0.197369174713476,0.840419144074273,0.846741189948881,1.01572555930326,1.05902876858781,0.465884955773762,1.62139528915451,1.47347658156355,2.17456945128983,1.05246344755625,2.08753003264876,0.332706910390444,0.805765416526121,0.879498126439802,2.05565156465806,0.530123673881927,0.504244617511550,0.790686245553730,0.416369874153813,1.08994676938820,0.978958760569113,0.171577967303555,0.376878312569948,1.00028702414490,0.395581358439945,1.16504406715232,-0.196698879780730,-0.354877031758726,0.911061481664713,0.555811880435252,0.134657207388087,1.74833542027882,1.47461897889551,-0.300991849276460,0.785726368436654,0.854090539413190,1.35784879192114,0.685125988227126,0.215289772002443,1.08658137768042,1.07972827850065,0.388891160876784,0.317854092883709,0.876908650502805,1.86234178774766,0.791405414818497,0.983599585718801,1.49165348292694,1.02723732277374,0.619917212457385,0.855495566501820,0.447424704551598,1.56566752393885,0.461274709554902,1.30453710033477,-0.0187370794611071,0.917063137473408,0.528955317748687,1.27256395385792,0.497908078668445,0.618263181922135,1.24656451487193,1.29853316918507,2.02145739133004,0.678830565274483,0.0300703476770143,1.13979030422469,-0.346026053600992,-0.213841634746707,1.72772376655277,0.689450199586051,0.770728122060895,1.19631688660874,-0.475986949197893,0.722671889065240,0.310741545490308,-0.0905184116971585,0.901331152216958,0.277182406518109,1.22319097978148,0.917993566204614,0.336549486480734,0.563823210585823,0.338403397472110,1.26748025058316,1.24914114340607,0.600121292171614,0.000599349877236088,0.288183547303699,1.86736240318596,1.42300574263572,-0.311545347309318,0.799845454543431,1.77672871652827,1.04728136272863,0.192203323015651,1.51665240219762,-0.750484100720113,0.694640471850734,0.0266853365663901,0.681090098223521,1.15077013867516,0.643175706102967,1.45138957291854,2.05893337322988,2.21371363686353,0.463780305989926,0.784800500355734,0.910337448414861,1.57848235791109,-0.0823620745279087,1.07060113389007,0.978982766785247,0.915166939697124,2.51665655661725,0.207841119501063,1.22102485854536,1.14528965505438,0.552270180033793,1.74848827746538,0.162414977957761,0.801456128223735,0.908347589428711,0.500955420006785,0.134143195268509,1.71437561683798,-0.570835452101027,0.355347876803348,1.74391861779308,-0.0254749166174195,-1.04322303189870,1.19155292489141,2.05851664501869,-0.250335919973808,0.553847721426111,0.794270883741650,-0.0181991686551390,0.199017769918627,1.14497225300163,0.383937194600811,0.667372048922489,0.575342096026052,0.766788894217268,2.42754253683587,1.30601788624045};
    vector<float> var,predprobs,h,one_h,r1,r2,assign1,assign2;
    int all = matrix.size()-1;
    matrix[0][0] = 1;
    // trial run
    for (int t = 0; t < n; ++t){
            // MATLAB
            // predprobs = studentpdf(X(t), muT, ...
            //                  betaT.*(kappaT+1)./(alphaT.*kappaT), ...
            //                  2 * alphaT);
            // C++ 
        
        var = cal_var(alphaT,betaT,kappaT);        
        predprobs = studentpdf(x[t],muT,var,all_mul(alphaT,2));
        if (t==300)
            printstudentpdf(x[t],muT,var,all_mul(alphaT,2));

        // cout << endl;
        // print(predprobs);

            // MATLAB
            //   H = hazard_func([1:t]');
            // C++
        vector<float> ht(t+1,0);
        h = hazard(ht,200);
        
            // MATLAB
            // R(2:t+1,t+1) = R(1:t,t) .* predprobs .* (1-H);
            // R(1,t+1) = sum( R(1:t,t) .* predprobs .* H );
            // R(:,t+1) = R(:,t+1) ./ sum(R(:,t+1));
            // C++
        one_h = sub_all(h,1);        
        // r1 = select(matrix,0,t,t);
        assign1 = ele_mul3(select(matrix,0,t,t),predprobs,one_h);
        int j = 0;
        for (int i=1;i<=t+1;++i){
            matrix[i][t+1] = assign1[j];
            j++;
        }
        // matrix_assign(matrix,assign1,1,t+1,t+1);
        
        matrix[0][t+1] = sum_ele_mul3(select(matrix,0,t,t),predprobs,h);
        
        r2 = select(matrix,0,all,t+1);
        assign2 = all_div(r2,sumvec(r2));
        
        // matrix_assign(matrix,assign2,0,all,t+1);
        j = 0;
        for (int i=0;i<=all;++i){
            matrix[i][t+1] = assign2[j];
            j++;
        }
        
        // if (t==67)
        //     // print_matrix(matrix);
        //     printvec(assign2);
            
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
        betaT = cal_beta(x[t],betaT,kappaT,muT);
        muT = cal_mu(x[t],muT,kappaT);
        kappaT = cal_kappa(kappaT);
        alphaT = cal_alpha(alphaT);
            // MATLAB
            // maxes(t) = find(R(:,t)==max(R(:,t)));
            // C++
        maxes[t] = max_index(select(matrix,0,all,t));
    }
    // print_matrix(matrix);
    // printvec(maxes);
    return 0;
}
