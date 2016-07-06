//
//  main.cpp
//  E2E
//
//  Created by 王青龙 on 6/30/16.
//  Copyright © 2016 王青龙. All rights reserved.
//

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <math.h>

using namespace std;

#define PROB_SMOOTH 1e-7

void read_index_word(string str, map<string, int>& word2index){
    ifstream fin(str);
    string line;
    string word;
    int index;

    while (getline(fin, line)) {
        istringstream linestream(line);
        linestream>>index;
        linestream>>word;
        word2index.insert({word,index});
    }
}

void read_corpus(string str, string str2, map<string, int>& word2index){
    ifstream fin(str);
    ofstream fout(str2);
    string line;
    string word;
    int index;
    bool firstword=true;

    while (getline(fin, line)) {
        firstword=true;
        istringstream linestream(line);
        while (linestream>>word) {
            if (firstword) {
                firstword=false;
            }else{
                fout<<" ";
            }
            if (word2index.find(word)==word2index.end()) {
                fout<<-1;
            }else{
                index=word2index[word];
                fout<<index;
            }
        }
        fout<<endl;
    }
    fout.close();
}


void cre_sur(double** (&sur), const int& max_line_en1, const int& max_line_en2){
    const int len1=max_line_en1;
    const int len2=max_line_en2;
    sur=new double* [len2];
    for (int i=0; i!=len2; ++i) {
        sur[i]=new double[len1];
    }
    for (int i=0; i!=len2; ++i) {
        for (int j=0; j!=len1; ++j) {
            sur[i][j]=0.0;
        }
    }

}

void compute_sur(double** (&sur), vector<string>& v_en1_original, vector<string>& v_en2_original, double rho){

    const int len1=(int)v_en1_original.size();
    const int len2=(int)v_en2_original.size();

    for (int i=0; i!=len2; ++i) {
        for (int j=0; j!=len1; ++j) {
            int len_word_min=(int)min(v_en2_original[i].size(), v_en1_original[j].size());
            int k=0;
            while (k<len_word_min) {
                if (v_en2_original[i][k]==v_en1_original[j][k]) {
                    ++k;
                }else{
                    break;
                }
            }
            double p=(double)k/max(v_en2_original[i].size(), v_en1_original[j].size());
            sur[i][j]=pow(M_E, rho*(p-1));
        }
    }
}

void cre_sem(double** (&sem), const int& max_line_en1, const int& max_line_en2){
    const int len1=max_line_en1;
    const int len2=max_line_en2;
    sem=new double* [len2];
    for (int i=0; i!=len2; ++i) {
        sem[i]=new double[len1];
    }
    for (int i=0; i!=len2; ++i) {
        for (int j=0; j!=len1; ++j) {
            sem[i][j]=0.0;
        }
    }
}

double lookup(map<int, map<int,double>>& dict, int first, int second)
{
    if(dict.find(first)!=dict.end())
        if(dict[first].find(second)!=dict[first].end())
            return dict[first][second];
    return 0;
}


void compute_sem(double** (&sem), vector<int>& v_ch, vector<int>& v_en1, vector<int>& v_en2, map<int, map<int, double>>& t_c2e, map<int, map<int, double>>& t_e2c){
    const int len1=(int)v_en1.size();
    const int len2=(int)v_en2.size();
//    const int lenc=(int)v_ch.size();
    for (int i=0; i!=len2; ++i) {

        for (int j=0; j!=len1; ++j) {
//            sem[v_en2[i]][v_en1[j]]=0.0;
            sem[i][j]=0.0;
//            for (int k=0; k!=lenc; ++k) {
//                double e2c = lookup(t_e2c, v_en2[i], v_ch[k]);
//                double c2e = lookup(t_c2e, v_ch[k], v_en1[j]);
//                sem[i][j]+= e2c * c2e;
//            }

            for (auto iter=t_e2c[v_en2[i]].begin(); iter!=t_e2c[v_en2[i]].end(); ++iter) {
                double e2c = iter->second;
                double c2e = lookup(t_c2e, iter->first, v_en1[j]);
                sem[i][j] += e2c * c2e;
            }

        }
    }
}



/*
void compute_sem(double** (&sem), vector<int>& v_en1, vector<int>& v_en2, map<int, map<int, double>>& t_c2e, map<int, map<int, double>>& t_e2c){
    const int len1=(int)v_en1.size();
    const int len2=(int)v_en2.size();
    const int lenc=(int)v_ch.size();
    for (int i=0; i!=len2; ++i) {
        for (int j=0; j!=len1; ++j) {
            //            sem[v_en2[i]][v_en1[j]]=0.0;
            sem[i][j]=0.0;
            for (int k=0; k!=lenc; ++k) {
                if (t_e2c.find(v_en2[i])!=t_e2c.end()&&t_e2c[v_en2[i]].find(v_ch[k])!=t_e2c[v_en2[i]].end()&&t_c2e.find(v_ch[k])!=t_c2e.end()&&t_c2e[v_ch[k]].find(v_en1[j])!=t_c2e[v_ch[k]].end()) {
                    sem[i][j]+=t_e2c[v_en2[i]][v_ch[k]]*t_c2e[v_ch[k]][v_en1[j]];
                }else{
                    if (t_e2c.find(v_en2[i])!=t_e2c.end()&&t_e2c[v_en2[i]].find(v_ch[k])!=t_e2c[v_en2[i]].end()) {
                        sem[i][j]+=t_e2c[v_en2[i]][v_ch[k]]*PROB_SMOOTH;
                    }else{
                        if (t_c2e.find(v_ch[k])!=t_c2e.end()&&t_c2e[v_ch[k]].find(v_en1[j])!=t_c2e[v_ch[k]].end()) {
                            sem[i][j]+=PROB_SMOOTH*t_c2e[v_ch[k]][v_en1[j]];
                        }else{
                            sem[i][j]+=PROB_SMOOTH*PROB_SMOOTH;
                        }
                    }
                }

            }
            
        }
    }
}*/

void cre_simi(double** (&simi), const int& max_line_en1, const int& max_line_en2){
    const int len1=max_line_en1;
    const int len2=max_line_en2;
    simi=new double* [len2];
    for (int i=0; i!=len2; ++i) {
        simi[i]=new double[len1];
    }
    for (int i=0; i!=len2; ++i) {
        for (int j=0; j!=len1; ++j) {
            simi[i][j]=0.0;
        }
    }
}

void compute_simi(double** (&simi), double** (&sem), double** (&sur), double alpha, const vector<int>& v_en1, const vector<int>& v_en2){
    const int len1=(int)v_en1.size();
    const int len2=(int)v_en2.size();

    for (int i=0; i!=len2; ++i) {
        for (int j=0; j!=len1; ++j) {
            simi[i][j]=alpha*sem[i][j]+(1-alpha)*sur[i][j];
        }
    }
}

void cre_net_n(double** (&net_n), const int& max_line_en1, const int& max_line_en2){//翻译概率
    const int I=2*max_line_en2;
    const int J=max_line_en1;
    net_n=new double* [I];
    for (int i=0; i!=I; ++i) {
        net_n[i]=new double[J];
    }
    for (int i=0; i!=I; ++i) {
        for (int j=0; j!=J; ++j) {
            net_n[i][j]=0.0;
        }
    }
}

void ini_net_n(double** (&net_n), double** (&simi), vector<int>& v_en1, vector<int>& v_en2){
    const int I=2*(int)v_en2.size();
    const int J=(int)v_en1.size();
    double sum=0.0;
    for (int j=0; j!=J; ++j) {
        sum=0;
        for (int i=0; i<I/2; ++i) {
            net_n[i][j]=simi[i][j];
            sum+=net_n[i][j];

        }
        for (int i=I/2; i<I; ++i) {
            net_n[i][j]=PROB_SMOOTH;//暂定
            sum+=net_n[i][j];
        }
        if (sum) {
            for (int i=0; i<I; ++i) {
                net_n[i][j]/=sum;
            }
        }
    }
}

void cre_net_e(double** (&net_e), const int& max_line_en2){
    const int I=2*max_line_en2;
    net_e=new double* [I];
    for (int i=0; i!=I; ++i) {
        net_e[i]=new double[I];
    }
    for (int i=0; i!=I; ++i) {
        for (int j=0; j!=I; ++j) {
            net_e[i][j]=0.0;
        }
    }
}

double compute_c(int d, int k){
    double res;
    res=pow(1+abs(abs(d)-1), -k);
    return res;
}

void init_net_e(double** (&net_e), vector<int>& v_en2, int k, double p0){
    const int I=2*(int)v_en2.size();
    for (int i=0; i!=I; ++i) {
        for (int j=0; j!=I; ++j) {
            net_e[i][j]=p0;
        }
    }
    for (int i=0; i!=I/2; ++i) {
        double sum=0.0;
        for (int j=0; j!=I/2; ++j) {
            sum+=compute_c(i-j,k);
        }
        for (int j=0; j!=I/2; ++j) {
            net_e[i][j]=(double)compute_c(i-j,k)/sum;
        }
    }
}

void read_t(string str, map<int, map<int, double>>& t){
    ifstream fin(str);
    string line;
    int index_ch;
    int index_en;
    double prop;
    while (getline(fin, line)) {
        istringstream linestream(line);
        linestream>>index_ch;
        if (t.find(index_ch)==t.end()) {
            t[index_ch]=map<int, double>();
        }
        linestream>>index_en;
        linestream>>prop;
        if(prop<0.001)
            continue;
        t[index_ch].insert({index_en, prop});
    }
}

void read_index_ch(string str, vector<vector<int>>& vv){
    ifstream fin(str);
    string line;
    int index;
    int i=0;

    while (getline(fin, line)) {
        istringstream linestream(line);
        vv.push_back(vector<int>());
        vv[i].push_back(0);//在每一行开头增加一个空位置
        while (linestream>>index) {
            vv[i].push_back(index);
        }
        ++i;
    }
}

void read_index_en(string str, vector<vector<int>>& vv){
    ifstream fin(str);
    string line;
    int index;
    int i=0;

    while (getline(fin, line)) {
        istringstream linestream(line);
        vv.push_back(vector<int>());
        while (linestream>>index) {
            vv[i].push_back(index);
        }
        ++i;
    }
}

void read_original_en(string str, vector<vector<string>>& vv){
    ifstream fin(str);
    string line;
    string word;
    int i=0;

    while (getline(fin, line)) {
        istringstream linestream(line);
        vv.push_back(vector<string>());
        while (linestream>>word) {
            vv[i].push_back(word);
        }
        ++i;
    }
}

void HMMRealViterbi(ofstream& fout, vector<int>& viterbi_alignment, vector<int>& vitar, vector<int>& v_en1, vector<int>& v_en2, double** (&net_n), double** (&net_e)){
    const int I=2*(int)v_en2.size();
    const int J=(int)v_en1.size();
    const int N=I*J;
    vector<double> alpha(N, -1);
    vector<double*> bp(N, (double*)0);
    vitar.resize(J);

    vector<double> alphainit(I, 1.0);//这里只是占坑，真正的初始化在下面
    //这里的初始化情况仅仅对于第一次有效
    double sum_alphainit=0.0;
    for (int i=0; i<I; ++i) {
        alphainit[i]=(i<I/2)?1:(2.0/I);//第一遍初始化分为前半段和后半段
        sum_alphainit+=alphainit[i];
    }
    for (int i=0; i<I; ++i) {
        alphainit[i]/=sum_alphainit;
    }

    for (int i=0; i<I; ++i) {
        alpha[i]=alphainit[i]*net_n[i][0];
        //        cout << I << ":" << J << "   " << i << "    beta    " << beta[i] <<"   betainit    "<<betainit[i]<<"    net_n    "<<net_n[i][0]<< endl;
        if (i>I/2) {//总感觉这里有问题，应该是>=
            alpha[i]=0;
        }
        bp[i]=0;
    }

    auto cur_alpha=alpha.begin()+I;
    //    auto cur_bp=bp.begin()+I;
    double** cur_bp=(&*bp.begin())+I;
    for (int j=1; j<J; ++j) {
        for (int ti=0; ti<I; ++ti, ++cur_alpha, ++cur_bp) {
            double* prev_alpha=&*(alpha.begin()+I*(j-1));

            //            auto prev_alpha=alpha.begin()+I*(j-1);
            double this_node=net_n[ti][j];//翻译概率
            for (int pi=0; pi<I; ++pi,++prev_alpha) {
                //                cout<<*prev_alpha<<endl;
                //                double test=*prev_alpha;
                double alpha_increment=*prev_alpha*net_e[pi][ti]*this_node;
                if (alpha_increment>*cur_alpha) {
                    (*cur_alpha)=alpha_increment;
                    (*cur_bp)=prev_alpha;//存放是指针
                }
            }
        }
    }

    vector<double> betainit(I, 1.0);
    for (int i=0; i<I; ++i) {
        alpha[N-I+i]*=betainit[i];
    }
    int j=J-1;
    cur_alpha=alpha.begin()+j*I;
    vitar[J-1]=int(max_element(cur_alpha, cur_alpha+I)-cur_alpha);
    while (bp[vitar[j]+j*I]) {//bp里面放的是此点的前一个点的指针
        cur_alpha-=I;
        vitar[j-1]=int(bp[vitar[j]+j*I]-(&*cur_alpha));//这样减掉本列的初始值确实得到本列的偏移值
        --j;
    }

    viterbi_alignment.resize(J);
    for (int j=1; j<=J; ++j) {
        viterbi_alignment[j]=vitar[j-1]+1;
        if (viterbi_alignment[j]>I) {
            viterbi_alignment[j]=0;
        }
    }

    bool firstword=true;
    for (int j=1; j<=J; ++j) {
        if (viterbi_alignment[j]) {
            if (firstword) {
                firstword=false;
            }else{
                fout<<" ";
            }
//            fout<<viterbi_alignment[j]-1<<"-"<<j-1;
            fout<<j-1<<"-"<<viterbi_alignment[j]-1;//保证前后顺序一样
        }
    }
    fout<<endl;

}

string dir_name="/Users/wangql/Desktop/E2E/";
map<string, int> word2index_ch, word2index_en1;
map<int, map<int, double>> t_c2e, t_e2c;

double** net_n;
double** net_e;
double** sur;
double** sem;
double** simi;

void loadVocab()
{
    //把index2Word转储成Word2index
    read_index_word(dir_name+"corpus.ch.vcb", word2index_ch);
    read_index_word(dir_name+"corpus.en.vcb", word2index_en1);
    //read_index_word(dir_name+"corpus.en.vcb", word2index_en2);
}

void loadTable()
{

    read_t(dir_name+"s2t64.t3.final", t_c2e);
    read_t(dir_name+"t2s64.t3.final", t_e2c);
}

void align(string file1, string file2, string output)
{
    //把Word转储成index形式
    read_corpus(dir_name+"src.lctok", dir_name+"src.lctok.index", word2index_ch);
    read_corpus(dir_name+file1+".lctok", dir_name+file1+".lctok.index", word2index_en1);
    read_corpus(dir_name+file2+".lctok", dir_name+file2+".lctok.index", word2index_en1);

    vector<vector<int>> vv_ch, vv_en1, vv_en2;
    read_index_ch(dir_name+"src.lctok.index", vv_ch);
    read_index_en(dir_name+file1+".lctok.index", vv_en1);
    read_index_en(dir_name+file2+".lctok.index", vv_en2);

    vector<vector<string>> vv_en1_original, vv_en2_original;
    read_original_en(dir_name+file1+".lctok", vv_en1_original);
    read_original_en(dir_name+file2+".lctok", vv_en2_original);


    vector<int> vit;
    vector<int> viterbi_alignment;

    ofstream fout(dir_name+output);
    for (int i=0; i<vv_en1.size(); ++i) {

        compute_sem(sem, vv_ch[i], vv_en1[i], vv_en2[i], t_c2e, t_e2c);
        compute_sur(sur, vv_en1_original[i], vv_en2_original[i], 3.0);
        compute_simi(simi, sem, sur, 0.3, vv_en1[i], vv_en2[i]);
        ini_net_n(net_n, simi, vv_en1[i], vv_en2[i]);
        init_net_e(net_e, vv_en2[i], 2, PROB_SMOOTH);
        HMMRealViterbi(fout, viterbi_alignment, vit, vv_en1[i], vv_en2[i], net_n, net_e);
    }
    fout.close();

}

int main(int argc, const char * argv[]) {

    loadVocab();
    loadTable();

    cre_net_n(net_n, 300, 300);

    cre_net_e(net_e, 300);

    cre_sur(sur, 300, 300);

    cre_sem(sem, 300, 300);

    cre_simi(simi, 300, 300);

    string file1 = "ref2";
    string file2 = "ref1";
    //align(file1, file2, file1+"-"+file2+".align");

    file1 = "ref4";
    file2 = "ref1";
    align(file1, file2, file1+"-"+file2+".align");

//    file1 = "ref1";
//    file2 = "ref3";
//    align(file1, file2, file1+"-"+file2+".align");
//
//    file1 = "ref1";
//    file2 = "ref4";
//    align(file1, file2, file1+"-"+file2+".align");
//
//    file1 = "ref2";
//    file2 = "ref1";
//    align(file1, file2, file1+"-"+file2+".align");
//
//    file1 = "ref3";
//    file2 = "ref1";
//    align(file1, file2, file1+"-"+file2+".align");
//
//    file1 = "ref4";
//    file2 = "ref1";
//    align(file1, file2, file1+"-"+file2+".align");



}
