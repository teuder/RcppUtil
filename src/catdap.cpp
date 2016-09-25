#include <Rcpp.h>
#include "RcppUtility.h"
using namespace Rcpp;


// [[Rcpp::plugins("cpp11")]]


template <typename T1, typename T2>
double catdap_base(const T1& y,const T2& x){
    using namespace std;

    //capdapの計算を行う
    // http://jasp.ism.ac.jp/ism/catdap/

    //データ数 n
    const R_xlen_t n = y.length();
    if(n!=x.length()) {
        stop("The length of vectors x and y must be equal.");
    }
    //Rcout << n << "\n";

    const T1 y_levels = unique(y);
    const T2 x_levels = unique(x);

    // show(y_levels);
    // show(x_levels);

    //目的変数と説明変数のカテゴリ数 ny nx
    const unsigned int ny = y_levels.size();
    const unsigned int nx = x_levels.size();

    //目的変数と説明変数による２次元分割表 Nyx
    vector<vector<R_xlen_t>> Nyx (ny, vector<R_xlen_t>(nx, 0));

    //目的変数と説明変数の各カテゴリを巡回するためのインデックス
    const IntegerVector yy = seq_len(ny)-1;
    const IntegerVector xx = seq_len(nx)-1;

    //目的変数と説明変数の各カテゴリのレコード数
    std::vector<R_xlen_t> Ny(ny);
    std::vector<R_xlen_t> Nx(nx);


    // x と y の２次元度数分布表 Nyx に値を入れる
    for(R_xlen_t i=0; i<n; ++i){
        for(int xi : xx){
            for(int yi : yy){
                if(y[i]==y_levels[yi] && x[i]==x_levels[xi]){
                    ++(Nyx[yi][xi]);
                    goto endloop;
                }
            }
        }
        endloop:;
    }

    // x と y の各カテゴリ値の度数分布に値を入れる
    for(int yi : yy){
        for(int xi : xx){
            Ny[yi] += Nyx[yi][xi];
            Nx[xi] += Nyx[yi][xi];
        }
    }
    // show(Ny);
    // show(Nx);


    //xとyが独立であると仮定したモデルのAIC_0
    //xとyが独立ではないと仮定したときのAIC_1
    // diffAIC = AIC_1 - AIC_0
    double diffAIC = 0.0;
    for(int xi : xx){
        for(int yi : yy){
            //1.0をかけているのは結果を実数に変換するため
            diffAIC += Nyx[yi][xi] * log(1.0*n*Nyx[yi][xi]/(Nx[xi]*Ny[yi]));
        }
    }

    return 2.0*((ny-1)*(nx-1) - diffAIC);
}



// [[Rcpp::export]]
double catdap(RObject y, RObject x){
    if(is<IntegerVector>(y)){
        if(is<IntegerVector>(x)) return catdap_base(as<IntegerVector>(y),as<IntegerVector>(x));
        else if(is<NumericVector>(x)) return catdap_base(as<IntegerVector>(y),as<NumericVector>(x));
        else if(is<CharacterVector>(x)) return catdap_base(as<IntegerVector>(y),as<CharacterVector>(x));
        else stop("Invalid type x");
    } else if(is<NumericVector>(y)){
        if(is<IntegerVector>(x)) return catdap_base(as<NumericVector>(y),as<IntegerVector>(x));
        else if(is<NumericVector>(x)) return catdap_base(as<NumericVector>(y),as<NumericVector>(x));
        else if(is<CharacterVector>(x)) return catdap_base(as<NumericVector>(y),as<CharacterVector>(x));
        else stop("Invalid type x");
    } else if(is<CharacterVector>(y)){
        if(is<IntegerVector>(x)) return catdap_base(as<CharacterVector>(y),as<IntegerVector>(x));
        else if(is<NumericVector>(x)) return catdap_base(as<CharacterVector>(y),as<NumericVector>(x));
        else if(is<CharacterVector>(x)) return catdap_base(as<CharacterVector>(y),as<CharacterVector>(x));
        else stop("Invalid type x");
    } else {
        stop("Invalid type y");
    }

    return 0.0;
}

//' apply_catdap
//' @param y : vector
//' @param df : data.frame
//' @export
// [[Rcpp::export]]
NumericVector apply_catdap(IntegerVector y, DataFrame df){
    R_xlen_t n = df.size();
    NumericVector out(n);
    for(R_xlen_t i=0; i<n; ++i){
        out[i] = catdap(y, df[i]);
    }
    return out;
}


