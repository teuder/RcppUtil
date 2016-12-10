#include <Rcpp.h>
#include "RcppUtility.h"
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

// catdap の算出に使う関数
inline double phi_catdap(double x){
    return x * std::log(x) - 1.0;
}

// yi= 目的
// xj= 説明

template <typename T1, typename T2>
double catdap_base(const T1& y,const T2& x){
    
    // capdapの計算を行う
    // http://jasp.ism.ac.jp/ism/catdap/
    
    // データ数 n
    const R_xlen_t n = y.size();
    if(n != x.size()) {
        stop("The length of vectors x and y must be equal.");
    }
    
    // 入力データをユニークな値だけにする
    const T1 y_levels = unique(y);
    const T2 x_levels = unique(x);
    
    //目的変数と説明変数のカテゴリ数 Cy Cx
    R_xlen_t Cy = y_levels.size();
    R_xlen_t Cx = x_levels.size();
    
    //目的変数と説明変数による２次元分割表 Nyx
    std::vector<std::vector<double>> Nyx (Cy, std::vector<double>(Cx, 0.0));
    
    // ２次元分割表 Nyx に値を入れる
    for(R_xlen_t i=0; i<n; ++i){ //各データ y[i] x[i] について
        // 値の組み合わせを判定する
        for(R_xlen_t yi=0; yi<Cy; ++yi){
            for(R_xlen_t xi=0; xi<Cx; ++xi){
                if(y[i]==y_levels[yi] && x[i]==x_levels[xi]){
                    ++(Nyx[yi][xi]);
                    goto endloop;
                }
            }
        }
        endloop:;
    }
    
    // catdap の算出のために
    // Nyx[i][j]が 0 である場合は 1/e に置き換える
    for(R_xlen_t yi=0; yi<Cy; ++yi){
        for(R_xlen_t xi=0; xi<Cx; ++xi){
            if(Nyx[yi][xi]==0.0){
                Nyx[yi][xi] = 1.0/exp(1.0);
            }
        }
    }
    
    //目的変数と説明変数の各カテゴリ毎のレコード数
    std::vector<double> Ny(Cy, 0.0);
    std::vector<double> Nx(Cx, 0.0);
    // 値をセット
    for(R_xlen_t yi=0; yi<Cy; ++yi){
        for(R_xlen_t xi=0; xi<Cx; ++xi){
            Ny[yi] += Nyx[yi][xi];
            Nx[xi] += Nyx[yi][xi];
        }
    }
    
    // catdapのためのAIC算出
    double PHIyx(0.0),PHIx(0.0),PHIy(0.0),PHIn(0.0);
    
    for(int xi=0; xi<Cx; ++xi){ 
        PHIx += phi_catdap(Nx[xi]);
        for(int yi=0; yi<Cy; ++yi){
            PHIyx += phi_catdap(Nyx[yi][xi]);
            PHIn += Nyx[yi][xi];
        }
    }
    PHIn = phi_catdap(PHIn);
    
    for(int yi=0; yi<Cy; ++yi){
        PHIy += phi_catdap(Ny[yi]);
    }
    
    // xとyが独立ではないと仮定したモデルの AICyx
    // xとyが独立であると仮定したモデルの AICy
    // double AICyx = -2.0 * (PHIyx - PHIx);
    // double AICy  = -2.0 * (PHIy  - PHIn);
    // double AICdiff = AICyx - AICy;
    
    // このAICの差分が負であれば、説明変数 x の値によって
    // 目的変数 y の値の比率は異なっていると判定される
    const double AICdiff = -2.0 * (PHIyx - PHIx - PHIy + PHIn);
    return AICdiff;
}




// [[Rcpp::export]]
double catdap(RObject y, RObject x){
    if(is<IntegerVector>(y)){
        if     (is<IntegerVector>(x))   return catdap_base(IntegerVector(y.get__()),IntegerVector(x.get__()));
        else if(is<NumericVector>(x))   return catdap_base(IntegerVector(y.get__()), NumericVector(x.get__()));
        else if(is<CharacterVector>(x)) return catdap_base(IntegerVector(y.get__()), CharacterVector(x.get__()));
        else {warning("Invalid type x");return 9999.0;}
    } else if(is<NumericVector>(y)){
        if     (is<IntegerVector>(x))   return catdap_base(NumericVector(y.get__()),IntegerVector(x.get__()));
        else if(is<NumericVector>(x))   return catdap_base(NumericVector(y.get__()), NumericVector(x.get__()));
        else if(is<CharacterVector>(x)) return catdap_base(NumericVector(y.get__()), CharacterVector(x.get__()));
        else {warning("Invalid type x");return 9999.0;}
    } else if(is<CharacterVector>(y)){
        if     (is<IntegerVector>(x))   return catdap_base(CharacterVector(y.get__()),IntegerVector(x.get__()));
        else if(is<NumericVector>(x))   return catdap_base(CharacterVector(y.get__()), NumericVector(x.get__()));
        else if(is<CharacterVector>(x)) return catdap_base(CharacterVector(y.get__()), CharacterVector(x.get__()));
        else {warning("Invalid type x");return 9999.0;}
    } else {
        stop("Invalid type y");
    }
    
    return 0.0;
}

//' apply_catdap
//' apply catdap to all columns in a data.frame
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


