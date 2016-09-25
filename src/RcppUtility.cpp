#include <Rcpp.h>
#include "RcppUtility.h"
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]




DataFrame convertList2DataFrame(List L){
    // リストをデータフレームに変換する
    // この方法を使うと新たな要素の値の
    // コピーが発生しないため効率が良い

    // 出力用データにセットする行名・列名を作成します
    R_xlen_t nrows = NumericVector(L[0]).length();
    R_xlen_t ncols = L.length();
    CharacterVector rownames(nrows);
    CharacterVector colnames(ncols);

    // 行名
    rownames = as<CharacterVector>(IntegerVector(seq_len(nrows)));

    // 列名
    if(L.attr("names")==R_NilValue) {
        for(int i=0; i<ncols; ++i){
            colnames[i] = (String("V") += String(i+1));
        }
    } else {
        colnames = L.attr("names");
    }

    // これらの属性に値をセットして返すと
    // R ではデータフレームとして扱われます
    L.attr("row.names") = rownames;
    L.attr("names")     = colnames;
    L.attr("class")     = "data.frame";


    // 型の変換
    // 上の属性をセットしてから行えば
    // 要素の値のコピーが発生しないようです
    return as<DataFrame>(L);

}


//' Convert vector to factor
//' @param x : vector
//' @export
// [[Rcpp::export]]
SEXP fast_factor( SEXP x ) {
    if(Rf_isFactor(x)){ return x;}
    switch( TYPEOF(x) ) {
    case INTSXP:  return fast_factor_template<INTSXP>(x);
    case REALSXP: return fast_factor_template<REALSXP>(x);
    case STRSXP:  return fast_factor_template<STRSXP>(x);
    }
    return R_NilValue;
}



