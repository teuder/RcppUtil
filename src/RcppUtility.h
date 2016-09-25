#ifndef RcppUtility_h
#define RcppUtility_h

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

template <typename V>
void show(std::string s, V x){
    Rcout << s << " ";
    for(auto xx : x) Rcout << xx << " ";
    Rcout << "\n";
}


template <typename V>
void show_v(std::string s, V x){
    Rcout << s << " ";
    for(auto xx : x) Rcout << xx << " ";
    Rcout << "\n";

}


template <typename S>
void show_s(std::string s, S v){
    Rcout << s << " " << v << "\n";
}

template <int RTYPE>
IntegerVector fast_factor_template( const Vector<RTYPE>& x ) {
    Vector<RTYPE> levs = sort_unique(x);
    IntegerVector out = match(x, levs);
    out.attr("levels") = as<CharacterVector>(levs);
    out.attr("class") = "factor";
    return out;
}


DataFrame convertList2DataFrame(List L);
SEXP fast_factor( SEXP x );


#endif
