#include <Rcpp.h>
#include "RcppUtility.h"
using namespace Rcpp;

#include <thread>

// [[Rcpp::plugins("cpp11")]]

std::vector<std::vector<R_xlen_t>> make_boundary(const R_xlen_t len, const int n){
    //長さlenのベクターをn個の区画に分割する境界

    R_xlen_t A = len/n; //商
    R_xlen_t B = len%n; //剰余

    //Rcout << "A " << A << "\n";
    //Rcout << "B " << B << "\n";

    std::vector<R_xlen_t> num_elements(n+1,A);
    for(std::size_t i=0; B!=0; --B, ++i){
        num_elements[i+1] += 1;
    }
    num_elements[0]=0;

    //show("boundary", boundary);

    std::vector<std::vector<R_xlen_t>> bound(n);
    for(auto& b : bound) b.resize(2);

    for(R_xlen_t b=0, e=0, i=0; i<n; ++i){
        b += num_elements[i];
        e  = num_elements[i+1] + b;
        bound[i][0]=b;
        bound[i][1]=e;
        //Rcout << "b " << b << " e " << e << "\n";
    }

    return bound;
}



NumericVector myfun(NumericVector v){
    NumericVector out(v.length());
    std::transform(v.begin(),v.end(), out.begin(), ::sqrt);
    return out;
}

std::vector<double> myfun2(std::vector<double>& v){
    std::vector<double> out(v.size());
    std::transform(v.begin(),v.end(), out.begin(), ::sqrt);
    return out;
}

NumericVector myfun3(NumericVector v){
    NumericVector out(v.size());
    std::transform(v.begin(),v.end(), out.begin(), ::sqrt);
    return out;
}

//' apply funciton to data.frame column
//'
//' @param df : data.frame
//' @param n : number of cores for parallel operation
//' @export
// [[Rcpp::export]]
DataFrame thread_function(DataFrame df, int n) {

    R_xlen_t len = df.length();

    std::vector<std::vector<R_xlen_t>> boundary = make_boundary(len,n);

    List out(len);

    //個々のスレッドで実行する処理
    auto worker = [&](R_xlen_t bgn, R_xlen_t end){
        Rcout << "thread id: " << std::this_thread::get_id() << "\n";
        // １つのスレッドでは列番号 b～e-1までの列を処理する
        for(R_xlen_t i = bgn; i < end; ++i){
            SEXP p = df[i];
            out[i] = myfun(p);
        }
    };

    //スレッドの作成と実行
    std::vector<std::thread> threads;
    for(R_xlen_t i=0; i<n; ++i){
        Rcout << "boundary[i][0] " << boundary[i][0] << " boundary[i][1] " << boundary[i][1] << "\n";
        threads.push_back(std::thread(worker, boundary[i][0], boundary[i][1]));
    }

    //全てのスレッドの終了を待つ
    for (std::thread &th : threads) {
        th.join();
    }

    return convertList2DataFrame(out);
}


//' apply funciton to data.frame column
//'
//' @param df : data.frame
//' @param n : number of cores for parallel operation
//' @export
// [[Rcpp::export]]
DataFrame thread_function2(std::vector<std::vector<double>> df, int n) {

    R_xlen_t len = df.size();

    std::vector<std::vector<R_xlen_t>> boundary = make_boundary(len,n);

    std::vector<std::vector<double>> out(len);

    //個々のスレッドで実行する処理
    auto worker = [&](R_xlen_t bgn, R_xlen_t end){
        //Rcout << "thread id: " << std::this_thread::get_id() << "\n";
        // １つのスレッドでは列番号 b～e-1までの列を処理する
        for(R_xlen_t i = bgn; i < end; ++i){
            //SEXP p = df[i];
            out[i] = myfun2(df[i]);
        }
    };

    //スレッドの作成と実行
    std::vector<std::thread> threads;
    for(R_xlen_t i=0; i<n; ++i){
        //Rcout << "boundary[i][0] " << boundary[i][0] << " boundary[i][1] " << boundary[i][1] << "\n";
        threads.push_back(std::thread(worker, boundary[i][0], boundary[i][1]));
    }

    //全てのスレッドの終了を待つ
    for (std::thread &th : threads) {
        th.join();
    }

    //return convertList2DataFrame(out);
    return wrap(out);
}

//' apply funciton to data.frame column
//'
//' @param df : data.frame
//' @param n : number of cores for parallel operation
//' @export
// [[Rcpp::export]]
DataFrame thread_function3(DataFrame df, int n) {

    R_xlen_t len = df.size();

    std::vector<std::vector<R_xlen_t>> boundary = make_boundary(len,n);

    List out(len);

    //個々のスレッドで実行する処理
    auto worker = [&](R_xlen_t bgn, R_xlen_t end){
        //Rcout << "thread id: " << std::this_thread::get_id() << "\n";
        // １つのスレッドでは列番号 b～e-1までの列を処理する
        for(R_xlen_t i = bgn; i < end; ++i){
            SEXP p = df[i];
            out[i] = myfun3(p);
        }
    };

    //スレッドの作成と実行
    std::vector<std::thread> threads;
    for(R_xlen_t i=0; i<n; ++i){
        //Rcout << "boundary[i][0] " << boundary[i][0] << " boundary[i][1] " << boundary[i][1] << "\n";
        threads.push_back(std::thread(worker, boundary[i][0], boundary[i][1]));
    }

    //全てのスレッドの終了を待つ
    for (std::thread &th : threads) {
        th.join();
    }

    return convertList2DataFrame(out);
    //return (out);
}







// [[Rcpp::export]]
List makeDummyList(SEXP x,
                     bool fullrank = true,
                     const std::string prefix = "")
{
    //１つの factor をダミー変数化したリストを返す

    IntegerVector v = fast_factor(x);

    R_xlen_t n = v.length();
    CharacterVector levels = v.attr("levels");
    int L = levels.length();

    List out(L);

    for(int i = 0; i<L; ++i){
        out[i] = NumericVector(n);
    }

    for(R_xlen_t i = 0; i<n; ++i){
        for(int j = 0; j <= L; ++j){
            if(v[i]==j) {
                NumericVector x = out[j-1];
                x[i] = 1.0;
                break;
            }
        }
    }
    //カテゴリ値を要素名に設定する
    CharacterVector names(L);
    for(R_xlen_t i=0; i<L; ++i){
        names[i] = String(prefix + levels[i]);
    }
    out.names()=names;


    if(fullrank){
        out.erase(0);
    }

    return out;
}


// [[Rcpp::export]]
DataFrame makeDummyDf(SEXP x,
                   bool fullrank = true,
                   const std::string prefix = ""){
    return convertList2DataFrame(makeDummyList(x, fullrank, prefix));
}




class MakeDummy {
public:

    MakeDummy(const DataFrame& df_, const List& LL_, const bool& fullrank_, const std::string sep_cat_)
    //:df(df_),LL(LL_),fullrank(fullrank_),sep_cat(sep_cat_)
    {
        df = df_;
        LL=LL_;
        fullrank=fullrank_;
        sep_cat=sep_cat_;

    };

    // 関数オブジェクト
    void operator()(const R_xlen_t bgn, const  R_xlen_t end){
        //std::cout << "thread id: " << std::this_thread::get_id() << std::endl;
        for(R_xlen_t i = bgn; i < end; ++i){
            //std::cout << "thread id: " << std::this_thread::get_id() << " col " << i << std::endl;
            std::cout  << i << std::endl;
            IntegerVector v = *(df.begin()+i);

            R_xlen_t n = v.length();
            CharacterVector levels = v.attr("levels");
            int L = levels.length();

            List out(L);

            for(int i = 0; i<L; ++i){
                *(out.begin()+i)= NumericVector(n);
                //out[i] = NumericVector(n);
            }

            for(R_xlen_t i = 0; i<n; ++i){
                for(int j = 0; j <= L; ++j){
                    if(v[i]==j) {
                        NumericVector x = out[j-1];
                        x[i] = 1.0;
                        break;
                    }
                }
            }
            //カテゴリ値を要素名に設定する
            CharacterVector names(L);
            for(R_xlen_t i=0; i<L; ++i){
                names[i] = String(sep_cat + levels[i]);
            }
            out.names()=names;


            if(fullrank){
                out.erase(0);
            }
            *(LL.begin()+i)=out;
            //LL[i] = out;
        }
    };



private:
    DataFrame df;
    List LL;
    bool fullrank;
    std::string sep_cat;
};













List convertLL2LV(List LL){

    //変数名
    CharacterVector varnames = LL.names();

    //Lの各要素の長さをカウントする
    R_xlen_t ncols = 0;
    for(List l : LL){
        ncols += l.length();
    }

    //LLの各要素を１つのリストにまとめて
    //ベクトルを要素とするリスト LV を作成する
    List LV(ncols);//値
    CharacterVector names(ncols); //要素名
    //show_v("names",names);

    R_xlen_t k=0;
    for(R_xlen_t i=0, Nll=LL.length(); i<Nll; ++i){
        List lli = LL[i];//
        CharacterVector catnames = lli.names();//変数 i のカテゴリ値
        for(R_xlen_t j=0, Nlli = lli.length(); j<Nlli; ++j){

            //要素に値をセット
            LV[k] = lli[j];

            //変数名とカテゴリ値を合体して要素名を作成
            String var = varnames[i];
            var += catnames[j];
            names[k] = var;

            ++k;

        }
    }
    LV.names() = names;
    return LV;
}





// [[Rcpp::export]]
DataFrame make_dummy_p( DataFrame df
                            ,const int cores = 2
                            ,const std::string sep_cat = "_"
                            ,const bool fullrank = true
) {


    CharacterVector names=df.attr("names");
    //LL.names()=df.names();
    //show_v("names",names);

    //データフレームのカラム数
    R_xlen_t len = df.length();

    //１つのカテゴリ変数をダミー変数に変換した結果（）を
    //１つの要素
    List LL(len);


    //長さlenのベクターを cores 個の区画に分割する境界を作成
    //std::vector<R_xlen_t> boundary = make_boundary(len, cores);
    std::vector<std::vector<R_xlen_t>> boundary = make_boundary(len, cores);


    //個々のスレッドで実行する処理をラムダ式で記述
    //makeDummyList()はリストを返すので
    //outはリストを要素とするリストとなる
    auto f = [&](const R_xlen_t bgn, const  R_xlen_t end, const bool fullrank=true){
        std::cout << "thread id: " << std::this_thread::get_id() << std::endl;
        for(R_xlen_t i = bgn; i < end; ++i){
            SEXP p = df[i];
            LL[i] = makeDummyList(p, fullrank, sep_cat);
        }
    };

    //スレッドの作成と実行
    std::vector<std::thread> threads;
    threads.reserve(cores);
    for(R_xlen_t i=0; i<cores; ++i){
        //Rcout << "b " << boundary[i][0] << " e " << boundary[i][1] << "\n";

        if((i==0) & (!fullrank)){
            threads.push_back(std::thread(f, boundary[i][0], boundary[i][1], false));
        }else{
            threads.push_back(std::thread(f, boundary[i][0], boundary[i][1], true));
        }

    }

    // //全てのスレッドの終了を待つ
    for (std::thread &th : threads) {
        th.join();
    }

    // リストを要素とするリストを
    // ベクターを要素とするリストに変換する
    LL.names() = df.names();
    List LV = convertLL2LV(LL);
    CharacterVector name = LV.names();

    // LV をデータフレームに変換して返す
    return convertList2DataFrame(LV);
    //return df;
}


// [[Rcpp::export]]
DataFrame make_dummy_pp( DataFrame df
                             ,const int cores = 2
                             ,const std::string sep_cat = "_"
                             ,const bool fullrank = true
) {


    CharacterVector names=df.attr("names");
    //LL.names()=df.names();
    //show_v("names",names);

    //データフレームのカラム数
    R_xlen_t len = df.length();

    //１つのカテゴリ変数をダミー変数に変換した結果（）を
    //１つの要素
    List LL(len);

    //長さlenのベクターを cores 個の区画に分割する境界を作成
    std::vector<std::vector<R_xlen_t> > boundary = make_boundary(len, cores);
    //return wrap(boundary);

    //スレッドの作成と実行
    std::vector<std::thread> threads;
    threads.reserve(cores);
    for(R_xlen_t i=0; i<cores; ++i){
        //Rcout << "b " << boundary[i][0] << " e " << boundary[i][1] << "\n";
        if((i==0) & (!fullrank)){
            threads.push_back(std::thread(MakeDummy(df, LL, false, sep_cat), boundary[i][0], boundary[i][1]));
        }else{
            threads.push_back(std::thread(MakeDummy(df, LL,  true, sep_cat), boundary[i][0], boundary[i][1]));
        }
    }

    // //全てのスレッドの終了を待つ
    for (std::thread &th : threads) {
        th.join();
    }


    // リストを要素とするリストを
    // ベクターを要素とするリストに変換する
    LL.names() = df.names();
    List LV = convertLL2LV(LL);
    CharacterVector name = LV.names();


    // LV をデータフレームに変換して返す
    return convertList2DataFrame(LV);
    //return df;
}






