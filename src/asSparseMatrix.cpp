//' @useDynLib RcppUtil
//' @importFrom Rcpp evalCpp
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::plugins("cpp11")]]

//' Create sparse matrix from data.frame
//' @param df : data.frame
//' @param checkNA : set false if you can gualantee there is no NAs in the df.
//' @export
// [[Rcpp::export]]
S4 asSparseMatrix( DataFrame df, bool checkNA = true){
  // df : データフレーム、値は数値のみである前提
  // na  df にNAが含まれる場合は 1、 含まれていない場合は 0 にセット

  //データの行数と列数
  int nrow = df.nrows();
  int ncol = df.length();

  // R の Vector では要素の追加（push_back）の効率が悪いので
  // 標準 C++ の vector を使用する。
  std::vector<R_xlen_t>    rows;
  std::vector<R_xlen_t>    cols;
  std::vector<double>      vals;

  //データフレームの全ての要素にアクセスして、
  //0じゃない要素の位置と値を保存する
  if(checkNA){
    //データがNAを含む可能性がある場合はNAのチェックを行う
    for(R_xlen_t col=0; col<ncol; ++col){
      NumericVector column = df[col];
      for(R_xlen_t row=0; row<nrow; ++row){
        // NA の要素は 0 とみなす
        if(NumericVector::is_na(column[row])) {
          stop("NA must be eliminated from the input df. NA exists at row %i column %i.", row+1, col+1);
        } else if(column[row]!= 0.0){
          rows.push_back(row+1);
          cols.push_back(col+1);
          vals.push_back(column[row]);
        }
      }
    }
  } else {
    //データがNAを含まないと保証された場合はNAのチェックを行わない
    for(R_xlen_t col=0; col<ncol; ++col){
      NumericVector column = noNA(NumericVector(df[col]));
      for(R_xlen_t row=0; row<nrow; ++row){
        if(column[row]!= 0.0){
          rows.push_back(row+1);
          cols.push_back(col+1);
          vals.push_back(column[row]);
        }
      }
    }
  }

  // Matrix パッケージから sparseMatrix 関数を呼び出す
  Environment env = Environment::namespace_env("Matrix");
  Function sparseMatrix = env["sparseMatrix"];
  S4 sm = sparseMatrix( Named("i") = wrap(rows),
                        Named("j") = wrap(cols),
                        Named("x") = wrap(vals),
                        Named("dims")=NumericVector::create(nrow,ncol));
  //疎行列に行名と列名をセットする
  //行名には NULL をセットする
  List dimnames = List::create(R_NilValue, df.names());
  sm.attr("Dimnames") = dimnames;

  return sm;
}


S4 asSparseMatrixold( DataFrame df, bool na = false){
    // df : データフレーム、値は数値のみである前提
    // na  df にNAが含まれる場合は 1、 含まれていない場合は 0 にセット

    //データの行数と列数
    int nrow = df.nrows();
    int ncol = df.length();

    // R の Vector では要素の追加（push_back）の効率が悪いので
    // 標準 C++ の list を使用する。
    std::vector<R_xlen_t>    rows;
    std::vector<R_xlen_t>    cols;
    std::vector<double>      vals;

    // rows.reserve(100000);
    // cols.reserve(100000);
    // vals.reserve(100000);

    //データフレームの全ての要素にアクセスして、
    //0じゃない要素の位置と値を保存する
    if(na){//NAありの場合
        //データがNAを含む可能性がある場合は
        //NAのチェックを行う
        for(R_xlen_t col=0; col<ncol; ++col){
            NumericVector column = df[col];
            for(R_xlen_t row=0; row<nrow; ++row){
                // NA の要素は 0 とみなす
                if(NumericVector::is_na(column[row])) continue;
                else if(column[row]!=0.0){
                    rows.push_back(row+1);
                    cols.push_back(col+1);
                    vals.push_back(column[row]);
                }
            }
        }
    } else {
        //データがNAを含まないと保証された場合は
        //NAのチェックを行わない
        for(R_xlen_t col=0; col<ncol; ++col){
            NumericVector column = noNA(NumericVector(df[col]));
            for(R_xlen_t row=0; row<nrow; ++row){
                if(column[row]!=0.0){
                    rows.push_back(row+1);
                    cols.push_back(col+1);
                    vals.push_back(column[row]);
                }
            }
        }
    }

    // Matrix パッケージから sparseMatrix 関数を読み込み
    // Matrix::sparseMatrix() という形式で呼び出したのと同じ
    Environment env = Environment::namespace_env("Matrix");
    Function sparseMatrix = env["sparseMatrix"];
    S4 sm = sparseMatrix( Named("i") = wrap(rows),
                          Named("j") = wrap(cols),
                          Named("x") = wrap(vals),
                          Named("dims")=NumericVector::create(nrow,ncol));
    //疎行列に列名をセットする
    //行名には NULL をセットする
    List dimnames = List::create(R_NilValue, df.names());
    sm.attr("Dimnames") = dimnames;

    return sm;
}


//' @export
// [[Rcpp::export]]
List sparseMatrixToDataFrame(S4 sm){

    // Matrixパッケージの疎行列（dgCMatrix）ではないデータが
    // 渡された場合は拒否します
    if(as<std::string>(sm.attr("class"))!="dgCMatrix")
        stop("Input must be dgCMatrix.");

    // 行数と列数
    //要素数や要素番号の型には R_xlen_t を使います
    IntegerVector dim = sm.slot("Dim");
    R_xlen_t nrow = dim[0];
    R_xlen_t ncol = dim[1];

    // dgCMatrix から値を取り出します
    // dgCMatrix の疎行列の形式は圧縮列格納方式となっています
    NumericVector X = sm.slot("x");//非ゼロ要素の値
    IntegerVector I = sm.slot("i");//非ゼロ要素の行番号（0ベース）
    IntegerVector P = sm.slot("p");//非ゼロ要素の列番号の開始位置（0ベース）

    // P から非ゼロ要素の列番号 J (0ベース) を再生します
    IntegerVector J(X.length());
    IntegerVector repeat = diff(P); //各列番号が連続する数
    R_xlen_t n = repeat.length();
    for(R_len_t j=0; j<n; ++j){
        if(repeat[j]==0) continue; //非ゼロ値が存在しない列は飛ばします
        // J に 各非ゼロ値に対応する列番号 j をセットします
        J[seq(R_xlen_t(P[j]),R_xlen_t(P[j+1]-1))] = rep(j, repeat[j]);
    }

    // 出力用データを作成します
    // ここでは列数が ncol のデータフレームを作成したいですが
    // DataFrame df(ncol)；とすると、１列目の値が ncol である
    // データフレームができてしまいます、。
    // そこで要素数が ncol であるリストを作成してから
    // 最後にデータフレームに変換します
    List L(ncol);

    // リストを初期化するため
    // 要素の値が 0 で長さが nrow のベクトルを作成して
    // そのベクトルをリストの各要素として代入します
    for(R_xlen_t i=0;i<ncol;++i){
        L[i] = NumericVector(nrow);
    }

    // 出力用データに値をセットします
    // 行番号 I、列番号 J、の要素に値 X を代入します
    n = I.length();
    NumericVector column;
    for(R_xlen_t k=0; k<n; ++k){
        // column は L の J[k] 番目の要素（列）への参照となります
        column = L[J[k]];
        // column の I[k] 番目の要素に X[k] の値を代入します
        column[I[k]] = X[k];
    }

    // 出力用データにセットする行名・列名を作成します
    List Dimnames = sm.slot("Dimnames");
    CharacterVector rownames(nrow);
    CharacterVector colnames(ncol);
    // 行名
    if(Dimnames[0]==R_NilValue){
        //行名がNULLの場合は、1 2 ... という行名にします
        IntegerVector v = seq_len(nrow);
        rownames = as<CharacterVector>(v);
    } else {
        rownames=Dimnames[0];
    }
    // 列名
    if(Dimnames[1]==R_NilValue){
        //列名がNULLの場合は、V1 V2 ... という列名にします
        //R の paste() 関数を呼んだら遅かったのでやめました
        for(int i=0; i<ncol; ++i){
            colnames[i] = (String("V") += String(i+1));
        }
    } else {
        colnames=Dimnames[1];
    }

    // これらの属性に値をセットして返すと
    // R ではデータフレームとして扱われます
    L.attr("row.names") = rownames;
    L.attr("names")     = colnames;
    L.attr("class")     = "data.frame";

    return L;
}
