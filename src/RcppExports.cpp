// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// ZI_GenerateAlphaMNMetro
NumericVector ZI_GenerateAlphaMNMetro(NumericVector Par, NumericVector Y, NumericVector Z, NumericMatrix Covar, NumericVector propsigma, NumericVector priormu, NumericVector priorSigmas);
RcppExport SEXP _ZIPBayes_ZI_GenerateAlphaMNMetro(SEXP ParSEXP, SEXP YSEXP, SEXP ZSEXP, SEXP CovarSEXP, SEXP propsigmaSEXP, SEXP priormuSEXP, SEXP priorSigmasSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Par(ParSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Covar(CovarSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type propsigma(propsigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type priormu(priormuSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type priorSigmas(priorSigmasSEXP);
    rcpp_result_gen = Rcpp::wrap(ZI_GenerateAlphaMNMetro(Par, Y, Z, Covar, propsigma, priormu, priorSigmas));
    return rcpp_result_gen;
END_RCPP
}
// ZI_GenerateBetaMetro
NumericVector ZI_GenerateBetaMetro(NumericVector Par, NumericVector U, NumericMatrix Covar, NumericVector propsigma, NumericVector priorgamma);
RcppExport SEXP _ZIPBayes_ZI_GenerateBetaMetro(SEXP ParSEXP, SEXP USEXP, SEXP CovarSEXP, SEXP propsigmaSEXP, SEXP priorgammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Par(ParSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type U(USEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Covar(CovarSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type propsigma(propsigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type priorgamma(priorgammaSEXP);
    rcpp_result_gen = Rcpp::wrap(ZI_GenerateBetaMetro(Par, U, Covar, propsigma, priorgamma));
    return rcpp_result_gen;
END_RCPP
}
// ZI_GenerateBigJoint
NumericVector ZI_GenerateBigJoint(int Yistar, int Ui1bound, int Ui2bound, double mu1i, double mu2i, double muZplusi, double muZminusi);
RcppExport SEXP _ZIPBayes_ZI_GenerateBigJoint(SEXP YistarSEXP, SEXP Ui1boundSEXP, SEXP Ui2boundSEXP, SEXP mu1iSEXP, SEXP mu2iSEXP, SEXP muZplusiSEXP, SEXP muZminusiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type Yistar(YistarSEXP);
    Rcpp::traits::input_parameter< int >::type Ui1bound(Ui1boundSEXP);
    Rcpp::traits::input_parameter< int >::type Ui2bound(Ui2boundSEXP);
    Rcpp::traits::input_parameter< double >::type mu1i(mu1iSEXP);
    Rcpp::traits::input_parameter< double >::type mu2i(mu2iSEXP);
    Rcpp::traits::input_parameter< double >::type muZplusi(muZplusiSEXP);
    Rcpp::traits::input_parameter< double >::type muZminusi(muZminusiSEXP);
    rcpp_result_gen = Rcpp::wrap(ZI_GenerateBigJoint(Yistar, Ui1bound, Ui2bound, mu1i, mu2i, muZplusi, muZminusi));
    return rcpp_result_gen;
END_RCPP
}
// ZI_GenerateJoint
NumericVector ZI_GenerateJoint(int Yistar, int minusbound, double phii, double mu2i, double muZplusi, double muZminusi);
RcppExport SEXP _ZIPBayes_ZI_GenerateJoint(SEXP YistarSEXP, SEXP minusboundSEXP, SEXP phiiSEXP, SEXP mu2iSEXP, SEXP muZplusiSEXP, SEXP muZminusiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type Yistar(YistarSEXP);
    Rcpp::traits::input_parameter< int >::type minusbound(minusboundSEXP);
    Rcpp::traits::input_parameter< double >::type phii(phiiSEXP);
    Rcpp::traits::input_parameter< double >::type mu2i(mu2iSEXP);
    Rcpp::traits::input_parameter< double >::type muZplusi(muZplusiSEXP);
    Rcpp::traits::input_parameter< double >::type muZminusi(muZminusiSEXP);
    rcpp_result_gen = Rcpp::wrap(ZI_GenerateJoint(Yistar, minusbound, phii, mu2i, muZplusi, muZminusi));
    return rcpp_result_gen;
END_RCPP
}
// ZI_GeneratePoiPar_Binary
NumericVector ZI_GeneratePoiPar_Binary(NumericVector Par, NumericMatrix covariates, NumericVector Outcome, NumericVector priorgamma);
RcppExport SEXP _ZIPBayes_ZI_GeneratePoiPar_Binary(SEXP ParSEXP, SEXP covariatesSEXP, SEXP OutcomeSEXP, SEXP priorgammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Par(ParSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type covariates(covariatesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Outcome(OutcomeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type priorgamma(priorgammaSEXP);
    rcpp_result_gen = Rcpp::wrap(ZI_GeneratePoiPar_Binary(Par, covariates, Outcome, priorgamma));
    return rcpp_result_gen;
END_RCPP
}
// ZI_GenerateU1
NumericVector ZI_GenerateU1(NumericVector Y, NumericVector U2, NumericVector mu1);
RcppExport SEXP _ZIPBayes_ZI_GenerateU1(SEXP YSEXP, SEXP U2SEXP, SEXP mu1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type U2(U2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu1(mu1SEXP);
    rcpp_result_gen = Rcpp::wrap(ZI_GenerateU1(Y, U2, mu1));
    return rcpp_result_gen;
END_RCPP
}
// ZI_GenerateU2
NumericVector ZI_GenerateU2(NumericVector Y, NumericVector U1, NumericVector mu2);
RcppExport SEXP _ZIPBayes_ZI_GenerateU2(SEXP YSEXP, SEXP U1SEXP, SEXP mu2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type U1(U1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu2(mu2SEXP);
    rcpp_result_gen = Rcpp::wrap(ZI_GenerateU2(Y, U1, mu2));
    return rcpp_result_gen;
END_RCPP
}
// ZI_GenerateV
NumericVector ZI_GenerateV(NumericVector Y, NumericVector Zminus, NumericVector muZminus, NumericMatrix covminus);
RcppExport SEXP _ZIPBayes_ZI_GenerateV(SEXP YSEXP, SEXP ZminusSEXP, SEXP muZminusSEXP, SEXP covminusSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Zminus(ZminusSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type muZminus(muZminusSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type covminus(covminusSEXP);
    rcpp_result_gen = Rcpp::wrap(ZI_GenerateV(Y, Zminus, muZminus, covminus));
    return rcpp_result_gen;
END_RCPP
}
// ZI_GenerateZpZmJoint
NumericVector ZI_GenerateZpZmJoint(int Yistar, int Yi, double muZplusi, double muZminusi);
RcppExport SEXP _ZIPBayes_ZI_GenerateZpZmJoint(SEXP YistarSEXP, SEXP YiSEXP, SEXP muZplusiSEXP, SEXP muZminusiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type Yistar(YistarSEXP);
    Rcpp::traits::input_parameter< int >::type Yi(YiSEXP);
    Rcpp::traits::input_parameter< double >::type muZplusi(muZplusiSEXP);
    Rcpp::traits::input_parameter< double >::type muZminusi(muZminusiSEXP);
    rcpp_result_gen = Rcpp::wrap(ZI_GenerateZpZmJoint(Yistar, Yi, muZplusi, muZminusi));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ZIPBayes_ZI_GenerateAlphaMNMetro", (DL_FUNC) &_ZIPBayes_ZI_GenerateAlphaMNMetro, 7},
    {"_ZIPBayes_ZI_GenerateBetaMetro", (DL_FUNC) &_ZIPBayes_ZI_GenerateBetaMetro, 5},
    {"_ZIPBayes_ZI_GenerateBigJoint", (DL_FUNC) &_ZIPBayes_ZI_GenerateBigJoint, 7},
    {"_ZIPBayes_ZI_GenerateJoint", (DL_FUNC) &_ZIPBayes_ZI_GenerateJoint, 6},
    {"_ZIPBayes_ZI_GeneratePoiPar_Binary", (DL_FUNC) &_ZIPBayes_ZI_GeneratePoiPar_Binary, 4},
    {"_ZIPBayes_ZI_GenerateU1", (DL_FUNC) &_ZIPBayes_ZI_GenerateU1, 3},
    {"_ZIPBayes_ZI_GenerateU2", (DL_FUNC) &_ZIPBayes_ZI_GenerateU2, 3},
    {"_ZIPBayes_ZI_GenerateV", (DL_FUNC) &_ZIPBayes_ZI_GenerateV, 4},
    {"_ZIPBayes_ZI_GenerateZpZmJoint", (DL_FUNC) &_ZIPBayes_ZI_GenerateZpZmJoint, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_ZIPBayes(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
