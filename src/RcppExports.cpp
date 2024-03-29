// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// model_landsepi
void model_landsepi(Rcpp::List time_param, Rcpp::NumericVector area_vector, Rcpp::NumericMatrix rotation_matrix, Rcpp::NumericMatrix croptypes_cultivars_prop, Rcpp::List dispersal, Rcpp::List inits, int seed, Rcpp::List cultivars_param, Rcpp::List basic_patho_param, Rcpp::List genes_param, Rcpp::List treatment_param);
RcppExport SEXP _landsepi_model_landsepi(SEXP time_paramSEXP, SEXP area_vectorSEXP, SEXP rotation_matrixSEXP, SEXP croptypes_cultivars_propSEXP, SEXP dispersalSEXP, SEXP initsSEXP, SEXP seedSEXP, SEXP cultivars_paramSEXP, SEXP basic_patho_paramSEXP, SEXP genes_paramSEXP, SEXP treatment_paramSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type time_param(time_paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type area_vector(area_vectorSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type rotation_matrix(rotation_matrixSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type croptypes_cultivars_prop(croptypes_cultivars_propSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type dispersal(dispersalSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type inits(initsSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type cultivars_param(cultivars_paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type basic_patho_param(basic_patho_paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type genes_param(genes_paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type treatment_param(treatment_paramSEXP);
    model_landsepi(time_param, area_vector, rotation_matrix, croptypes_cultivars_prop, dispersal, inits, seed, cultivars_param, basic_patho_param, genes_param, treatment_param);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_landsepi_model_landsepi", (DL_FUNC) &_landsepi_model_landsepi, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_landsepi(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
