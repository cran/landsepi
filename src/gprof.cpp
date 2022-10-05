/*
#if defined __linux__ && defined DEBUG

#include <Rcpp.h>
#include <gperftools/profiler.h>

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
SEXP start_profiler(SEXP str) {
  ProfilerStart(as<const char*>(str));
  return R_NilValue;
}

//' @export
// [[Rcpp::export]]
SEXP stop_profiler() {
  ProfilerStop();
  return R_NilValue;
}
#endif
*/
