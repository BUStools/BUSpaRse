// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// EC2gene_export
Rcpp::List EC2gene_export(Rcpp::DataFrame tr2g, std::string kallisto_out_path, bool verbose);
RcppExport SEXP _BUSpaRse_EC2gene_export(SEXP tr2gSEXP, SEXP kallisto_out_pathSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type tr2g(tr2gSEXP);
    Rcpp::traits::input_parameter< std::string >::type kallisto_out_path(kallisto_out_pathSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(EC2gene_export(tr2g, kallisto_out_path, verbose));
    return rcpp_result_gen;
END_RCPP
}
// fill_cell_gene
List fill_cell_gene(std::string fn, std::string kallisto_out_path, DataFrame tr2g, int est_ncells, int est_ngenes, std::vector<std::string> whitelist, bool gene_count, bool tcc, bool single_gene, bool verbose, int progress_unit);
RcppExport SEXP _BUSpaRse_fill_cell_gene(SEXP fnSEXP, SEXP kallisto_out_pathSEXP, SEXP tr2gSEXP, SEXP est_ncellsSEXP, SEXP est_ngenesSEXP, SEXP whitelistSEXP, SEXP gene_countSEXP, SEXP tccSEXP, SEXP single_geneSEXP, SEXP verboseSEXP, SEXP progress_unitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fn(fnSEXP);
    Rcpp::traits::input_parameter< std::string >::type kallisto_out_path(kallisto_out_pathSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type tr2g(tr2gSEXP);
    Rcpp::traits::input_parameter< int >::type est_ncells(est_ncellsSEXP);
    Rcpp::traits::input_parameter< int >::type est_ngenes(est_ngenesSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type whitelist(whitelistSEXP);
    Rcpp::traits::input_parameter< bool >::type gene_count(gene_countSEXP);
    Rcpp::traits::input_parameter< bool >::type tcc(tccSEXP);
    Rcpp::traits::input_parameter< bool >::type single_gene(single_geneSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type progress_unit(progress_unitSEXP);
    rcpp_result_gen = Rcpp::wrap(fill_cell_gene(fn, kallisto_out_path, tr2g, est_ncells, est_ngenes, whitelist, gene_count, tcc, single_gene, verbose, progress_unit));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BUSpaRse_EC2gene_export", (DL_FUNC) &_BUSpaRse_EC2gene_export, 3},
    {"_BUSpaRse_fill_cell_gene", (DL_FUNC) &_BUSpaRse_fill_cell_gene, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_BUSpaRse(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
