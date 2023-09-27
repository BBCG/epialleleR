// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rcpp_call_methylation_genome
Rcpp::List rcpp_call_methylation_genome(std::string in_fn, std::string out_fn, Rcpp::List& genome, std::string tag, int nthreads);
RcppExport SEXP _epialleleR_rcpp_call_methylation_genome(SEXP in_fnSEXP, SEXP out_fnSEXP, SEXP genomeSEXP, SEXP tagSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type in_fn(in_fnSEXP);
    Rcpp::traits::input_parameter< std::string >::type out_fn(out_fnSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type genome(genomeSEXP);
    Rcpp::traits::input_parameter< std::string >::type tag(tagSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_call_methylation_genome(in_fn, out_fn, genome, tag, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_check_bam
Rcpp::List rcpp_check_bam(std::string fn);
RcppExport SEXP _epialleleR_rcpp_check_bam(SEXP fnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fn(fnSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_check_bam(fn));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_cx_report
Rcpp::DataFrame rcpp_cx_report(Rcpp::DataFrame& df, Rcpp::LogicalVector& pass, std::string ctx);
RcppExport SEXP _epialleleR_rcpp_cx_report(SEXP dfSEXP, SEXP passSEXP, SEXP ctxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector& >::type pass(passSEXP);
    Rcpp::traits::input_parameter< std::string >::type ctx(ctxSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_cx_report(df, pass, ctx));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_extract_patterns
Rcpp::DataFrame rcpp_extract_patterns(Rcpp::DataFrame& df, unsigned int target_rname, unsigned int target_start, unsigned int target_end, signed int min_overlap, std::string& ctx, double min_ctx_freq, bool clip, unsigned int reverse_offset, Rcpp::IntegerVector& hlght);
RcppExport SEXP _epialleleR_rcpp_extract_patterns(SEXP dfSEXP, SEXP target_rnameSEXP, SEXP target_startSEXP, SEXP target_endSEXP, SEXP min_overlapSEXP, SEXP ctxSEXP, SEXP min_ctx_freqSEXP, SEXP clipSEXP, SEXP reverse_offsetSEXP, SEXP hlghtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type target_rname(target_rnameSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type target_start(target_startSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type target_end(target_endSEXP);
    Rcpp::traits::input_parameter< signed int >::type min_overlap(min_overlapSEXP);
    Rcpp::traits::input_parameter< std::string& >::type ctx(ctxSEXP);
    Rcpp::traits::input_parameter< double >::type min_ctx_freq(min_ctx_freqSEXP);
    Rcpp::traits::input_parameter< bool >::type clip(clipSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type reverse_offset(reverse_offsetSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type hlght(hlghtSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_extract_patterns(df, target_rname, target_start, target_end, min_overlap, ctx, min_ctx_freq, clip, reverse_offset, hlght));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_fep
std::vector<double> rcpp_fep(Rcpp::DataFrame& df, std::vector<std::string> colnames);
RcppExport SEXP _epialleleR_rcpp_fep(SEXP dfSEXP, SEXP colnamesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type colnames(colnamesSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_fep(df, colnames));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_get_base_freqs
Rcpp::NumericMatrix rcpp_get_base_freqs(Rcpp::DataFrame& df, std::vector<bool> pass, Rcpp::DataFrame& vcf);
RcppExport SEXP _epialleleR_rcpp_get_base_freqs(SEXP dfSEXP, SEXP passSEXP, SEXP vcfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< std::vector<bool> >::type pass(passSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame& >::type vcf(vcfSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_get_base_freqs(df, pass, vcf));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_get_xm_beta
std::vector<double> rcpp_get_xm_beta(Rcpp::DataFrame& df, std::string ctx_meth, std::string ctx_unmeth);
RcppExport SEXP _epialleleR_rcpp_get_xm_beta(SEXP dfSEXP, SEXP ctx_methSEXP, SEXP ctx_unmethSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< std::string >::type ctx_meth(ctx_methSEXP);
    Rcpp::traits::input_parameter< std::string >::type ctx_unmeth(ctx_unmethSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_get_xm_beta(df, ctx_meth, ctx_unmeth));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_match_amplicon
std::vector<int> rcpp_match_amplicon(Rcpp::DataFrame& df, Rcpp::DataFrame& bed, int tolerance);
RcppExport SEXP _epialleleR_rcpp_match_amplicon(SEXP dfSEXP, SEXP bedSEXP, SEXP toleranceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame& >::type bed(bedSEXP);
    Rcpp::traits::input_parameter< int >::type tolerance(toleranceSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_match_amplicon(df, bed, tolerance));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_match_capture
std::vector<int> rcpp_match_capture(Rcpp::DataFrame& df, Rcpp::DataFrame& bed, signed int min_overlap);
RcppExport SEXP _epialleleR_rcpp_match_capture(SEXP dfSEXP, SEXP bedSEXP, SEXP min_overlapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame& >::type bed(bedSEXP);
    Rcpp::traits::input_parameter< signed int >::type min_overlap(min_overlapSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_match_capture(df, bed, min_overlap));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_mhl_report
Rcpp::DataFrame rcpp_mhl_report(Rcpp::DataFrame& df, std::string ctx, int hmax);
RcppExport SEXP _epialleleR_rcpp_mhl_report(SEXP dfSEXP, SEXP ctxSEXP, SEXP hmaxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< std::string >::type ctx(ctxSEXP);
    Rcpp::traits::input_parameter< int >::type hmax(hmaxSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_mhl_report(df, ctx, hmax));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_read_bam_paired
Rcpp::DataFrame rcpp_read_bam_paired(std::string fn, int min_mapq, int min_baseq, bool skip_duplicates, int nthreads);
RcppExport SEXP _epialleleR_rcpp_read_bam_paired(SEXP fnSEXP, SEXP min_mapqSEXP, SEXP min_baseqSEXP, SEXP skip_duplicatesSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fn(fnSEXP);
    Rcpp::traits::input_parameter< int >::type min_mapq(min_mapqSEXP);
    Rcpp::traits::input_parameter< int >::type min_baseq(min_baseqSEXP);
    Rcpp::traits::input_parameter< bool >::type skip_duplicates(skip_duplicatesSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_read_bam_paired(fn, min_mapq, min_baseq, skip_duplicates, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_read_bam_single
Rcpp::DataFrame rcpp_read_bam_single(std::string fn, int min_mapq, int min_baseq, bool skip_duplicates, int nthreads);
RcppExport SEXP _epialleleR_rcpp_read_bam_single(SEXP fnSEXP, SEXP min_mapqSEXP, SEXP min_baseqSEXP, SEXP skip_duplicatesSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fn(fnSEXP);
    Rcpp::traits::input_parameter< int >::type min_mapq(min_mapqSEXP);
    Rcpp::traits::input_parameter< int >::type min_baseq(min_baseqSEXP);
    Rcpp::traits::input_parameter< bool >::type skip_duplicates(skip_duplicatesSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_read_bam_single(fn, min_mapq, min_baseq, skip_duplicates, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_read_genome
Rcpp::List rcpp_read_genome(std::string fn, int nthreads);
RcppExport SEXP _epialleleR_rcpp_read_genome(SEXP fnSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fn(fnSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_read_genome(fn, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_simulate_bam
int rcpp_simulate_bam(std::vector<std::string> header, Rcpp::DataFrame& fields, Rcpp::DataFrame& i_tags, Rcpp::DataFrame& s_tags, std::string out_fn);
RcppExport SEXP _epialleleR_rcpp_simulate_bam(SEXP headerSEXP, SEXP fieldsSEXP, SEXP i_tagsSEXP, SEXP s_tagsSEXP, SEXP out_fnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type header(headerSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame& >::type fields(fieldsSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame& >::type i_tags(i_tagsSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame& >::type s_tags(s_tagsSEXP);
    Rcpp::traits::input_parameter< std::string >::type out_fn(out_fnSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_simulate_bam(header, fields, i_tags, s_tags, out_fn));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_threshold_reads
std::vector<bool> rcpp_threshold_reads(Rcpp::DataFrame& df, std::string ctx_meth, std::string ctx_unmeth, std::string ooctx_meth, std::string ooctx_unmeth, unsigned int min_n_ctx, double min_ctx_meth_frac, double max_ooctx_meth_frac);
RcppExport SEXP _epialleleR_rcpp_threshold_reads(SEXP dfSEXP, SEXP ctx_methSEXP, SEXP ctx_unmethSEXP, SEXP ooctx_methSEXP, SEXP ooctx_unmethSEXP, SEXP min_n_ctxSEXP, SEXP min_ctx_meth_fracSEXP, SEXP max_ooctx_meth_fracSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< std::string >::type ctx_meth(ctx_methSEXP);
    Rcpp::traits::input_parameter< std::string >::type ctx_unmeth(ctx_unmethSEXP);
    Rcpp::traits::input_parameter< std::string >::type ooctx_meth(ooctx_methSEXP);
    Rcpp::traits::input_parameter< std::string >::type ooctx_unmeth(ooctx_unmethSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type min_n_ctx(min_n_ctxSEXP);
    Rcpp::traits::input_parameter< double >::type min_ctx_meth_frac(min_ctx_meth_fracSEXP);
    Rcpp::traits::input_parameter< double >::type max_ooctx_meth_frac(max_ooctx_meth_fracSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_threshold_reads(df, ctx_meth, ctx_unmeth, ooctx_meth, ooctx_unmeth, min_n_ctx, min_ctx_meth_frac, max_ooctx_meth_frac));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_epialleleR_rcpp_call_methylation_genome", (DL_FUNC) &_epialleleR_rcpp_call_methylation_genome, 5},
    {"_epialleleR_rcpp_check_bam", (DL_FUNC) &_epialleleR_rcpp_check_bam, 1},
    {"_epialleleR_rcpp_cx_report", (DL_FUNC) &_epialleleR_rcpp_cx_report, 3},
    {"_epialleleR_rcpp_extract_patterns", (DL_FUNC) &_epialleleR_rcpp_extract_patterns, 10},
    {"_epialleleR_rcpp_fep", (DL_FUNC) &_epialleleR_rcpp_fep, 2},
    {"_epialleleR_rcpp_get_base_freqs", (DL_FUNC) &_epialleleR_rcpp_get_base_freqs, 3},
    {"_epialleleR_rcpp_get_xm_beta", (DL_FUNC) &_epialleleR_rcpp_get_xm_beta, 3},
    {"_epialleleR_rcpp_match_amplicon", (DL_FUNC) &_epialleleR_rcpp_match_amplicon, 3},
    {"_epialleleR_rcpp_match_capture", (DL_FUNC) &_epialleleR_rcpp_match_capture, 3},
    {"_epialleleR_rcpp_mhl_report", (DL_FUNC) &_epialleleR_rcpp_mhl_report, 3},
    {"_epialleleR_rcpp_read_bam_paired", (DL_FUNC) &_epialleleR_rcpp_read_bam_paired, 5},
    {"_epialleleR_rcpp_read_bam_single", (DL_FUNC) &_epialleleR_rcpp_read_bam_single, 5},
    {"_epialleleR_rcpp_read_genome", (DL_FUNC) &_epialleleR_rcpp_read_genome, 2},
    {"_epialleleR_rcpp_simulate_bam", (DL_FUNC) &_epialleleR_rcpp_simulate_bam, 5},
    {"_epialleleR_rcpp_threshold_reads", (DL_FUNC) &_epialleleR_rcpp_threshold_reads, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_epialleleR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
