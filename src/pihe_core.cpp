// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// Per-event-time summary for the PiHE two-pass computation. Defined at file
// scope so the pointer-based implementation can take a reusable buffer of it
// from the caller (the fused analysis loop) and avoid reallocating per cell.
struct PiheEvSummary {
  double n1k;
  double n0k;
  double OTk;
  double Ok;
};

// Forward declaration of the pointer-based implementation (defined below).
void pihe_core_impl(const double*, const int*, const int*, const int*, int,
                    std::vector<PiheEvSummary>&, double*);

//' Core PiHE hazard ratio computation (C++ backend)
//'
//' @description
//' Internal C++ function that computes the Pike-Halley Estimator (PiHE) for
//' the hazard ratio. Accepts pooled sorted vectors plus an integer group
//' indicator, performs group splitting and the two-pointer merge scan
//' entirely in C++, and returns the quantities needed for the Halley
//' correction and Wald interval. Not intended to be called directly by
//' users; use \code{coxph_fast()} instead.
//'
//' Implementation notes: Pass 1 accumulates pooled scalars
//' (O_T, O_C, E_T, E_C) and stores per-distinct-event-time summaries in a
//' single struct array (better cache locality than four parallel vectors)
//' with a single \code{reserve(n)} call (no event-count prepass over the
//' input). Pass 2 walks the saved summaries once at the Pike anchor
//' theta_0 to compute U_0, I_0, and J_0.
//'
//' When \code{strata_sorted} is supplied (non-empty), the data are assumed
//' sorted by \code{(strata_sorted, time_sorted)}: Pass 1 resets the at-risk
//' counts at each stratum boundary so that risk sets never cross strata,
//' while \code{O_T}, \code{O_C}, \code{E_T}, \code{E_C} keep accumulating
//' across strata, giving the pooled Pike anchor. Pass 2 is unaffected, since
//' the saved per-event summaries are already stratum-local.
//'
//' @param time_sorted A numeric vector of pooled follow-up times sorted in
//'   ascending order (within stratum, when \code{strata_sorted} is given).
//' @param event_sorted An integer vector of event indicators (1 = event,
//'   0 = censored), aligned with \code{time_sorted}.
//' @param j_sorted An integer vector of group indicators (1 = treatment,
//'   0 = control), aligned with \code{time_sorted}.
//' @param strata_sorted An integer vector of stratum codes aligned with
//'   \code{time_sorted}, assumed sorted so that rows sharing a stratum are
//'   contiguous. Pass a zero-length vector to indicate no stratification.
//'
//' @return A numeric vector of length 4: \code{c(theta_0, U_0, I_0, J_0)}.
//'   Returns a length-4 vector of \code{NA_real_} when the estimate cannot
//'   be computed.
//'
//' @keywords internal
// [[Rcpp::export]]
NumericVector pihe_core(
    const NumericVector& time_sorted,
    const IntegerVector& event_sorted,
    const IntegerVector& j_sorted,
    const IntegerVector& strata_sorted
) {
  const int n = time_sorted.size();
  std::vector<PiheEvSummary> ev;
  double out[4];
  const int* strata_ptr = strata_sorted.size() > 0 ? strata_sorted.begin()
                                                    : nullptr;
  pihe_core_impl(time_sorted.begin(), event_sorted.begin(),
                 j_sorted.begin(), strata_ptr, n, ev, out);
  return NumericVector::create(out[0], out[1], out[2], out[3]);
}

// Pointer-based implementation with external linkage. Writes
// theta_0, U_0, I_0, J_0 into out[0..3], or four NA values when the estimate
// cannot be computed. The per-event-time summaries are stored in the
// caller-supplied buffer ev, which is cleared on entry and reused across cells
// to avoid a reserve(n) allocation per call. Algorithm identical to the
// exported wrapper above.
//
// strata_sorted may be nullptr, meaning "no stratification" (a single
// implicit stratum spanning the whole input) -- this collapses the outer
// stratum-segmenting loop below to a single iteration over [0, n), reducing
// to the original unstratified algorithm exactly.
void pihe_core_impl(
    const double* time_sorted,
    const int* event_sorted,
    const int* j_sorted,
    const int* strata_sorted,
    int n,
    std::vector<PiheEvSummary>& ev,
    double* out
) {
  out[0] = NA_REAL; out[1] = NA_REAL; out[2] = NA_REAL; out[3] = NA_REAL;
  if (n == 0) return;

  ev.clear();
  if ((int) ev.capacity() < n) ev.reserve(n);

  double O_T = 0.0, O_C = 0.0, E_T = 0.0, E_C = 0.0;

  // Outer loop: segment the (strata, time)-sorted input into contiguous
  // per-stratum blocks. With strata_sorted == nullptr this runs exactly once
  // over the full range [0, n), identical to the original single-stratum
  // scan.
  int strat_start = 0;
  while (strat_start < n) {
    int strat_end = strat_start;
    if (strata_sorted == nullptr) {
      strat_end = n;
    } else {
      const int s = strata_sorted[strat_start];
      while (strat_end < n && strata_sorted[strat_end] == s) ++strat_end;
    }

    // Stratum-local at-risk counts: totals per group within this segment
    int n1 = 0, n0 = 0;
    for (int k = strat_start; k < strat_end; ++k) {
      if (j_sorted[k] == 1) ++n1; else ++n0;
    }

    // Pass 1: single left-to-right scan over this stratum's sorted data
    int i = strat_start;
    while (i < strat_end) {
      const double t = time_sorted[i];

      // Consume the tied block at time t, bounded to this stratum's segment
      // so a tied block can never cross a stratum boundary
      int d1 = 0, d0 = 0, c1 = 0, c0 = 0;
      int j = i;
      while (j < strat_end && time_sorted[j] == t) {
        if (j_sorted[j] == 1) {
          ++c1;
          if (event_sorted[j] == 1) ++d1;
        } else {
          ++c0;
          if (event_sorted[j] == 1) ++d0;
        }
        ++j;
      }

      const int d  = d1 + d0;
      const int nj = n1 + n0;

      if (d > 0 && nj > 0) {
        const double dn1 = (double)n1;
        const double dn0 = (double)n0;
        const double dnj = (double)nj;
        const double dd  = (double)d;
        const double dd1 = (double)d1;

        O_T += dd1;
        O_C += dd - dd1;
        E_T += dd * dn1 / dnj;
        E_C += dd * dn0 / dnj;

        ev.push_back(PiheEvSummary{dn1, dn0, dd1, dd});
      }

      // Decrement at-risk counts by block size
      n1 -= c1;
      n0 -= c0;
      i   = j;
    }

    strat_start = strat_end;
  }

  if (O_T == 0.0 || O_C == 0.0 || E_T == 0.0 || E_C == 0.0) return;

  // Pike anchor
  const double theta_0 = (O_T * E_C) / (O_C * E_T);

  // Pass 2: walk saved summaries to compute U_0, I_0, J_0 at theta_0
  double U_0 = 0.0, I_0 = 0.0, J_0 = 0.0;
  const std::size_t K = ev.size();
  for (std::size_t k = 0; k < K; ++k) {
    const PiheEvSummary& e = ev[k];
    const double denom = e.n0k + e.n1k * theta_0;
    if (denom == 0.0) return;

    const double p_k = e.n1k * theta_0 / denom;
    const double q_k = 1.0 - p_k;
    const double pq  = p_k * q_k;
    const double Opq = e.Ok * pq;

    U_0 += e.OTk - e.Ok * p_k;
    I_0 += Opq;
    J_0 += Opq * (1.0 - 2.0 * p_k);
  }

  if (!std::isfinite(I_0) || I_0 == 0.0) return;

  out[0] = theta_0; out[1] = U_0; out[2] = I_0; out[3] = J_0;
}
