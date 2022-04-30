#include <Rcpp.h>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>



//[[Rcpp::plugins("cpp11")]]

using namespace Rcpp;

template<typename T> using vector2d = std::vector<std::vector<T>>;
template<typename T> using vector3d = std::vector<std::vector<std::vector<T>>>;



//' Below function filters positions so that they're all at least `min_dist`
//' from each other, and it does this within sequences.
//[[Rcpp::export]]
LogicalVector min_dist_filter(CharacterVector seq,
                              IntegerVector pos,
                              const int& min_dist) {
    size_t n_pos = pos.size();
    if (seq.size() != n_pos) stop("seq.size() != pos.size()");
    LogicalVector out(n_pos, 0);
    if (n_pos == 0) return out;
    out[0] = true;
    int lpos = pos[0];
    String lseq = seq[0];
    size_t i = 0;
    while (true) {
        while (i < n_pos && lseq == seq[i] && (pos[i] - lpos) < min_dist) {
            i++;
        }
        if (i >= n_pos) break;
        out[i] = true;
        lpos = pos[i];
        if (lseq != seq[i]) lseq = seq[i];
        i++;
    }
    return out;
}



//' Similar to above, but it instead groups positions using an integer vector.
//[[Rcpp::export]]
IntegerVector min_dist_grouper(StringVector seq,
                               IntegerVector pos,
                               const int& min_dist) {
    size_t n_pos = pos.size();
    if (seq.size() != n_pos) stop("seq.size() != pos.size()");
    IntegerVector out(n_pos);
    if (n_pos == 0) return out;
    int grp = 1;
    out[0] = grp;
    int lpos = pos[0];
    String lseq = seq[0];
    for (size_t i = 1; i < n_pos; i++) {
        if (lseq != seq[i]) {
            grp++;
            lseq = seq[i];
            lpos = pos[i];
        } else if ((pos[i] - lpos) >= min_dist) {
            grp++;
            lpos = pos[i];
        }
        out[i] = grp;
    }
    return out;
}

//' Similar to above, but it instead provides starts and ends for groups.
//[[Rcpp::export]]
IntegerMatrix min_dist_group_matrix(StringVector seq,
                                    IntegerVector pos,
                                    const int& min_dist) {
    size_t n_pos = pos.size();
    if (seq.size() != n_pos) stop("seq.size() != pos.size()");
    if (n_pos == 0) stop("seq.size() == 0");

    std::vector<int> starts;
    std::vector<int> ends;
    starts.reserve(n_pos);
    ends.reserve(n_pos);

    starts.push_back(1);
    int lpos = pos[0];
    String lseq = seq[0];
    for (size_t i = 1; i < n_pos; i++) {
        if (lseq != seq[i]) {
            ends.push_back(i);
            starts.push_back(i+1);
            lseq = seq[i];
            lpos = pos[i];
        } else if ((pos[i] - lpos) >= min_dist) {
            ends.push_back(i);
            starts.push_back(i+1);
            lpos = pos[i];
        }
    }
    ends.push_back(n_pos);

    size_t n_grps = starts.size();
    if (n_grps != ends.size()) {
        Rcout << "starts size = " << starts.size() << std::endl;
        Rcout << "ends size = " << ends.size() << std::endl;
        stop("starts.size() != ends.size()");
    }

    IntegerMatrix out(n_grps, 2);
    for (size_t i = 0; i < n_grps; i++) {
        out( i , 0 ) = starts[i];
        out( i , 1 ) = ends[i];
    }


    return out;
}





//' Get the indices for the alleles that have counts > 0 in these pools.
//' Change biallelic to `false` if there are >2 of these bc MakeTree
//' only allows biallelic SNPs.
void get_allele_indices(bool& biallelic, size_t& ref_idx, size_t& alt_idx,
                        const vector3d<unsigned>& allele_counts,
                        const size_t& n_pools,
                        const size_t& allele_idx) {
    biallelic = true;
    ref_idx = 99;
    alt_idx = 99;
    std::vector<unsigned> sums(6, 0);
    unsigned n_nonzeros = 0;
    for (unsigned k = 0; k < sums.size(); k++) {
        for (size_t i = 0; i < n_pools; i++) {
            sums[k] += allele_counts[i][allele_idx][k];
        }
        if (sums[k] > 0) {
            n_nonzeros++;
            if (ref_idx == 99) {
                ref_idx = k;
            } else if (alt_idx == 99) {
                alt_idx = k;
            }
        }
        if (n_nonzeros > 2) {
            biallelic = false;
            break;
        }
    }
    if (n_nonzeros < 2) biallelic = false;

    return;
}


void make_strings(vector2d<std::string>& out,
                  const bool& biallelic,
                  const size_t& ref_idx,
                  const size_t& alt_idx,
                  const vector3d<unsigned>& allele_counts,
                  const size_t& n_pools,
                  const size_t& allele_idx) {

    std::string s = "";

    if (biallelic) {
        if (ref_idx > 5 || alt_idx > 5) {
            Rcout << ref_idx << ' ' << alt_idx << std::endl;
            stop("improper indices");
        }
        for (size_t i = 0; i < n_pools; i++) {
            const std::vector<unsigned>& ac(allele_counts[i][allele_idx]);
            s = std::to_string(ac[ref_idx]);
            s += ',';
            s += std::to_string(ac[alt_idx]);
            out[i].push_back(s);
        }
    } else {
        for (size_t i = 0; i < n_pools; i++) {
            out[i].push_back(s);
        }
    }

    return;
}


//[[Rcpp::export]]
vector2d<std::string> maketree_format(const vector3d<unsigned>& allele_counts) {

    // allele_counts[pool][allele][base]

    size_t n_pools = allele_counts.size();
    size_t n_alleles = allele_counts[0].size();

    vector2d<std::string> out(n_pools);
    if (n_pools == 0) return out;
    if (n_alleles == 0) return out;
    for (std::vector<std::string>& x : out) x.reserve(n_alleles);

    bool biallelic;
    size_t ref_idx, alt_idx;

    for (size_t j = 0; j < n_alleles; j++) {
        get_allele_indices(biallelic, ref_idx, alt_idx, allele_counts,
                           n_pools, j);
        make_strings(out, biallelic, ref_idx, alt_idx, allele_counts,
                     n_pools, j);
    }

    return out;
}






//' Get allele frequencies from string outputs from SNAPE
//[[Rcpp::export]]
std::vector<double> get_snape_af(const std::vector<std::string>& strings) {

    std::vector<double> out(strings.size());

    // Calculate number of sub-strings separated by ':':
    size_t nc = 0;
    for (const char& c : strings[0]) if (c == ':') nc++;
    nc++;

    std::string tmp0;
    std::vector<std::string> tmp1(nc);

    for (size_t i = 0; i < strings.size(); i++) {
        const std::string& s(strings[i]);
        std::istringstream iss(s);
        size_t j = 0;
        while (std::getline(iss, tmp0, ':')) {
            tmp1[j] = tmp0;
            j++;
        }
        out[i] = std::stod(tmp1.back());
    }

    return out;
}

// Get allele counts from gSYNC outputs
//[[Rcpp::export]]
List split_sync_strings(const std::vector<std::string>& strings) {

    List out(strings.size());

    // Calculate number of sub-strings separated by ':':
    size_t nc = 0;
    for (const char& c : strings[0]) if (c == ':') nc++;
    nc++;

    std::string tmp0;
    std::vector<int> tmp1(nc);

    for (size_t i = 0; i < strings.size(); i++) {
        const std::string& s(strings[i]);
        std::istringstream iss(s);
        size_t j = 0;
        while (std::getline(iss, tmp0, ':')) {
            tmp1[j] = std::stoi(tmp0);
            j++;
        }
        out[i] = tmp1;
    }

    return out;
}


//' This function calculates the number of alleles at each locus
//[[Rcpp::export]]
std::vector<int> count_alleles(const vector3d<unsigned>& allele_counts) {

    size_t n_pools = allele_counts.size();
    size_t n_alleles = allele_counts[0].size();
    if (n_pools == 0) return std::vector<int>();
    if (n_alleles == 0) return std::vector<int>();

    std::vector<int> out;
    out.reserve(n_alleles);

    for (size_t j = 0; j < n_alleles; j++) {
        out.push_back(0);
        int& n_nonzeros(out.back());
        for (size_t k = 0; k < 6; k++) {
            unsigned sum_k = 0;
            for (size_t i = 0; i < n_pools; i++) {
                sum_k += allele_counts[i][j][k];
            }
            if (sum_k > 0) n_nonzeros++;
        }
    }

    return out;
}


//' This function determines which alleles are non-zero across all pools.
//' Returns two-column integer matrix.
//'
//[[Rcpp::export]]
IntegerMatrix nonzero_alleles(const vector3d<unsigned>& allele_counts) {

    size_t n_pools = allele_counts.size();
    size_t n_alleles = allele_counts[0].size();

    IntegerMatrix out(n_alleles, 2);

    if (n_pools == 0 || n_alleles == 0) return out;

    bool biallelic;
    size_t ref_idx, alt_idx;

    for (size_t j = 0; j < n_alleles; j++) {
        get_allele_indices(biallelic, ref_idx, alt_idx, allele_counts,
                           n_pools, j);
        if (! biallelic) stop("nonzero_alleles should only contain biallelic loci");
        // Because we're returning this to R:
        ref_idx++;
        alt_idx++;
        out(j, 0) = ref_idx;
        out(j, 1) = alt_idx;
    }

    return out;
}





inline int one_which_allele(const std::vector<unsigned>& cnt,
                            const double& frq) {

    double cnt_sum = static_cast<double>(std::accumulate(cnt.begin(),
                                                         cnt.end(), 0));
    double cnt_frq;
    std::vector<double> cnt_frq_diffs;
    cnt_frq_diffs.reserve(cnt.size());
    for (const unsigned& x : cnt) {
        cnt_frq = static_cast<double>(x) / cnt_sum;
        cnt_frq_diffs.push_back(std::abs(cnt_frq - frq));
    }

    int mdi = std::min_element(cnt_frq_diffs.begin(), cnt_frq_diffs.end()) -
        cnt_frq_diffs.begin();
    // Because output should be 1-based indices.
    mdi++;
    return mdi;

}


//' This function determines which allele the frequencies must be referring to.
//' This is a sanity check to make sure they don't vary between pools.
//'
//' Only recommended for biallelic loci!
//'
//[[Rcpp::export]]
IntegerMatrix which_allele(const vector3d<unsigned>& allele_counts,
                           NumericMatrix allele_freqs) {

    size_t n_pools = allele_counts.size();
    size_t n_alleles = allele_counts[0].size();

    IntegerMatrix out(n_alleles, n_pools);

    if (n_pools == 0) return out;
    if (n_alleles == 0) return out;
    if (allele_freqs.ncol() != n_pools) {
        stop("allele_freqs.ncol() != n_pools");
    }
    if (allele_freqs.nrow() != n_alleles) {
        stop("allele_freqs.nrow() != n_alleles");
    }

    int tmp;

    for (size_t i = 0; i < n_pools; i++) {
        for (size_t j = 0; j < n_alleles; j++) {
            tmp = one_which_allele(allele_counts[i][j],
                                   allele_freqs(j, i));
            out(j, i) = tmp;
        }
    }

    return out;
}



