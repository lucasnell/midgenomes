#include <Rcpp.h>

#include <vector>
#include <string>
#include <deque>

#include <fstream>
#include <iostream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

/*
 Note that this also requires installation of boost
 because `gzip_decompressor()` isn't header-only
 */



#include <progress.hpp>  // for the progress bar

// [[Rcpp::depends(RcppProgress,BH)]]
// [[Rcpp::plugins(cpp11,openmp)]]

using namespace Rcpp;



inline void expand_path(std::string& file_name) {
    // Obtain environment containing function
    Environment base("package:base");
    // Make function callable from C++
    Function pe_r = base["path.expand"];
    // Call the function and receive its list output
    SEXP fai_file_exp = pe_r(file_name);
    file_name = as<std::string>(fai_file_exp);
    return;
}



//' For prettier long error messages.
//'
//'
inline void str_stop(const std::vector<std::string>& err_msg_vec) {
    std::string err_msg = "";
    for (const std::string& err : err_msg_vec) err_msg += err;
    throw(Rcpp::exception(err_msg.c_str(), false));
}


//'
//'
// [[Rcpp::export]]
std::vector<std::string> get_read_names(std::string in_fn,
                                        const uint32_t& n_reads) {

    std::vector<std::string> read_names;
    read_names.reserve(n_reads);

    expand_path(in_fn);

    Progress prog_bar(n_reads, true);

    // Input
    std::ifstream in_file(in_fn.c_str(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> in_buf;
    in_buf.push(boost::iostreams::gzip_decompressor());
    in_buf.push(in_file);
    std::istream in_stream(&in_buf);

    std::string line;
    std::string name_i;
    uint32_t i = 0;
    while (std::getline(in_stream, line)) {
        if (i > 0 && i % 1000 == 0) {
            Rcpp::checkUserInterrupt();
            prog_bar.increment(250);
        }
        if (i % 4 == 0) {
            std::string::size_type spc = line.find(' ', 2);
            if (spc == std::string::npos) {
                spc = line.size();
            } else {
                // Don't include the space:
                spc--;
            }
            name_i = line.substr(1, spc);
            read_names.push_back(name_i);
        }
        i++;
    }

    return read_names;

}



//'
//' `read_nums` is the read numbers of the reads you want in the output FASTQ.
//'
// [[Rcpp::export]]
void filter_fastq(std::string in_fn,
                  std::string out_fn,
                  std::deque<uint32_t> read_nums) {

    expand_path(in_fn);
    expand_path(out_fn);

    Progress prog_bar(read_nums.size(), true);

    // Input
    std::ifstream in_file(in_fn.c_str(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> in_buf;
    in_buf.push(boost::iostreams::gzip_decompressor());
    in_buf.push(in_file);
    std::istream in_stream(&in_buf);

    Rcpp::checkUserInterrupt();

    // Output
    std::ofstream out_file(out_fn.c_str(), std::ios_base::out | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::output> out_buf;
    out_buf.push(boost::iostreams::gzip_compressor());
    out_buf.push(out_file);
    std::ostream out_stream(&out_buf);

    Rcpp::checkUserInterrupt();

    // Iterate lines
    std::string line;
    line.reserve(100000);
    bool to_out = false;
    uint32_t found_reads = 0;
    uint32_t line_i = 0;
    uint32_t read_i = 0;
    while (std::getline(in_stream, line)) {
        if (read_nums.empty()) break;
        if (line_i % 1000 == 0) Rcpp::checkUserInterrupt();
        if (line_i % 4 == 0) {
            to_out = read_i == read_nums.front();
            read_i++;
            if (to_out) {
                found_reads++;
                read_nums.pop_front();
                prog_bar.increment();
            }
        }
        if (to_out) out_stream << line << std::endl;
        line_i++;
    }

    // Cleanup
    in_file.close();
    boost::iostreams::close(out_buf);
    out_file.close();

    if (read_nums.size() > 0) {
        Rcout << static_cast<uint32_t>(read_nums.size()) <<
            " reads left!\n";
        Rcout << "Found " << found_reads << " reads in total.\n";
    }

    return;


}

