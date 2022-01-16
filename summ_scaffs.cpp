// Provides some basic stats on scaffolds.

/*

 On linux, compile with
 g++ -std=c++11 summ_scaffs.cpp -o summ_scaffs

 On your mac, compile with
 g++ -std=c++11 -isysroot $(xcrun --show-sdk-path) summ_scaffs.cpp -o summ_scaffs

 */


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>



typedef unsigned long int uint32;


int main(int argc, char* argv[]){

    if (argc == 1) {
        std::cerr << "No file name provided." << std::endl;
        return 1;
    }
    if (argc > 2) {
        std::cerr << "Too many arguments provided." << std::endl;
        return 1;
    }

    std::cout << std::endl << std::endl;
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl << std::endl;
    std::cout << "File " << argv[1] << std::endl << std::endl;

    std::fstream newfile;
    newfile.open(argv[1], std::ios::in);

    if (newfile.is_open()){
        std::string line;
        // some of these programs produce single-line FASTA files, which
        // results in very long lines :/
        line.reserve(5000000);
        std::vector<uint32> sizes;
        sizes.reserve(10000);
        uint32 total_size = 0;
        uint32 total_N = 0;
        while (std::getline(newfile, line)){
            if (line[0] == '>') {
                sizes.push_back(0);
            } else {
                if (sizes.empty()) {
                    std::cerr << "File does not start with header" << std::endl;
                    return 1;
                }
                sizes.back() += line.size();
                total_size += line.size();
                for (uint32 i = 0; i < line.size(); i++) {
                    if (line[i] == 'N' || line[i] == 'n') total_N++;
                }
            }
        }
        newfile.close();

        // Sort in decreasing order:
        std::sort(sizes.rbegin(), sizes.rend());
        // Calculate N50:
        uint32 n50_threshold = static_cast<uint32>(
            std::round(static_cast<double>(total_size) / 2.0));
        uint32 cum_sum = 0;
        uint32 i = 0;
        for (; i < sizes.size(); i++) {
            cum_sum += sizes[i];
            if (cum_sum >= n50_threshold) break;
        }
        uint32 n50 = sizes[i];

        std::cout << "size = " << total_size << std::endl;
        std::cout << sizes.size() << " contigs" << std::endl;
        std::cout << "N50 = " << n50 << std::endl;
        std::cout << "min = " << sizes.back() << std::endl;
        std::cout << "max = " << sizes.front() << std::endl;
        std::cout << "total N = " << total_N << std::endl;

    } else {

        std::cerr << "Cannot open file" << argv[1] << '!' << std::endl;

        return 1;

    }

    std::cout << std::endl << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
    std::cout << std::endl << std::endl << std::endl;

    return 0;

}
