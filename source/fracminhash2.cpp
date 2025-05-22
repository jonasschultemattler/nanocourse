#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include "util.hpp"
#include "hyperloglog2.hpp"


static double jaccard(const double containment, const uint64_t size_a, const uint64_t size_b) {
    return containment*size_a/(size_a+size_b-containment*size_a);
}


// Note: How can you compute the intersection of two vectors in C++?
// One possibility: use std::set_intersection on sorted vectors
// e.g. std::sort(A.begin(), A.end());
// std::vector<uint64_t> intersection;
// std::set_intersection(A.begin(), A.end(), B.begin(), B.end(), back_inserter(intersection));


void fracminsketch(const std::filesystem::path &filepath, std::vector<uint64_t> &frac_sketch,
                   const double s, uint64_t (*hashFunc)(uint64_t)=wyhash)
{
    auto fin = seqan3::sequence_file_input<my_traits>{filepath};
    auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k});
    for(auto & record : fin) {
        for(auto && kmer : record.sequence() | kmer_view) {
            uint64_t hash = hashFunc(kmer);
            if(hash <= UINT64_MAX*s)
                frac_sketch.push_back(hash);
        }
    }
    std::sort(frac_sketch.begin(), frac_sketch.end());
}


double fracMinHash(const std::vector<uint64_t> &frac_sketch_a, const std::vector<uint64_t> &frac_sketch_b, const uint64_t size_a,
                   const double s)
{
    std::vector<uint64_t> intersection;
    std::set_intersection(frac_sketch_a.begin(), frac_sketch_a.end(), frac_sketch_b.begin(),
                          frac_sketch_b.end(), back_inserter(intersection));
    
    return (double) intersection.size()/(frac_sketch_a.size()*(1-std::pow(1-s, size_a)));
}


void fracminhash_similarities(const std::vector<std::filesystem::path> &filepaths, double matrix[n][n],
    const double s=0.1)
{
    uint64_t sizes[n];
    for(int i = 0; i < n; i++) {
        sizes[i] = hyperloglog(filepaths[i]);
    }
    for(int i = 0; i < n; i++) {
        matrix[i][i] = 1;
        std::vector<uint64_t> fracsketch_i;
        fracminsketch(filepaths[i], fracsketch_i, s);
        for(int j = i+1; j < n; j++) {
            std::vector<uint64_t> fracsketch_j;
            fracminsketch(filepaths[j], fracsketch_j, s);
            const double containment = fracMinHash(fracsketch_i, fracsketch_j, sizes[i], s);
            matrix[i][j] = matrix[j][i] = jaccard(containment, sizes[i], sizes[j]);
        }
    }
}


int main(int argc, char** argv)
{
    double matrix[n][n];

    if(argc == 1) {
        fracminhash_similarities(files, matrix);
    }
    else if(argc == 2) {
        double s = std::stod(argv[1]);
        fracminhash_similarities(files, matrix, s);
    }
    else {
        std::cout << "usage: optionally provide scaling factor\n";
        return -1;
    }

    print_matrix(matrix);
}
