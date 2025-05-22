#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include "util.hpp"
#include "hyperloglog2.hpp"


static double jaccard(const double containment, const uint64_t size_a, const uint64_t size_b) {
    return containment*size_a/(size_a+size_b-containment*size_a);
}


void fracminsketch(const std::filesystem::path &filepath, std::vector<uint64_t> &frac_sketch,
                   const double s, uint64_t (*hashFunc)(uint64_t)=wyhash)
{
    // TODO: implement FracMinSketch here
    auto fin = seqan3::sequence_file_input<my_traits>{filepath};
    auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k});
    for(auto & record : fin) {
        for(auto && kmer : record.sequence() | kmer_view) {
            
        }
    }
}


double fracMinHash(const std::vector<uint64_t> &frac_sketch_a, const std::vector<uint64_t> &frac_sketch_b, const uint64_t size_a,
                   const double s)
{
    // TODO: implement FracMinHashing here
    return 0.0;
}


void fracminhash_similarities(const std::vector<std::filesystem::path> &filepaths, double matrix[n][n],
    const double s=0.1)
{
    // TODO: fill similarity matrix here
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
