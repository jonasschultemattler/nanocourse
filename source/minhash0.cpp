#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include "hashfunctions.hpp"
#include "util.hpp"

const uint8_t k = 31;

const uint64_t prime = (1UL << 61) - 1;
const int seed = 1;


struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};


std::vector<uint64_t> multi_hash(uint64_t input, size_t N) {
    std::vector<uint64_t> hashes;
    for (size_t i = 0; i < N; ++i) {
        // Combine input with index as salt
        uint64_t combined = input ^ (0x9e3779b97f4a7c15ULL * i);
        hashes.push_back(util::Fingerprint(combined));
    }
    return hashes;
}


void minhash_sketch(const std::filesystem::path &filepath, uint64_t permutations, uint64_t minhashs[],
                    const uint64_t a[], const uint64_t b[], uint64_t (*hashFunc)(uint64_t))
{
    // TODO: implement MinHashSketch here
    auto fin = seqan3::sequence_file_input<my_traits>{filepath};
    auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k});
    for(auto & record : fin) {
        for(auto && kmer : record.sequence() | kmer_view) {
            
        }
    }
}


double minhash_similarity(const uint64_t minhashs_a[], const uint64_t minhashs_b[], const int permutations)
{
    // TODO: implement MinHashing here
    return 0.0;
}


void minhash_similarities(const std::vector<std::filesystem::path> &filepaths, double matrix[n][n],
    const int permutations=100, uint64_t (*hashFunc)(uint64_t)=wyhash) 
{
    // TODO: implement MinHashing here

    for(int i = 0; i < n; i++) {
        for(int j = i+1; j < n; j++) {
            matrix[i][j] = matrix[j][i] = 0.0;
        }
    }
}


int main(int argc, char** argv)
{
    double matrix[n][n];

    if(argc == 1) {
        minhash_similarities(files, matrix);
    }
    else if(argc == 2) {
        const int permutations = std::stoi(argv[1]);
        minhash_similarities(files, matrix, permutations);
    }
    else if(argc == 3) {
        const int permutations = std::stoi(argv[1]);
        if(std::string(argv[2]) == "stdhash")
            minhash_similarities(files, matrix, permutations, stdhash);
        else if(std::string(argv[2]) == "wyhash")
            minhash_similarities(files, matrix, permutations, wyhash);
        else if(std::string(argv[2]) == "farmhash")
            minhash_similarities(files, matrix, permutations, farmhash);
        else
            std::cout << "no such hash function\n";
    }
    else {
        std::cout << "usage: optionally provide number of permutations and hash function\n";
        return -1;
    }
    
    print_matrix(matrix);

}
