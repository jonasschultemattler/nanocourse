#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

#include "hashfunctions.hpp"


const std::vector<std::filesystem::path> files = {
        "data/ecoli1_k31_ust.fa.gz",
        "data/ecoli2_k31_ust.fa.gz",
        "data/ecoli4_k31_ust.fa.gz",
        "data/salmonella_100_k31_ust.fa.gz"};
const int n = 4;
const uint8_t k = 31;

const uint64_t prime = (1UL << 61) - 1;
const int seed = 1;


struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};

void print_matrix(double matrix[n][n]) {
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << '\n';
    }
}

void minhash_sketch(const std::filesystem::path &filepath, uint64_t permutations, uint64_t minhashs[],
                    const uint64_t a[], const uint64_t b[], uint64_t (*hashFunc)(uint64_t))
{
    // TODO: implement MinHashSketch here
    for(int i = 0; i < permutations; i++)
        minhashs[i] = UINT64_MAX;

    auto fin = seqan3::sequence_file_input<my_traits>{filepath};
    auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k});
    for(auto & record : fin) {
        for(auto && kmer : record.sequence() | kmer_view) {
            uint64_t hash = hashFunc(kmer);
            for(int i = 0; i < permutations; i++) {
                minhashs[i] = std::min((hash*a[i]+b[i]) % prime, minhashs[i]);
            }
        }
    }
}


double minhash_similarity(const uint64_t minhashs_a[], const uint64_t minhashs_b[], const int permutations)
{
    // TODO: implement MinHashing here
    int y = 0;
    for(int i = 0; i < permutations; i++)
        y += minhashs_a[i] == minhashs_b[i];

    return (double) y/permutations;
}


void minhash_similarities(const std::vector<std::filesystem::path> &filepaths, double matrix[n][n],
    const int permutations=100, uint64_t (*hashFunc)(uint64_t)=wyhash) 
{
    uint64_t a[permutations];
    uint64_t b[permutations];
    for(int i = 0; i < permutations; i++) {
        a[i] = (1+std::rand()) % prime;
        b[i] = std::rand() % prime;
    }

    for(int i = 0; i < n; i++) {
        matrix[i][i] = 1;
        uint64_t minhashs_i[permutations];
        minhash_sketch(filepaths[i], permutations, minhashs_i, a, b, hashFunc);
        for(int j = i+1; j < n; j++) {
            uint64_t minhashs_j[permutations];
            minhash_sketch(filepaths[j], permutations, minhashs_j, a, b, hashFunc);
            matrix[i][j] = matrix[j][i] = minhash_similarity(minhashs_i, minhashs_j, permutations);
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
