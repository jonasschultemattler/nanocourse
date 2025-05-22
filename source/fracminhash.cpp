#include "hyperloglog.hpp"

const std::vector<std::filesystem::path> files = {
        "data/ecoli1_k31_ust.fa.gz",
        "data/ecoli2_k31_ust.fa.gz",
        "data/ecoli4_k31_ust.fa.gz",
        "data/salmonella_100_k31_ust.fa.gz"};
const int n = 4;


void print_matrix(double matrix[n][n]) {
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << '\n';
    }
}

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
    // TODO: implement FracMinHashing here
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
