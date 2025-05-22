#include <unordered_set>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include "hashfunctions.hpp"
#include "util.hpp"

const uint8_t k = 31;


struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};

static size_t set_intersection_size(const std::unordered_set<uint64_t> &set_a,
                                    const std::unordered_set<uint64_t> &set_b)
{
    if(set_b.size() < set_a.size()) {
        return set_intersection_size(set_b, set_a);
    }
    size_t cardinality = 0;
    for(std::unordered_set<uint64_t>::const_iterator it = set_a.begin(); it != set_a.end(); it++) {
        if(set_b.find(*it) != set_b.end())
            cardinality++;
    }
    return cardinality;
}

static size_t set_union_size(const std::unordered_set<uint64_t> &set_a,
                             const std::unordered_set<uint64_t> &set_b)
{
    std::unordered_set<uint64_t> set_union;
    set_union.insert(set_a.begin(), set_a.end());
    set_union.insert(set_b.begin(), set_b.end());
    return set_union.size();
}


void fill_ht(const std::filesystem::path &filepath,
             std::unordered_set<uint64_t> &kmerset, const double percentage)
{
    // TODO: implement naive jaccard here
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k});
    for(auto & record : stream) {
        for(auto && kmer : record.sequence() | kmer_view) {
            
        }
    }
}


double jaccard_similarity(const std::unordered_set<uint64_t> &kmerset_a,
                          const std::unordered_set<uint64_t> &kmerset_b)
{
    // TODO: implement naive jaccard here
    return 0.0;
}


void jaccard_similarities(const std::vector<std::filesystem::path> &filepaths, double matrix[n][n], const double percentage)
{
    // TODO: implement naive jaccard here
    for(int i = 0; i < n; i++) {
        for(int j = i+1; j < n; j++) {
            matrix[i][j] = matrix[j][i] = 0.0;
        }
    }
}


int main(int argc, char** argv)
{
    if(argc != 2) {
        std::cout << "usage: provide sample probability in [0,1]\n";
        return -1;
    }
    double percentage = std::stod(argv[1]);
    double matrix[n][n];

    jaccard_similarities(files, matrix, percentage);
    print_matrix(matrix);
}
