#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include "hashfunctions.hpp"

const uint8_t k = 31;

struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};


static inline constexpr uint8_t leading_zeros(const uint64_t x) {
    return std::countl_zero(x);
}


uint64_t flajolet_martin(const std::filesystem::path &filepath, uint64_t (*hashFunc)(uint64_t)=wyhash)
{
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k});

    uint8_t l = 0;

    for(auto & record : stream) {
        for(auto && kmer : record.sequence() | kmer_view) {
            uint64_t hash = hashFunc(kmer);
            uint8_t zeros = leading_zeros(hash);
            l = std::max(l, zeros);
        }
    }

    return 1 << l;
}


int main(int argc, char** argv)
{
    const std::filesystem::path file = argv[1];
    std::cout << file << '\n';
    
    uint64_t count = flajolet_martin(file);

    std::cout << "Distinct kmers: " << count << '\n';
}
