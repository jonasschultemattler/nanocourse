#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <random>
#include "hashfunctions.hpp"

const uint8_t k = 31;

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


static inline constexpr uint8_t leading_zeros(const uint64_t x) {
    return std::countr_zero(x);
}


uint64_t fmplus(const std::filesystem::path &filepath)
{
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k});

    const uint64_t s = 100;

    uint8_t registers[s];
    for(uint64_t i=0; i < s; i++)
        registers[i] = 0;

    for(auto & record : stream) {
        for(auto && kmer : record.sequence() | kmer_view) {
            auto hashes = multi_hash(kmer, s);
            for(uint64_t i=0; i < s; i++) {
                const uint8_t zeros = leading_zeros(hashes[i]);
                registers[i] = std::max(registers[i], zeros);
            }
        }
    }

    double sum = 0.0;
    for(uint64_t j=0; j < s; ++j) {
        sum += 1.0 / (1ULL << registers[j]);
    }
    double phi = 0.77351;
    return s/(sum*phi);

}


int main(int argc, char** argv)
{
    const std::filesystem::path file = argv[1];
    std::cout << file << '\n';
    
    uint64_t count = fmplus(file);

    std::cout << "Distinct kmers: " << count << '\n';
}
