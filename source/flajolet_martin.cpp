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
    // TODO: implement Flajolet-Martin’s algorithm here
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


// uint64_t flajolet_martin(const std::filesystem::path &filepath, std::vector<uint64_t (*)(uint64_t)> hashFuncs)
// {
//     // TODO: implement Flajolet-Martin’s algorithm here
//     auto stream = seqan3::sequence_file_input<my_traits>{filepath};
//     auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k});

//     std::vector<uint8_t> l;
//     for(int i = 0; i < hashFuncs.size(); i++)
//         l.push_back(0);

//     for(auto & record : stream) {
//         for(auto && kmer : record.sequence() | kmer_view) {
//             for(int i = 0; i < hashFuncs.size(); i++) {
//                 uint64_t hash = hashFuncs[i](kmer);
//                 uint8_t zeros = leading_zeros(hash);
//                 l[i] = std::max(l[i], zeros);
//             }
//         }
//     }
//     uint8_t l_mean = 0;
//     for(int i = 0; i < l.size(); i++)
//         l_mean += l[i];
//     l_mean /= l.size();

//     return 1 << l_mean;
// }


int main(int argc, char** argv)
{
    const std::filesystem::path file = argv[1];
    std::cout << file << '\n';
    // uint64_t count;
    // if(argc == 3) {
    //     if(argv[2] == "stdhash")
    //         count = flajolet_martin(file, stdhash);
    //     if(argv[2] == "wyhash")
    //         count = flajolet_martin(file, wyhash);
    // }
    // else {
    //     count = flajolet_martin(file);
    // }
    uint64_t count = flajolet_martin(file);

    std::cout << "Distinct kmers: " << count << '\n';

}
