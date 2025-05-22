#include <unordered_set>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>


const uint8_t k = 31;

struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};


uint64_t naive_couting(const std::filesystem::path &filepath)
{
    // TODO: implement naive counting here
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k});
    for(auto & record : stream) {
        for(auto && kmer : record.sequence() | kmer_view) {
        }
    }
    return 0;
}



int main(int argc, char** argv)
{
    const std::filesystem::path file = argv[1];
    std::cout << file << '\n';

    uint64_t count = naive_couting(file);
    std::cout << "Distinct kmers: " << count << '\n';

}
