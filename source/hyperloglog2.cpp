#include "hyperloglog2.hpp"

// TODO: implement HyperLogLog in hperloglog2.hpp

int main(int argc, char** argv)
{
    const std::filesystem::path file = argv[1];
    std::cout << file << '\n';

    uint64_t count;
    if(argc == 3) {
        const uint8_t precision = atoi(argv[2]);
        count = hyperloglog(file, precision);
    }
    else {
        count = hyperloglog(file);
    }
    
    std::cout << "Distinct kmers: " << count << '\n';

}
