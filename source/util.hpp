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