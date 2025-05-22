add_executable(naivecounting naive_counting.cpp)
target_link_libraries (naivecounting seqan3::seqan3)

add_executable(flajoletmartin flajolet_martin.cpp)
target_link_libraries (flajoletmartin seqan3::seqan3)

add_executable(hyperloglog hyperloglog.cpp)
target_link_libraries (hyperloglog seqan3::seqan3)

add_executable(jaccard jaccard.cpp)
target_link_libraries (jaccard seqan3::seqan3)

add_executable(jaccard_sample jaccard_sample.cpp)
target_link_libraries (jaccard_sample seqan3::seqan3)

add_executable(minhash minhash.cpp)
target_link_libraries (minhash seqan3::seqan3)

add_executable(fracminhash fracminhash.cpp)
target_link_libraries (fracminhash seqan3::seqan3)