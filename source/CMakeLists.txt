add_executable(naivecounting0 naive_counting0.cpp)
target_link_libraries (naivecounting0 seqan3::seqan3)

add_executable(naivecounting1 naive_counting1.cpp)
target_link_libraries (naivecounting1 seqan3::seqan3)

add_executable(flajoletmartin0 flajolet_martin0.cpp)
target_link_libraries (flajoletmartin0 seqan3::seqan3)

add_executable(flajoletmartin1 flajolet_martin1.cpp)
target_link_libraries (flajoletmartin1 seqan3::seqan3)

add_executable(fmplus0 fmplus0.cpp)
target_link_libraries (fmplus0 seqan3::seqan3)

add_executable(fmplus1 fmplus1.cpp)
target_link_libraries (fmplus1 seqan3::seqan3)

add_executable(hyperloglog0 hyperloglog0.cpp)
target_link_libraries (hyperloglog0 seqan3::seqan3)

add_executable(hyperloglog1 hyperloglog1.cpp)
target_link_libraries (hyperloglog1 seqan3::seqan3)

add_executable(hyperloglog2 hyperloglog2.cpp)
target_link_libraries (hyperloglog2 seqan3::seqan3)

add_executable(jaccard0 jaccard0.cpp)
target_link_libraries (jaccard0 seqan3::seqan3)

add_executable(jaccard1 jaccard1.cpp)
target_link_libraries (jaccard1 seqan3::seqan3)

add_executable(jaccard_sample0 jaccard_sample0.cpp)
target_link_libraries (jaccard_sample0 seqan3::seqan3)

add_executable(jaccard_sample1 jaccard_sample1.cpp)
target_link_libraries (jaccard_sample1 seqan3::seqan3)

add_executable(minhash0 minhash0.cpp)
target_link_libraries (minhash0 seqan3::seqan3)

add_executable(minhash1 minhash1.cpp)
target_link_libraries (minhash1 seqan3::seqan3)

add_executable(fracminhash0 fracminhash0.cpp)
target_link_libraries (fracminhash0 seqan3::seqan3)

add_executable(fracminhash1 fracminhash1.cpp)
target_link_libraries (fracminhash1 seqan3::seqan3)

add_executable(fracminhash2 fracminhash2.cpp)
target_link_libraries (fracminhash2 seqan3::seqan3)