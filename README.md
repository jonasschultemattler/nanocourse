# Genomic analysis using sketching techniques


## Setup

Requirements:

 - gcc >= 12 or clang >=17
 - cmake >= 3.20
 - git
 - python3

Checkout
```
git clone --recurse-submodules https://github.com/jonasschultemattler/cardinalities.git
```

See https://docs.seqan.de/seqan3/main_user/setup.html for compiler setup. Compile with

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release .. -D CMAKE_CXX_COMPILER=g++-14
make
```

Test
```
./build/source/naivecounting1 "data/ecoli1_k31_ust.fa.gz"
```

Create a virtual python environment
```
python3 -m venv venv
source venv/bin/activate
python3 -m pip install -r requirements.txt
```

Test
```
jupyter notebook course.ipynb
```


