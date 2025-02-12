# bgen2igd

Convert BGEN files to [IGD](https://github.com/aprilweilab/picovcf). Uses the
[limix/bgen](https://github.com/limix/bgen) library for BGEN reading. Converts
to hard calls by taking the allele with the maximum probability (ties are
broken by allele order in the file).

Mixed phasedness is not supported. Unphased data is restricted to haploid or
diploid. These are not fundamental restrictions, just limitations of the
current implementation.

## Building

```
git clone --recursive https://github.com/dcdehaas/bgen2igd.git
cd bgen2igd && mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j
```

Or use the provide Docker image:
```
docker build . -t bgen2igd:latest
```

## Using

```
bgen2igd <bgen filename> <output igd filename>
```

Or, via docker:
```
docker run -v $PWD:/pwd/ -it bgen2igd:latest bgen2igd /pwd/your_file.bgen /pwd/your_file.igd
```
