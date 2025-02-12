#include <iostream>
#include <sstream>
#include <algorithm>

#include "bgen/bstring.h"
#include "picovcf/picovcf.hpp"

extern "C" {
#include <bgen/bgen.h>
}

using namespace picovcf;

// See https://www.chg.ox.ac.uk/~gav/bgen_format/spec/latest.html
static bool colexicographicPairCmp(const std::pair<size_t, size_t>& pair1,
                                   const std::pair<size_t, size_t>& pair2) {
    if (pair1.second == pair2.second) {
        return pair1.first < pair2.first;
    }
    return pair1.second < pair2.second;
}

// Convert a BGEN file to an IGD file. Supports phased and unphased, though there are
// some ploidy restrictions on the unphased data.
inline void bgenToIGD(const std::string& bgenFilename,
                      const std::string& outFilename,
                      std::string description = "",
                      bool verbose = false) {
    bgen_file* bgenFile = bgen_file_open(bgenFilename.c_str());
    PICOVCF_ASSERT_OR_MALFORMED(nullptr != bgenFile, "Invalid bgen file " << bgenFilename);
    const uint32_t numIndividuals = bgen_file_nsamples(bgenFile);
    PICOVCF_ASSERT_OR_MALFORMED(numIndividuals > 0, "No samples in BGEN file");

    const uint32_t numVariants = bgen_file_nvariants(bgenFile);
    PICOVCF_ASSERT_OR_MALFORMED(numVariants > 0, "No variants in BGEN file");

    std::stringstream ssMetaFile;
    ssMetaFile << bgenFilename << ".meta";
    std::string metaFilename = ssMetaFile.str();
    bgen_metafile* metaFile = bgen_metafile_create(bgenFile, metaFilename.c_str(), 1, 0);
    PICOVCF_ASSERT_OR_MALFORMED(nullptr != metaFile, "Failed to create .meta file " << metaFilename);

    const bgen_partition* partition = bgen_metafile_read_partition(metaFile, 0);
    PICOVCF_ASSERT_OR_MALFORMED(nullptr != partition, "Failed to load partition 0");
    bgen_metafile_close(metaFile);
    metaFile = nullptr;
    PICOVCF_RELEASE_ASSERT(bgen_partition_nvariants(partition) == numVariants);

    const bgen_variant* const firstVariant = bgen_partition_get_variant(partition, 0);
    bgen_genotype* genotype = bgen_file_open_genotype(bgenFile, firstVariant->genotype_offset);
    PICOVCF_RELEASE_ASSERT(genotype != nullptr);
    const uint8_t minPloidy = bgen_genotype_min_ploidy(genotype);
    const uint8_t maxPloidy = bgen_genotype_max_ploidy(genotype);
    PICOVCF_ASSERT_OR_MALFORMED(minPloidy == maxPloidy, "bgen2igd only support single-ploidy files, currently");
    const size_t ploidy = static_cast<size_t>(maxPloidy);
    const bool isPhased = bgen_genotype_phased(genotype);

    std::vector<std::string> variantIds;
    std::ofstream outFile(outFilename, std::ios::binary);
    IGDWriter writer(ploidy, numIndividuals, isPhased);
    writer.writeHeader(outFile, bgenFilename, description);

    for (size_t i = 0; i < bgen_partition_nvariants(partition); i++) {
        const bgen_variant* const variant = bgen_partition_get_variant(partition, i);
        PICOVCF_RELEASE_ASSERT(variant != nullptr);

        bgen_genotype* genotype = bgen_file_open_genotype(bgenFile, variant->genotype_offset);
        PICOVCF_RELEASE_ASSERT(genotype != nullptr);
        {
            const uint8_t _minPloidy = bgen_genotype_min_ploidy(genotype);
            const uint8_t _maxPloidy = bgen_genotype_max_ploidy(genotype);
            PICOVCF_ASSERT_OR_MALFORMED(_minPloidy == _maxPloidy, "bgen2igd only support single-ploidy files, currently");
            PICOVCF_ASSERT_OR_MALFORMED(_minPloidy == ploidy, "bgen2igd only support single-ploidy files, currently");
        }
        PICOVCF_ASSERT_OR_MALFORMED(isPhased == bgen_genotype_phased(genotype), "bgen2igd only support same phasedness throughout file, currently");

        const size_t probsPerIndividual = bgen_genotype_ncombs(genotype);
        const size_t numProbabilities = probsPerIndividual * numIndividuals;
        PICOVCF_RELEASE_ASSERT(numProbabilities >= numIndividuals);
        std::vector<double> probabilities(numProbabilities);
        PICOVCF_ASSERT_OR_MALFORMED(0 == bgen_genotype_read(genotype, probabilities.data()),
                                    "Failed to read genotype for variant from BGEN");
        const size_t numAlleles = static_cast<size_t>(bgen_genotype_nalleles(genotype));
        PICOVCF_ASSERT_OR_MALFORMED(0 != numAlleles, "No alleles for BGEN genotype");

        IGDSampleList missingData;
        const size_t numSampleLists = isPhased ? numAlleles : (numAlleles * ploidy);
        std::vector<IGDSampleList> variantGtData(numSampleLists);

        // This restriction is just because enumerating the order of genotype probabilities gets
        // annoying for ploidy > 2, and we haven't had need for such data (yet).
        PICOVCF_ASSERT_OR_MALFORMED(isPhased || ploidy <= 2, "bgen2igd only support haploid or diploid data, for unphased inputs");
        SampleT sampleIndex = 0;
        if (isPhased) {
            const size_t numSamples = numIndividuals * ploidy;
            for (size_t sampleId = 0; sampleId < numSamples; sampleId++) {
                const size_t baseOffset = (numAlleles * sampleId);
                double maxProb = 0.0;
                size_t maxAllele = 0;
                double sumProb = 0.0;

                if (bgen_genotype_missing(genotype, sampleId / ploidy)) {
                    missingData.push_back(sampleId);
                } else {
                    for (size_t k = 0; k < numAlleles; k++) {
                        const size_t idx = baseOffset + k;
                        const double prob = probabilities.at(idx);
                        sumProb += prob;
                        if (prob > maxProb) {
                            maxProb = prob;
                            maxAllele = k;
                        }
                    }
                    // The BGEN spec leaves out the k'th allele, but the library we use recovers it, so
                    // this should be true.
                    PICOVCF_RELEASE_ASSERT(sumProb >= 0.99999 && sumProb <= 1.00001);

                    // We assume that the allele at position 0 is the reference allele.
                    if (maxAllele != 0) {
                        variantGtData.at(maxAllele).push_back(sampleId);
                    }
                }

            }
        } else {
            // There has to be an easier way to do this. We could do this more generically than
            // for just haploid/diploid, but it's not needed right now.
            std::vector<std::pair<size_t, size_t>> order;
            for (size_t k = 0; k < numAlleles; k++) {
                if (ploidy == 2) {
                    for (size_t l = k; l < numAlleles; l++) {
                        order.emplace_back(k, l);
                    }
                } else {
                    PICOVCF_RELEASE_ASSERT(ploidy == 1);
                    order.emplace_back(k, 0);
                }
            }
            std::sort(order.begin(), order.end(), colexicographicPairCmp);
            PICOVCF_RELEASE_ASSERT(order.size() == probsPerIndividual);

            for (size_t indivId = 0; indivId < numIndividuals; indivId++) {
                const size_t baseOffset = (probsPerIndividual * indivId);

                if (bgen_genotype_missing(genotype, indivId)) {
                    missingData.push_back(indivId);
                } else {
                    double maxProb = 0.0;
                    ssize_t maxAllele1 = 0;
                    ssize_t maxAllele2 = 0;
                    double sumProb = 0.0;

                    for (size_t k = 0; k < order.size(); k++) {
                        const auto& hapPair = order[k];
                        const size_t idx = baseOffset + k;
                        const double prob = probabilities.at(idx);
                        sumProb += prob;
                        if (prob > maxProb) {
                            maxProb = prob;
                            maxAllele1 = hapPair.first;
                            maxAllele2 = hapPair.second;
                        }
                    }
                    // The BGEN spec leaves out the k'th allele, but the library we use recovers it, so
                    // this should be true.
                    PICOVCF_RELEASE_ASSERT(sumProb >= 0.99999 && sumProb <= 1.00001);

                    // We assume that the allele at position 0 is the reference allele.
                    if (ploidy == 1) {
                        if (maxAllele1 != 0) {
                            variantGtData.at(maxAllele1).push_back(indivId);
                        }
                    } else {
                        if (maxAllele1 == maxAllele2) {
                            if (maxAllele1 != 0) {
                                variantGtData.at(maxAllele1 + numAlleles).push_back(indivId);
                            }
                        } else {
                            if (maxAllele1 != 0) {
                                variantGtData.at(maxAllele1).push_back(indivId);
                            }
                            if (maxAllele2 != 0) {
                                variantGtData.at(maxAllele2).push_back(indivId);
                            }
                        }
                    }
                }
            }
        }

        std::string refAllele(bgen_string_data(variant->allele_ids[0]),
                              bgen_string_length(variant->allele_ids[0]));
        const auto position = variant->position;
        for (size_t i = 0; i < numAlleles; i++) {
            if (i == 0) {
                PICOVCF_RELEASE_ASSERT(variantGtData[i].empty());
                continue;
            }
            std::string alt(bgen_string_data(variant->allele_ids[i]),
                            bgen_string_length(variant->allele_ids[i]));
            const uint8_t numCopies = isPhased ? 0 : 1;
            writer.writeVariantSamples(outFile, position, refAllele, alt, variantGtData[i], false, numCopies);
            if (!isPhased && ploidy > 1) {
                const size_t copies2Index = i + numAlleles;
                const auto& sampleList = variantGtData[copies2Index];
                if (!sampleList.empty()) {
                    writer.writeVariantSamples(outFile,
                                               position,
                                               refAllele,
                                               alt,
                                               sampleList,
                                               false,
                                               /*numCopies=*/2);
                }
            }
        }
        if (!missingData.empty()) {
            writer.writeVariantSamples(outFile, position, refAllele, "", missingData, true);
        }

        bgen_genotype_close(genotype);
        genotype = nullptr;
    }

    writer.writeIndex(outFile);
    writer.writeVariantInfo(outFile);
    outFile.seekp(0);
    writer.writeHeader(outFile, bgenFilename, description);

    if (verbose) {
        std::cout << "Wrote " << writer.m_totalCount << " total variants" << std::endl;
        std::cout << "Of which " << writer.m_sparseCount << " were written sparsely" << std::endl;
    }
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: bgen2igd <bgen file> <output IGD file>" << std::endl;
        return 1;
    }
    std::string inFile = argv[1];
    std::string outFile = argv[2];
    bgenToIGD(inFile, outFile, "", true);
    return 0;
}
