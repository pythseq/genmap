#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

#include "common.hpp"
#include "genmap_helper.hpp"

template <typename TSpec, typename TLengthSum, unsigned LEVELS, unsigned WORDS_PER_BLOCK>
unsigned GemMapFastFMIndexConfig<TSpec, TLengthSum, LEVELS, WORDS_PER_BLOCK>::SAMPLING = 10;

using namespace seqan;

int main(int argc, char const ** argv)
{
    // Argument parser
    ArgumentParser parser("classify");
    sharedSetup(parser);
    addDescription(parser, "TODO");

    //addOption(parser, ArgParseOption("E", "errors", "Number of errors", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("K", "kmers", "Path to the kmer-file", ArgParseArgument::INPUT_FILE, "IN"));

    std::vector<std::string> const fastaFileTypes {"fsa", "fna", "fastq", "fasta", "fas", "fa"};
    addOption(parser, ArgParseOption("F", "reads", "Path to the read-file", ArgParseArgument::INPUT_FILE, "IN"));
    setValidValues(parser, "reads", fastaFileTypes);

    //addOption(parser, ArgParseOption("T", "threads", "Number of threads", ArgParseArgument::INTEGER, "INT"));
    //setDefaultValue(parser, "threads", omp_get_max_threads());

    //addOption(parser, ArgParseOption("v", "verbose", "Outputs some additional information."));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // Retrieve input parameters
    //unsigned errors, threads;
    CharString kmers_path, reads_path;
    //getOptionValue(errors, parser, "errors");
    //getOptionValue(threads, parser, "threads");
    getOptionValue(kmers_path, parser, "kmers");
    getOptionValue(reads_path, parser, "reads");
    //bool verbose = isSet(parser, "verbose");

    // read kmers
    typedef StringSet<Dna5String, Owner<ConcatDirect<SizeSpec_<uint32_t/*TSeqNo*/, uint32_t/*TSeqPos*/> > > > TKmers;
    TKmers kmers;
    std::ifstream kmerFile(toCString(kmers_path));
    std::string line;
    while (getline(kmerFile, line, '\n'))
    {
        size_t const pos = line.find('\t');
        if (pos != std::string::npos)
        {
            Dna5String const kmer = line.substr(pos + 1);
            appendValue(kmers, kmer);
        }
    }

    // build index
    typedef Index<TKmers, FMIndex<Nothing, TGemMapFastFMIndexConfig<uint32_t/*BWTLength*/> > > TIndex;
    TIndex index(kmers);
    indexCreate(index, FibreSALF());

    // read reads and search
    uint32_t kmer_size = length(kmers[0]);
    std::set<uint32_t> kmers_found;
    SeqFileIn seqFileIn(toCString(reads_path));
    CharString id;
    Dna5String seq; // TODO: transform Iupac to Dna5
    while (!atEnd(seqFileIn))
    {
        readRecord(id, seq, seqFileIn);
        //std::cout << id << '\t' << seq << '\n';
        if (length(seq) < kmer_size)
            continue;

        for (uint32_t i = 0; i < length(seq) - kmer_size + 1; ++i)
        {
            Dna5String kmer = infixWithLength(seq, i, kmer_size);
            Dna5String kmer_rc = kmer;
            reverseComplement(kmer_rc); // TODO: inefficient
            if (kmer > kmer_rc)
                kmer = kmer_rc;

            // search kmer in index
            reverse(kmer);
            Iter<TIndex, VSTree<TopDown<> > > iter(index);
            if (goDown(iter, kmer))
            {
                for (auto occ : getOccurrences(iter))
                {
                    kmers_found.insert(occ.i1);
                }
            }

            //std::cout << kmer << '\n';
        }
    }

    // output to file
    std::string output_path = std::string(toCString(reads_path)) + ".result";
    std::ofstream output(output_path);
    output << "\n\n0\n";
    for (uint32_t const kmer_id : kmers_found)
    {
        output << (kmer_id + 1) << "\t1.0\n";
    }
    output.close();

    return 0;
}
