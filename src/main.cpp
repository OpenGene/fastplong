#include <stdio.h>
#include "fastqreader.h"
#include <time.h>
#include "cmdline.h"
#include <sstream>
#include "util.h"
#include "options.h"
#include "processor.h"
#include "evaluator.h"
#include "sequence.h"

// TODO: code refactoring to remove these global variables
string command;
mutex logmtx;

int main(int argc, char* argv[]){
    // display version info if no argument is given
    if(argc == 1) {
        cerr << "fastplong: ultra-fast FASTQ preprocessing and quality control for long reads" << endl << "version " << FASTPLONG_VER << endl;
        //cerr << "fastplong --help to see the help"<<endl;
        //return 0;
    }
    if (argc == 2 && (strcmp(argv[1], "-v")==0 || strcmp(argv[1], "--version")==0)){
        cout << "fastplong " << FASTPLONG_VER << endl;
        return 0;
    }
    cmdline::parser cmd;
    // input/output
    cmd.add<string>("in", 'i', "read input file name", false, "");
    cmd.add<string>("out", 'o', "read output file name", false, "");
    cmd.add<string>("failed_out", 0, "specify the file to store reads that cannot pass the filters.", false, "");
    cmd.add<int>("compression", 'z', "compression level for gzip output (1 ~ 9). 1 is fastest, 9 is smallest, default is 4.", false, 4);
    cmd.add("stdin", 0, "input from STDIN.");
    cmd.add("stdout", 0, "stream passing-filters reads to STDOUT. This option will result in interleaved FASTQ output for paired-end output. Disabled by default.");
    cmd.add<int>("reads_to_process", 0, "specify how many reads/pairs to be processed. Default 0 means process all reads.", false, 0);
    cmd.add("dont_overwrite", 0, "don't overwrite existing files. Overwritting is allowed by default.");
    cmd.add("verbose", 'V', "output verbose log information (i.e. when every 1M reads are processed).");

    // adapter
    cmd.add("disable_adapter_trimming", 'A', "adapter trimming is enabled by default. If this option is specified, adapter trimming is disabled");
    cmd.add<string>("start_adapter", 's', "the adapter sequence at read start (5').", false, "auto");
    cmd.add<string>("end_adapter", 'e', "the adapter sequence at read end (3').", false, "auto");
    cmd.add<string>("adapter_fasta", 'a', "specify a FASTA file to trim both read by all the sequences in this FASTA file", false, "");
    cmd.add<double>("distance_threshold", 'd', "threshold of sequence-adapter-distance/adapter-length (0.0 ~ 1.0), greater value means more adapters detected", false, 0.25);
    cmd.add<int>("trimming_extension", 0, "when an adapter is detected, extend the trimming to make cleaner trimming, default 10 means trimming 10 bases more", false, 10);

    // trimming
    cmd.add<int>("trim_front", 'f', "trimming how many bases in front for read, default is 0", false, 0);
    cmd.add<int>("trim_tail", 't', "trimming how many bases in tail for read, default is 0", false, 0);
    
    // polyX tail trimming
    cmd.add("trim_poly_x", 'x', "enable polyX trimming in 3' ends.");
    cmd.add<int>("poly_x_min_len", 0, "the minimum length to detect polyX in the read tail. 10 by default.", false, 10);

    // cutting by quality at front or tail
    cmd.add("cut_front", '5', "move a sliding window from front (5') to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.");
    cmd.add("cut_tail", '3', "move a sliding window from tail (3') to front, drop the bases in the window if its mean quality < threshold, stop otherwise.");
    cmd.add<int>("cut_window_size", 'W', "the window size option shared by cut_front, cut_tail. Range: 1~1000, default: 4", false, 4);
    cmd.add<int>("cut_mean_quality", 'M', "the mean quality requirement option shared by cut_front, cut_tail. Range: 1~36 default: 20 (Q20)", false, 20);
    cmd.add<int>("cut_front_window_size", 0, "the window size option of cut_front, default to cut_window_size if not specified", false, 4);
    cmd.add<int>("cut_front_mean_quality", 0, "the mean quality requirement option for cut_front, default to cut_mean_quality if not specified", false, 20);
    cmd.add<int>("cut_tail_window_size", 0, "the window size option of cut_tail, default to cut_window_size if not specified", false, 4);
    cmd.add<int>("cut_tail_mean_quality", 0, "the mean quality requirement option for cut_tail, default to cut_mean_quality if not specified", false, 20);

    // mask low quality regions with N
    cmd.add("mask", 'N', "mask the low quality regions with N, these regions are detected by sliding window with mean quality < mask_mean_quality.");
    cmd.add<int>("mask_window_size", 0, "the size of the sliding window to evaluate the mean quality for N masking(5~1000000), default: 50", false, 50);
    cmd.add<int>("mask_mean_quality", 0, "the mean quality requirement for sliding window N masking (5~30), default: 10 (Q10)", false, 10);

    // break reads into high-quality fragments, and discard low-quality fragments
    cmd.add("break", 'b', "break the reads by discarding the low quality regions, these regions are detected by sliding window with mean quality < break_mean_quality.");
    cmd.add<int>("break_window_size", 0, "the size of the sliding window to evaluate the mean quality for sliding window breaking(5~1000000), default: 100", false, 100);
    cmd.add<int>("break_mean_quality", 0, "the mean quality requirement for sliding window breaking (5~30), default: 10 (Q10)", false, 10);

    // quality filtering
    cmd.add("disable_quality_filtering", 'Q', "quality filtering is enabled by default. If this option is specified, quality filtering is disabled");
    cmd.add<int>("qualified_quality_phred", 'q', "the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified.", false, 15);
    cmd.add<int>("unqualified_percent_limit", 'u', "how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40%", false, 40);
    cmd.add<int>("n_base_limit", 0, "if number of N base is >n_base_limit, then this read is discarded (0~1000000). 0 means no N allowed, default 1000000 means no N limit", false, 1000000);
    cmd.add<int>("n_percent_limit", 'n', "if one read's N base percentage is >n_percent_limit, then this read is discarded (0~100). Default 10 means 10%", false, 10);
    cmd.add<int>("mean_qual", 'm', "if one read's mean_qual quality score <mean_qual, then this read is discarded. Default 0 means no requirement", false, 0);

    // length filtering
    cmd.add("disable_length_filtering", 'L', "length filtering is enabled by default. If this option is specified, length filtering is disabled");
    cmd.add<int>("length_required", 'l', "reads shorter than length_required will be discarded, default is 20.", false, 20);
    cmd.add<int>("length_limit", 0, "reads longer than length_limit will be discarded, default 0 means no limitation.", false, 0);

    // low complexity filtering
    cmd.add("low_complexity_filter", 'y', "enable low complexity filter. The complexity is defined as the percentage of base that is different from its next base (base[i] != base[i+1]).");
    cmd.add<int>("complexity_threshold", 'Y', "the threshold for low complexity filter (0~100). Default is 30, which means 30% complexity is required.", false, 30);

    // reporting
    cmd.add<string>("json", 'j', "the json format report file name", false, "fastplong.json");
    cmd.add<string>("html", 'h', "the html format report file name", false, "fastplong.html");
    cmd.add<string>("report_title", 'R', "should be quoted with \' or \", default is \"fastplong report\"", false, "fastplong report");

    // threading
    cmd.add<int>("thread", 'w', "worker thread number, default is 3", false, 3);

    // split the output
    cmd.add<int>("split", 0, "split output by limiting total split file number with this option (2~999), a sequential number prefix will be added to output name ( 0001.out.fq, 0002.out.fq...), disabled by default", false, 0);
    cmd.add<long>("split_by_lines", 0, "split output by limiting lines of each file with this option(>=1000), a sequential number prefix will be added to output name ( 0001.out.fq, 0002.out.fq...), disabled by default", false, 0);
    cmd.add<int>("split_prefix_digits", 0, "the digits for the sequential number padding (1~10), default is 4, so the filename will be padded as 0001.xxx, 0 to disable padding", false, 4);
    
    cmd.parse_check(argc, argv);

    if(argc == 1) {
        cerr << cmd.usage() <<endl;
    }

    if(argc == 1) {
        return 0;
    }

    Options opt;

    // I/O
    opt.in = cmd.get<string>("in");
    opt.out = cmd.get<string>("out");
    opt.failedOut = cmd.get<string>("failed_out");
    opt.compression = cmd.get<int>("compression");
    opt.readsToProcess = cmd.get<int>("reads_to_process");
    opt.dontOverwrite = cmd.exist("dont_overwrite");
    opt.inputFromSTDIN = cmd.exist("stdin");
    opt.outputToSTDOUT = cmd.exist("stdout");
    opt.verbose = cmd.exist("verbose");


    // adapter cutting
    opt.adapter.enabled = !cmd.exist("disable_adapter_trimming");
    opt.adapter.sequenceStart = cmd.get<string>("start_adapter");
    opt.adapter.sequenceEnd = cmd.get<string>("end_adapter");
    opt.adapter.fastaFile = cmd.get<string>("adapter_fasta");
    opt.adapter.edMax = cmd.get<double>("distance_threshold");
    opt.adapter.trimmingExtension = cmd.get<int>("trimming_extension");

    // if the start adapter is specified and the end is not, use the reverse complement of start adapter
    if(opt.adapter.sequenceStart != "auto" && opt.adapter.sequenceEnd == "auto") {
        opt.adapter.sequenceEnd = Sequence::reverseComplement((&opt.adapter.sequenceStart));
    }

    if(!opt.adapter.fastaFile.empty()) {
        opt.loadFastaAdapters();
    }

    // trimming
    opt.trim.front = cmd.get<int>("trim_front");
    opt.trim.tail = cmd.get<int>("trim_tail");

    // polyX tail trimming
    if(cmd.exist("trim_poly_x")) {
        opt.polyXTrim.enabled = true;
    }
    opt.polyXTrim.minLen = cmd.get<int>("poly_x_min_len");


    // sliding window cutting by quality
    opt.qualityCut.enabledFront = cmd.exist("cut_front");
    opt.qualityCut.enabledTail = cmd.exist("cut_tail");

    opt.qualityCut.windowSizeShared = cmd.get<int>("cut_window_size");
    opt.qualityCut.qualityShared = cmd.get<int>("cut_mean_quality");

    if(cmd.exist("cut_front_window_size"))
        opt.qualityCut.windowSizeFront = cmd.get<int>("cut_front_window_size");
    else
        opt.qualityCut.windowSizeFront = opt.qualityCut.windowSizeShared;
    if(cmd.exist("cut_front_mean_quality"))
        opt.qualityCut.qualityFront = cmd.get<int>("cut_front_mean_quality");
    else
        opt.qualityCut.qualityFront = opt.qualityCut.qualityShared;

    if(cmd.exist("cut_tail_window_size"))
        opt.qualityCut.windowSizeTail = cmd.get<int>("cut_tail_window_size");
    else
        opt.qualityCut.windowSizeTail = opt.qualityCut.windowSizeShared;
    if(cmd.exist("cut_tail_mean_quality"))
        opt.qualityCut.qualityTail = cmd.get<int>("cut_tail_mean_quality");
    else
        opt.qualityCut.qualityTail = opt.qualityCut.qualityShared;

    // raise a warning if cutting option is not enabled but -W/-M is enabled
    if(!opt.qualityCut.enabledFront && !opt.qualityCut.enabledTail) {
        if(cmd.exist("cut_window_size") || cmd.exist("cut_mean_quality") 
            || cmd.exist("cut_front_window_size") || cmd.exist("cut_front_mean_quality") 
            || cmd.exist("cut_tail_window_size") || cmd.exist("cut_tail_mean_quality") )
            cerr << "WARNING: you specified the options for cutting by quality, but forgot to enable any of cut_front/cut_tail/cut_right. This will have no effect." << endl;
    }

    // quality filtering
    opt.qualfilter.enabled = !cmd.exist("disable_quality_filtering");
    opt.qualfilter.qualifiedQual = num2qual(cmd.get<int>("qualified_quality_phred"));
    opt.qualfilter.unqualifiedPercentLimit = cmd.get<int>("unqualified_percent_limit");
    opt.qualfilter.avgQualReq = cmd.get<int>("mean_qual");
    opt.qualfilter.nBasePercentLimit = cmd.get<int>("n_percent_limit");
    opt.qualfilter.nBaseLimit = cmd.get<int>("n_base_limit");

    // length filtering
    opt.lengthFilter.enabled = !cmd.exist("disable_length_filtering");
    opt.lengthFilter.requiredLength = cmd.get<int>("length_required");
    opt.lengthFilter.maxLength = cmd.get<int>("length_limit");

    // low complexity filter
    opt.complexityFilter.enabled = cmd.exist("low_complexity_filter");
    opt.complexityFilter.threshold = (min(100, max(0, cmd.get<int>("complexity_threshold")))) / 100.0;

    // N masking by quality
    opt.mask.enabled = cmd.exist("mask");
    opt.mask.windowSize = cmd.get<int>("mask_window_size");
    opt.mask.quality = cmd.get<int>("mask_mean_quality");

    // break reads into high-quality fragments, and discard low-quality fragments
    opt.breakOpt.enabled = cmd.exist("break");
    opt.breakOpt.windowSize = cmd.get<int>("break_window_size");
    opt.breakOpt.quality = cmd.get<int>("break_mean_quality");

    // threading
    opt.thread = cmd.get<int>("thread");

    // reporting
    opt.jsonFile = cmd.get<string>("json");
    opt.htmlFile = cmd.get<string>("html");
    opt.reportTitle = cmd.get<string>("report_title");

    // splitting
    opt.split.enabled = cmd.exist("split") || cmd.exist("split_by_lines");
    opt.split.digits = cmd.get<int>("split_prefix_digits");
    if(cmd.exist("split") && cmd.exist("split_by_lines")) {
        error_exit("You cannot set both splitting by file number (--split) and splitting by file lines (--split_by_lines), please choose either.");
    }
    if(cmd.exist("split")) {
        opt.split.number = cmd.get<int>("split");
        opt.split.needEvaluation = true;
        opt.split.byFileNumber = true;
    }
    if(cmd.exist("split_by_lines")) {
        long lines = cmd.get<long>("split_by_lines");
        if(lines % 4 != 0) {
            error_exit("Line number (--split_by_lines) should be a multiple of 4");
        }
        opt.split.size = lines / 4; // 4 lines per record
        opt.split.needEvaluation = false;
        opt.split.byFileLines = true;
    }

    if(opt.inputFromSTDIN || opt.in=="/dev/stdin") {
        if(opt.split.needEvaluation) {
            error_exit("Splitting by file number is not supported in STDIN mode");
        }
    }

    stringstream ss;
    for(int i=0;i<argc;i++){
        ss << argv[i] << " ";
    }
    command = ss.str();

    time_t t1 = time(NULL);

    bool supportEvaluation = !opt.inputFromSTDIN && opt.in!="/dev/stdin";

    Evaluator eva(&opt);
    if(supportEvaluation) {
        eva.evaluateSeqLenAndCheckRNA();
    }

    long readNum = 0;

    // using evaluator to guess how many reads in total
    if(opt.shallDetectAdapter()) {
        if(!supportEvaluation)
            cerr << "Adapter auto-detection is disabled for STDIN mode" << endl;
        else {
            eva.evalAdapterAndReadNum(&opt, readNum);
            cerr << endl;
        }
    }

    opt.validate();

    // using evaluator to guess how many reads in total
    if(opt.split.needEvaluation && supportEvaluation) {
        // if readNum is not 0, means it is already evaluated by other functions
        if(readNum == 0) {
            eva.evaluateReadNum(readNum);
        }
        opt.split.size = readNum / opt.split.number;
        // one record per file at least
        if(opt.split.size <= 0) {
            opt.split.size = 1;
            cerr << "WARNING: the input file has less reads than the number of files to split" << endl;
        }
    }

    Processor p(&opt);
    p.process();
    
    time_t t2 = time(NULL);

    cerr << endl << "JSON report: " << opt.jsonFile << endl;
    cerr << "HTML report: " << opt.htmlFile << endl;
    cerr << endl << command << endl;
    cerr << "fastplong v" << FASTPLONG_VER << ", time used: " << (t2)-t1 << " seconds" << endl;

    return 0;
}
