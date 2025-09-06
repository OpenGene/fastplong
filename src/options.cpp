#include "options.h"
#include "util.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include "fastareader.h"

Options::Options(){
    in = "";
    out = "";
    reportTitle = "fastplong report";
    thread = 3;
    compression = 4;
    dontOverwrite = false;
    inputFromSTDIN = false;
    outputToSTDOUT = false;
    readsToProcess = 0;
    verbose = false;
    writerBufferSize = 0x01L<<22; // 4M writer buffer for per output by default
    isRNA = false;
}

void Options::init() {
}


bool Options::adapterCuttingEnabled() {
    if(adapter.enabled){
        if(!adapter.sequenceStart.empty() || !adapter.sequenceEnd.empty())
            return true;
    }
    return false;
}

bool Options::polyXTrimmingEnabled() {
    return polyXTrim.enabled;
}

void Options::loadFastaAdapters() {
    if(adapter.fastaFile.empty()) {
        adapter.hasFasta = false;
        return;
    }

    check_file_valid(adapter.fastaFile);

    FastaReader reader(adapter.fastaFile);
    reader.readAll();

    map<string, string> contigs = reader.contigs();
    map<string, string>::iterator iter;
    for(iter = contigs.begin(); iter != contigs.end(); iter++) {
        if(iter->second.length()>=6) {
            adapter.seqsInFasta.push_back(iter->second);
        }
        else {
            cerr << "skip too short adapter sequence in " <<  adapter.fastaFile << " (6bp required): " << iter->second << endl;
        }
    }

    if(adapter.seqsInFasta.size() > 0) {
        adapter.hasFasta = true;
    } else {
        adapter.hasFasta = false;
    }
}

bool Options::validate() {
    if(in.empty()) {
        error_exit("read input should be specified by --in, or enable --stdin if you want to read STDIN");
    } else {
        check_file_valid(in);
    }

    if(outputToSTDOUT) {
        if(!out.empty()) {
            cerr << "In STDOUT mode, ignore the output filename " << out << endl;
            out = "";
        }
    }

    // if output to STDOUT, then...
    if(outputToSTDOUT) {
        if(split.enabled) {
            error_exit("splitting mode cannot work with stdout mode");
        }
    }


    if(!out.empty()) {
        //check_file_writable(out);
        if(dontOverwrite && file_exists(out)) {
            error_exit(out + " already exists and you have set to not rewrite output files by --dont_overwrite");
        }
    }

    if(!failedOut.empty()) {
        if(dontOverwrite && file_exists(failedOut)) {
            error_exit(failedOut + " already exists and you have set to not rewrite output files by --dont_overwrite");
        }
        if(failedOut == out)
            error_exit("--failed_out and --out shouldn't have same file name");
    }

    if(dontOverwrite) {
        if(file_exists(jsonFile)) {
            error_exit(jsonFile + " already exists and you have set to not rewrite output files by --dont_overwrite");
        }
        if(file_exists(htmlFile)) {
            error_exit(htmlFile + " already exists and you have set to not rewrite output files by --dont_overwrite");
        }
    }

    if(compression < 1 || compression > 9)
        error_exit("compression level (--compression) should be between 1 ~ 9, 1 for fastest, 9 for smallest");

    if(readsToProcess < 0)
        error_exit("the number of reads to process (--reads_to_process) cannot be negative");

    if(thread < 1) {
        thread = 1;
    } else if(thread > 16) {
        cerr << "WARNING: fastp uses up to 16 threads although you specified " << thread << endl;
        thread = 16;
    }

    if(trim.front < 0)
        error_exit("trim_front1 (--trim_front1) should be >0, suggest 0 ~ 100");

    if(trim.tail < 0 )
        error_exit("trim_tail1 (--trim_tail1) should be >0, suggest 0 ~ 100");

    if(qualfilter.qualifiedQual - 33 < 0 || qualfilter.qualifiedQual - 33 > 93)
        error_exit("qualitified phred (--qualified_quality_phred) should be 0 ~ 93, suggest 3 ~ 20");

    if(qualfilter.avgQualReq < 0 || qualfilter.avgQualReq  > 93)
        error_exit("average quality score requirement (--mean_qual) should be 0 ~ 93, suggest 5 ~ 30");

    if(qualfilter.unqualifiedPercentLimit < 0 || qualfilter.unqualifiedPercentLimit > 100)
        error_exit("unqualified percent limit (--unqualified_percent_limit) should be 0 ~ 100, suggest 20 ~ 60");

    if(qualfilter.nBasePercentLimit < 0 || qualfilter.nBasePercentLimit > 100)
        error_exit("N base percent limit (--n_percent_limit) should be 0 ~ 100, suggest 5 ~ 20");

    if(qualfilter.nBaseLimit < 0 or qualfilter.nBaseLimit > 1000000)
        error_exit("N base number limit (--n_base_limit) should be 0 ~ 1000000");

    if(lengthFilter.requiredLength < 0 )
        error_exit("length requirement (--length_required) should be >0, suggest >50");

    if(split.enabled ) {
        if(split.digits < 0 || split.digits > 10)
            error_exit("you have enabled splitting output to multiple files, the digits number of file name prefix (--split_prefix_digits) should be 0 ~ 10.");
        
        if(split.byFileNumber) {
            if(split.number < 2 || split.number >= 1000)
                error_exit("you have enabled splitting output by file number, the number of files (--split) should be 2 ~ 999.");
            // thread number cannot be more than the number of file to split
            if(thread > split.number)
                thread = split.number;
        }

        if(split.byFileLines) {
            if(split.size < 1000/4)
                error_exit("you have enabled splitting output by file lines, the file lines (--split_by_lines) should be >= 1000.");
        }
    }

    if(qualityCut.enabledFront || qualityCut.enabledTail) {
        if(qualityCut.windowSizeShared < 1 || qualityCut.windowSizeShared > 1000)
            error_exit("the sliding window size for cutting by quality (--cut_window_size) should be between 1~1000.");
        if(qualityCut.qualityShared < 1 || qualityCut.qualityShared > 30)
            error_exit("the mean quality requirement for cutting by quality (--cut_mean_quality) should be 1 ~ 30, suggest 15 ~ 20.");
        if(qualityCut.windowSizeFront < 1 || qualityCut.windowSizeFront > 1000)
            error_exit("the sliding window size for cutting by quality (--cut_front_window_size) should be between 1~1000.");
        if(qualityCut.qualityFront < 1 || qualityCut.qualityFront > 30)
            error_exit("the mean quality requirement for cutting by quality (--cut_front_mean_quality) should be 1 ~ 30, suggest 15 ~ 20.");
        if(qualityCut.windowSizeTail < 1 || qualityCut.windowSizeTail > 1000)
            error_exit("the sliding window size for cutting by quality (--cut_tail_window_size) should be between 1~1000.");
        if(qualityCut.qualityTail < 1 || qualityCut.qualityTail > 30)
            error_exit("the mean quality requirement for cutting by quality (--cut_tail_mean_quality) should be 1 ~ 30, suggest 13 ~ 20.");
    }

    if(adapter.sequenceStart!="auto" && !adapter.sequenceStart.empty()) {
        // validate adapter sequence for single end adapter trimming
        if(adapter.sequenceStart.length() <= 3)
            error_exit("the sequence of <adapter_sequence> should be longer than 3");

        // validate bases
        for(int i=0; i<adapter.sequenceStart.length(); i++) {
            char c = adapter.sequenceStart[i];
            if(c!='A' && c!='T' && c!='C' && c!='G') {
                error_exit("the adapter <adapter_sequence> can only have bases in {A, T, C, G}, but the given sequence is: " + adapter.sequenceStart);
            }
        }
    }

    if(adapter.edMax <0 || adapter.edMax > 1.0) {
        error_exit("the adapter <distance_threshold> should be 0.0 ~ 1.0, suggest 0.1 ~ 0.3");
    }
    if(adapter.trimmingExtension <0 || adapter.trimmingExtension > 100) {
        error_exit("the adapter <trimming_extension> should be 0 ~ 100, suggest 5 ~ 30");
    }


    return true;
}

bool Options::shallDetectAdapter() {
    if(!adapter.enabled)
        return false;

    return adapter.sequenceStart == "auto" || adapter.sequenceEnd == "auto" ;
}


vector<string> Options::makeListFromFileByLine(string filename) {
    vector<string> ret;
    ifstream file;
    file.open(filename.c_str(), ifstream::in);
    const int maxLine = 1000;
    char line[maxLine];
    cerr << "filter by index, loading " << filename << endl;
    while(file.getline(line, maxLine)){
        // trim \n, \r or \r\n in the tail
        int readed = strlen(line);
        if(readed >=2 ){
            if(line[readed-1] == '\n' || line[readed-1] == '\r'){
                line[readed-1] = '\0';
                if(line[readed-2] == '\r')
                    line[readed-2] = '\0';
            }
        }
        string linestr(line);
        for(int i=0; i<linestr.length(); i++) {
            if(linestr[i] != 'A' && linestr[i] != 'T' && linestr[i] != 'C' && linestr[i] != 'G') {
                error_exit("processing " + filename + ", each line should be one barcode, which can only contain A/T/C/G");
            }
        }
        cerr << linestr << endl;
        ret.push_back(linestr);
    }
    cerr << endl;
    return ret;
}

string Options::getReadStartAdapter(){
    if(adapter.sequenceStart == "" || adapter.sequenceStart == "auto")
        return "unspecified";
    else
        return adapter.sequenceStart;
}

string Options::getReadEndAdapter(){
    if(adapter.sequenceEnd == "" || adapter.sequenceEnd == "auto")
        return "unspecified";
    else
        return adapter.sequenceEnd;
}