#include "jsonreporter.h"

JsonReporter::JsonReporter(Options* opt){
    mOptions = opt;
}

JsonReporter::~JsonReporter(){
}

extern string command;
void JsonReporter::report(FilterResult* result, Stats* preStats1, Stats* postStats1) {
    ofstream ofs;
    ofs.open(mOptions->jsonFile, ifstream::out);
    ofs << "{" << endl;

    long pre_total_reads = preStats1->getReads();

    long pre_total_bases = preStats1->getBases();

    long pre_q20_bases = preStats1->getQ20();

    long pre_q30_bases = preStats1->getQ30();

    long pre_total_gc = preStats1->getGCNumber();

    long post_total_reads = postStats1->getReads();

    long post_total_bases = postStats1->getBases();

    long post_q20_bases = postStats1->getQ20();

    long post_q30_bases = postStats1->getQ30();

    long post_total_gc = postStats1->getGCNumber();

    // summary
    ofs << "\t" << "\"summary\": {" << endl;
    ofs << "\t\t" << "\"fastplong_version\": \""<< FASTPLONG_VER << "\"," << endl;
    ofs << "\t\t" << "\"before_filtering\": {" << endl;
    ofs << "\t\t\t" << "\"total_reads\":" << pre_total_reads << "," << endl; 
    ofs << "\t\t\t" << "\"total_bases\":" << pre_total_bases << "," << endl; 
    ofs << "\t\t\t" << "\"q20_bases\":" << pre_q20_bases << "," << endl; 
    ofs << "\t\t\t" << "\"q30_bases\":" << pre_q30_bases << "," << endl; 
    ofs << "\t\t\t" << "\"q20_rate\":" << (pre_total_bases == 0?0.0:(double)pre_q20_bases / (double)pre_total_bases) << "," << endl; 
    ofs << "\t\t\t" << "\"q30_rate\":" << (pre_total_bases == 0?0.0:(double)pre_q30_bases / (double)pre_total_bases) << "," << endl; 
    ofs << "\t\t\t" << "\"read_mean_length\":" << preStats1->getMeanLength() << "," << endl;
    ofs << "\t\t\t" << "\"gc_content\":" << (pre_total_bases == 0?0.0:(double)pre_total_gc / (double)pre_total_bases)  << endl; 
    ofs << "\t\t" << "}," << endl;

    ofs << "\t\t" << "\"after_filtering\": {" << endl;
    ofs << "\t\t\t" << "\"total_reads\":" << post_total_reads << "," << endl; 
    ofs << "\t\t\t" << "\"total_bases\":" << post_total_bases << "," << endl; 
    ofs << "\t\t\t" << "\"q20_bases\":" << post_q20_bases << "," << endl; 
    ofs << "\t\t\t" << "\"q30_bases\":" << post_q30_bases << "," << endl; 
    ofs << "\t\t\t" << "\"q20_rate\":" << (post_total_bases == 0?0.0:(double)post_q20_bases / (double)post_total_bases) << "," << endl; 
    ofs << "\t\t\t" << "\"q30_rate\":" << (post_total_bases == 0?0.0:(double)post_q30_bases / (double)post_total_bases) << "," << endl; 
    ofs << "\t\t\t" << "\"read_mean_length\":" << postStats1->getMeanLength() << "," << endl;
    ofs << "\t\t\t" << "\"gc_content\":" << (post_total_bases == 0?0.0:(double)post_total_gc / (double)post_total_bases)  << endl; 
    ofs << "\t\t" << "}";

    ofs << endl;

    ofs << "\t" << "}," << endl;

    if(result) {
        ofs << "\t" << "\"filtering_result\": " ;
        result -> reportJson(ofs, "\t");
    }

    if(result && mOptions->adapterCuttingEnabled()) {
        ofs << "\t" << "\"adapter_cutting\": " ;
        result -> reportAdapterJson(ofs, "\t");
    }

    if(result && mOptions->polyXTrimmingEnabled()) {
        ofs << "\t" << "\"polyx_trimming\": " ;
        result -> reportPolyXTrimJson(ofs, "\t");
    }

    if(preStats1) {
        ofs << "\t" << "\"read_before_filtering\": " ;
        preStats1 -> reportJson(ofs, "\t");
    }

    if(postStats1) {
        string name = "read_after_filtering";
        ofs << "\t" << "\"" << name << "\": " ;
        postStats1 -> reportJson(ofs, "\t");
    }

    ofs << "\t\"command\": " << "\"" << command << "\"" << endl;

    ofs << "}";
}