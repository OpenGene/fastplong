#include <stdlib.h>
#include "filterresult.h"
#include "stats.h"
#include "htmlreporter.h"
#include <memory.h>

FilterResult::FilterResult(Options* opt, bool paired){
    mOptions = opt;
    mTrimmedAdapterRead = 0;
    mTrimmedAdapterBases = 0;
    for(int i=0; i<FILTER_RESULT_TYPES; i++) {
        mFilterReadStats[i] = 0;
    }
    mCorrectionMatrix = new long[64];
    memset(mCorrectionMatrix, 0, sizeof(long)*64);
}

FilterResult::~FilterResult() {
    delete mCorrectionMatrix;
}

void FilterResult::addFilterResult(int result, int readNum) {
    if(result < PASS_FILTER || result >= FILTER_RESULT_TYPES)
        return ;
    mFilterReadStats[result] += readNum;
}

FilterResult* FilterResult::merge(vector<FilterResult*>& list) {
    if(list.size() == 0)
        return nullptr;
    FilterResult* result = new FilterResult(list[0]->mOptions);

    long* target = result->getFilterReadStats();
    // read stats
    for(int i=0; i<list.size(); i++) {
        long* current = list[i]->getFilterReadStats();
        for(int j=0; j<FILTER_RESULT_TYPES; j++) {
            target[j] += current[j];
        }
        result->mTrimmedAdapterRead += list[i]->mTrimmedAdapterRead;
        result->mTrimmedAdapterBases += list[i]->mTrimmedAdapterBases;

        for(int b=0; b<4; b++) {
          result->mTrimmedPolyXReads[b] += list[i]->mTrimmedPolyXReads[b];
          result->mTrimmedPolyXBases[b] += list[i]->mTrimmedPolyXBases[b];
        }

        // merge adapter stats
        map<string, long>::iterator iter;
        for(iter = list[i]->mAdapter.begin(); iter != list[i]->mAdapter.end(); iter++) {
            if(result->mAdapter.count(iter->first) > 0)
                result->mAdapter[iter->first] += iter->second;
            else
                result->mAdapter[iter->first] = iter->second;
        }
    }

    // sort adapters list by adapter length from short to long

    return result;
}

void FilterResult::addReadTrimmed(int bases) {
    mTrimmedAdapterBases  += bases;
    mTrimmedAdapterRead++;
}

void FilterResult::addAdapterTrimmed(string adapter ) {
    if(adapter.empty())
        return;
    
    if(mAdapter.count(adapter) >0 )
        mAdapter[adapter]++;
    else
        mAdapter[adapter] = 1;
}

void FilterResult::addPolyXTrimmed(int base, int length) {
    mTrimmedPolyXReads[base] += 1;
    mTrimmedPolyXBases[base] += length;
}

long FilterResult::getTotalPolyXTrimmedReads() {
  long sum_reads = 0;
  for(int b = 0; b < 4; b++)
    sum_reads += mTrimmedPolyXReads[b];
  return sum_reads;
}

long FilterResult::getTotalPolyXTrimmedBases() {
  long sum_bases = 0;
  for(int b = 0; b < 4; b++)
    sum_bases += mTrimmedPolyXBases[b];
  return sum_bases;
}


void FilterResult::print() {
    cerr <<  "reads passed filter: " << mFilterReadStats[PASS_FILTER] << endl;
    cerr <<  "reads failed due to low quality: " << mFilterReadStats[FAIL_QUALITY] << endl;
    cerr <<  "reads failed due to too many N: " << mFilterReadStats[FAIL_N_BASE] << endl;
    if(mOptions->lengthFilter.enabled) {
        cerr <<  "reads failed due to too short: " << mFilterReadStats[FAIL_LENGTH] << endl;
        if(mOptions->lengthFilter.maxLength > 0)
            cerr <<  "reads failed due to too long: " << mFilterReadStats[FAIL_TOO_LONG] << endl;
    }
    if(mOptions->complexityFilter.enabled) {
        cerr <<  "reads failed due to low complexity: " << mFilterReadStats[FAIL_COMPLEXITY] << endl;
    }
    if(mOptions->adapter.enabled) {
        cerr <<  "reads with adapter trimmed: " << mTrimmedAdapterRead << endl;
        cerr <<  "bases trimmed due to adapters: " << mTrimmedAdapterBases << endl;
    }
    if(mOptions->polyXTrim.enabled) {
        cerr <<  "reads with polyX in 3' end: " << getTotalPolyXTrimmedReads() << endl;
        cerr <<  "bases trimmed in polyX tail: " << getTotalPolyXTrimmedBases() << endl;
    }
}

void FilterResult::reportJson(ofstream& ofs, string padding) {
    ofs << "{" << endl;

    ofs << padding << "\t" << "\"passed_filter_reads\": " << mFilterReadStats[PASS_FILTER] << "," << endl;
    ofs << padding << "\t" << "\"low_quality_reads\": " << mFilterReadStats[FAIL_QUALITY] << "," << endl;
    ofs << padding << "\t" << "\"too_many_N_reads\": " << mFilterReadStats[FAIL_N_BASE] << "," << endl;
    if(mOptions->complexityFilter.enabled)
        ofs << padding << "\t" << "\"low_complexity_reads\": " << mFilterReadStats[FAIL_COMPLEXITY] << "," << endl;
    ofs << padding << "\t" << "\"too_short_reads\": " << mFilterReadStats[FAIL_LENGTH] << "," << endl;
    ofs << padding << "\t" << "\"too_long_reads\": " << mFilterReadStats[FAIL_TOO_LONG] << endl;

    ofs << padding << "}," << endl;
}

void FilterResult::outputAdaptersJson(ofstream& ofs, map<string, long, classcomp>& adapterCounts) {
    map<string, long>::iterator iter;

    long total = 0;
    for(iter = adapterCounts.begin(); iter!=adapterCounts.end(); iter++) {
        total += iter->second;
    }

    if(total == 0)
        return ;

    const double reportThreshold = 0.01;
    const double dTotal = (double)total;
    bool firstItem = true;
    long reported = 0;
    for(iter = adapterCounts.begin(); iter!=adapterCounts.end(); iter++) {
        if(iter->second /dTotal < reportThreshold )
            continue;

        if(!firstItem)
            ofs << ", ";
        else
            firstItem = false;
        ofs << "\"" << iter->first << "\":" << iter->second;

        reported += iter->second;
    }

    long unreported = total - reported;

    if(unreported > 0) {
        if(!firstItem)
            ofs << ", ";
        ofs << "\"" << "others" << "\":" << unreported;
    }
}

void FilterResult::reportAdapterJson(ofstream& ofs, string padding) {
    ofs << "{" << endl;

    ofs << padding << "\t" << "\"adapter_trimmed_reads\": " << mTrimmedAdapterRead << "," << endl;
    ofs << padding << "\t" << "\"adapter_trimmed_bases\": " << mTrimmedAdapterBases << "," << endl;
    ofs << padding << "\t" << "\"read_start_adapter\": \"" << mOptions->getReadStartAdapter() << "\"," << endl;
    ofs << padding << "\t" << "\"read_end_adapter\": \"" << mOptions->getReadEndAdapter() << "\"," << endl;

    ofs << padding << "\t" << "\"read_adapter_counts\": " << "{";
        outputAdaptersJson(ofs, mAdapter);
    ofs << "}";
    ofs << endl;

    ofs << padding << "}," << endl;
}

void writeBaseCountsJson(ofstream& ofs, string pad, string key, long total, long (&counts)[4]) {
  ofs << pad << "\t\"total_" << key << "\": " << total << "," << endl;
  ofs << pad << "\t\"" << key << "\":{";
  for (int b=0; b<4; b++) {
    if(b > 0)
      ofs << ", ";
    ofs << "\"" << ATCG_BASES[b] << "\": " << counts[b];
  }
  ofs << "}";
}

void FilterResult::reportPolyXTrimJson(ofstream& ofs, string padding) {
    ofs << padding << "{" << endl;
    writeBaseCountsJson(ofs, padding, "polyx_trimmed_reads", getTotalPolyXTrimmedReads(), mTrimmedPolyXReads);
    ofs << "," << endl;
    writeBaseCountsJson(ofs, padding, "polyx_trimmed_bases", getTotalPolyXTrimmedBases(), mTrimmedPolyXBases);
    ofs << endl << padding << "}," << endl;
}

/*void FilterResult::reportHtml(ofstream& ofs, long totalReads) {
    const int types = 4;
    const string divName = "filtering_result";
    string labels[4] = {"good_reads", "low_quality_reads", "too_many_N_reads", "too_short_reads"};
    long counts[4] = {mFilterReadStats[PASS_FILTER], mFilterReadStats[FAIL_QUALITY], mFilterReadStats[FAIL_N_BASE], mFilterReadStats[FAIL_LENGTH]};
    
    string json_str = "var data=[";
    json_str += "{values:[" + Stats::list2string(counts, types) + "],";
    json_str += "labels:['good_reads', 'low_quality_reads', 'too_many_N_reads', 'too_short_reads'],";
    json_str += "textinfo: 'none',";
    json_str += "type:'pie'}];\n";
    string title = "Filtering statistics of sampled " + to_string(totalReads) + " reads";
    json_str += "var layout={title:'" + title + "', width:800, height:400};\n";
    json_str += "Plotly.newPlot('" + divName + "', data, layout);\n";

    ofs << "<div class='figure' id='" + divName + "'></div>\n";
    ofs << "\n<script type=\"text/javascript\">" << endl;
    ofs << json_str;
    ofs << "</script>" << endl;
} */

void FilterResult::reportHtml(ofstream& ofs, long totalReads, long totalBases) {
    double total = (double)totalReads;
    ofs << "<table class='summary_table'>\n";
    HtmlReporter::outputRow(ofs, "reads passed filters:", HtmlReporter::formatNumber(mFilterReadStats[PASS_FILTER]) + " (" + to_string(mFilterReadStats[PASS_FILTER] * 100.0 / total) + "%)");

    HtmlReporter::outputRow(ofs, "reads with low quality:", HtmlReporter::formatNumber(mFilterReadStats[FAIL_QUALITY]) + " (" + to_string(mFilterReadStats[FAIL_QUALITY] * 100.0 / total) + "%)");
    HtmlReporter::outputRow(ofs, "reads with too many N:", HtmlReporter::formatNumber(mFilterReadStats[FAIL_N_BASE]) + " (" + to_string(mFilterReadStats[FAIL_N_BASE] * 100.0 / total) + "%)");
    if(mOptions->lengthFilter.enabled) {
        HtmlReporter::outputRow(ofs, "reads too short:", HtmlReporter::formatNumber(mFilterReadStats[FAIL_LENGTH]) + " (" + to_string(mFilterReadStats[FAIL_LENGTH] * 100.0 / total) + "%)");
        if(mOptions->lengthFilter.maxLength > 0)
            HtmlReporter::outputRow(ofs, "reads too long:", HtmlReporter::formatNumber(mFilterReadStats[FAIL_TOO_LONG]) + " (" + to_string(mFilterReadStats[FAIL_TOO_LONG] * 100.0 / total) + "%)");
    }
    if(mOptions->complexityFilter.enabled)
        HtmlReporter::outputRow(ofs, "reads with low complexity:", HtmlReporter::formatNumber(mFilterReadStats[FAIL_COMPLEXITY]) + " (" + to_string(mFilterReadStats[FAIL_COMPLEXITY] * 100.0 / total) + "%)");
    ofs << "</table>\n";
}

void FilterResult::reportAdapterHtml(ofstream& ofs, long totalBases) {
    ofs << "<div class='subsection_title' onclick=showOrHide('read1_adapters')>Adapter or bad ligation of read1</div>\n";
    ofs << "<div id='read1_adapters'>\n";
    outputAdaptersHtml(ofs, mAdapter, totalBases);
    ofs << "</div>\n";
}

void FilterResult::outputAdaptersHtml(ofstream& ofs, map<string, long, classcomp>& adapterCounts, long totalBases) {

    map<string, long>::iterator iter;

    long total = 0;
    long totalAdapterBases = 0;
    for(iter = adapterCounts.begin(); iter!=adapterCounts.end(); iter++) {
        total += iter->second;
        totalAdapterBases += iter->first.length() * iter->second;
    }

    double frac = (double)totalAdapterBases / (double)totalBases;


    if(frac < 0.01) {
        ofs << "<div class='sub_section_tips'>The input has little adapter percentage (~" << to_string(frac*100.0) << "%), probably it's trimmed before.</div>\n";
    }

    if(total == 0)
        return ;

    ofs << "<table class='summary_table'>\n";
    ofs << "<tr><td class='adapter_col' style='font-size:14px;color:#ffffff;background:#556699'>" << "Sequence" << "</td><td class='col2' style='font-size:14px;color:#ffffff;background:#556699'>" << "Occurrences" << "</td></tr>\n";

    const double reportThreshold = 0.01;
    const double dTotal = (double)total;
    long reported = 0;
    for(iter = adapterCounts.begin(); iter!=adapterCounts.end(); iter++) {
        if(iter->second /dTotal < reportThreshold )
            continue;

        ofs << "<tr><td class='adapter_col'>" << iter->first << "</td><td class='col2'>" << iter->second << "</td></tr>\n";

        reported += iter->second;
    }

    long unreported = total - reported;

    if(unreported > 0) {
        string tag = "other adapter sequences";
        if(reported == 0)
            tag = "all adapter sequences";
        ofs << "<tr><td class='adapter_col'>" << tag << "</td><td class='col2'>" << unreported << "</td></tr>\n";
    }
    ofs << "</table>\n";
}
