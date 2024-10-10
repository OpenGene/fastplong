#include "htmlreporter.h"
#include <chrono>
#include <memory.h>

extern string command;

HtmlReporter::HtmlReporter(Options* opt){
    mOptions = opt;
}

HtmlReporter::~HtmlReporter(){
}

void HtmlReporter::outputRow(ofstream& ofs, string key, long v) {
    ofs << "<tr><td class='col1'>" + key + "</td><td class='col2'>" + to_string(v) + "</td></tr>\n";
}

void HtmlReporter::outputRow(ofstream& ofs, string key, string v) {
    ofs << "<tr><td class='col1'>" + key + "</td><td class='col2'>" + v + "</td></tr>\n";
}

string HtmlReporter::formatNumber(long number) {
    double num = (double)number;
    string unit[6] = {"", "K", "M", "G", "T", "P"};
    int order = 0;
    while (num > 1000.0) {
        order += 1;
        num /= 1000.0;
    }

    if (order == 0)
        return to_string(number);
    else
        return to_string(num) + " " + unit[order];
}

string HtmlReporter::getPercents(long numerator, long denominator) {
    if(denominator == 0)
        return "0.0";
    else
        return to_string((double)numerator * 100.0 / (double)denominator);
}

void HtmlReporter::printSummary(ofstream& ofs, FilterResult* result, Stats* preStats1, Stats* postStats1) {

    ofs << endl;
    ofs << "<h3 style='text-align:left;'><a href='https://github.com/OpenGene/fastplong' target='_blank' style='color:#663355;text-decoration:none;'>" + mOptions->reportTitle + "</a><a href='https://github.com/OpenGene/fastplong' target='_blank' style='font-size:-2;text-decoration:none;'>(fastplong version v" + string(FASTPLONG_VER)+ ")</a></h3>"<<endl;
    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('summary')><a name='summary'>Summary</a> </div>\n";
    ofs << "<div id='summary'>\n";

    if(result) {
        ofs << "<div class='subsection_title'>Filtering result</div>\n";
        ofs << "<div id='filtering_result'>\n";
        result -> reportHtml(ofs, preStats1->getReads(), preStats1->getBases());
        ofs << "</div>\n";
    }

    ofs << "</div>\n";
    ofs << "</div>\n";

    /*if(result && mOptions->adapterCuttingEnabled()) {
        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('adapters')><a name='summary'>Adapters</a></div>\n";
        ofs << "<div id='adapters'>\n";

        result->reportAdapterHtml(ofs, pre_total_bases);

        ofs << "</div>\n";
        ofs << "</div>\n";
    }*/

}

void HtmlReporter::report(FilterResult* result, Stats* preStats1, Stats* postStats1) {
    ofstream ofs;
    ofs.open(mOptions->htmlFile, ifstream::out);

    printHeader(ofs);

    printSummary(ofs, result, preStats1, postStats1);

    // basic statistics
    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('basic_stat')><a name='summary'>Basic statistics</a></div>\n";
    ofs << "<table id='basic_stat' class='section_table'>\n";
    ofs << "<tr><td>\n";
    if(preStats1) {
        preStats1 -> reportHtmlBasicInfo(ofs, "Before filtering");
    }
    ofs << "</td><td>\n";
    if(postStats1) {
        postStats1 -> reportHtmlBasicInfo(ofs, "After filtering");
    }
    ofs << "</td></tr>\n";
    ofs << "</table>\n";
    ofs << "</div>\n";

    // median quality histogram
    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('median_qual_stat')><a name='summary'>Median qual histogram</a></div>\n";
    ofs << "<table id='median_qual_stat' class='section_table'>\n";
    ofs << "<tr><td>\n";
    if(preStats1) {
        preStats1 -> reporHtmlMedianQualHist(ofs, "Before filtering");
    }
    ofs << "</td><td>\n";
    if(postStats1) {
        postStats1 -> reporHtmlMedianQualHist(ofs, "After filtering");
    }
    ofs << "</td></tr>\n";
    ofs << "</table>\n";
    ofs << "</div>\n";

    // median quality length density
    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('median_qual_length_density')><a name='summary'>Median qual length density</a></div>\n";
    ofs << "<table id='median_qual_length_density' class='section_table'>\n";
    ofs << "<tr><td>\n";
    if(preStats1) {
        preStats1 -> reporHtmlMedianQualLengthDensity(ofs, "Before filtering");
    }
    ofs << "</td><td>\n";
    if(postStats1) {
        postStats1 -> reporHtmlMedianQualLengthDensity(ofs, "After filtering");
    }
    ofs << "</td></tr>\n";
    ofs << "</table>\n";
    ofs << "</div>\n";

    // quality statistics
    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('quality_stat')><a name='summary'>Quality statistics</a></div>\n";
    ofs << "<table id='quality_stat' class='section_table'>\n";
    ofs << "<tr><td>\n";
    if(preStats1) {
        preStats1 -> reportHtmlQuality(ofs, "Before filtering");
    }
    ofs << "</td><td>\n";
    if(postStats1) {
        postStats1 -> reportHtmlQuality(ofs, "After filtering");
    }
    ofs << "</td></tr>\n";
    ofs << "</table>\n";
    ofs << "</div>\n";

    // content statistics
    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('contents_stat')><a name='summary'>Base contents statistics</a></div>\n";
    ofs << "<table id='contents_stat' class='section_table'>\n";
    ofs << "<tr><td>\n";
    if(preStats1) {
        preStats1 -> reportHtmlContents(ofs, "Before filtering");
    }
    ofs << "</td><td>\n";
    if(postStats1) {
        postStats1 -> reportHtmlContents(ofs, "After filtering");
    }
    ofs << "</td></tr>\n";
    ofs << "</table>\n";
    ofs << "</div>\n";

    // kmer statistics
    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('kmer_stat')><a name='summary'>k-mer statistics</a></div>\n";
    ofs << "<table id='kmer_stat' class='section_table'>\n";
    ofs << "<tr><td>\n";
    if(preStats1) {
        preStats1 -> reportHtmlKMER(ofs, "Before filtering");
    }
    ofs << "</td><td>\n";
    if(postStats1) {
        postStats1 -> reportHtmlKMER(ofs, "After filtering");
    }
    ofs << "</td></tr>\n";
    ofs << "</table>\n";
    ofs << "</div>\n";

    printFooter(ofs);

}

void HtmlReporter::printHeader(ofstream& ofs){
    ofs << "<html><head><meta http-equiv=\"content-type\" content=\"text/html;charset=utf-8\" />";
    ofs << "<title>fastplong report at " + getCurrentSystemTime() + " </title>";
    printJS(ofs);
    printCSS(ofs);
    ofs << "</head>";
    ofs << "<body><div id='container'>";
}

void HtmlReporter::printCSS(ofstream& ofs){
    ofs << "<style type=\"text/css\">" << endl;
    ofs << "td {border:1px solid #dddddd;padding:5px;font-size:12px;}" << endl;
    ofs << "table {border:1px solid #999999;padding:2x;border-collapse:collapse;width:100%}" << endl;
    ofs << ".col1 {width:240px; font-weight:bold;}" << endl;
    ofs << ".adapter_col {width:500px; font-size:10px;}" << endl;
    ofs << "img {padding:30px;}" << endl;
    ofs << "#menu {font-family:Consolas, 'Liberation Mono', Menlo, Courier, monospace;}" << endl;
    ofs << "#menu a {color:#0366d6; font-size:18px;font-weight:600;line-height:28px;text-decoration:none;font-family:-apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol'}" << endl;
    ofs << "a:visited {color: #999999}" << endl;
    ofs << ".alignleft {text-align:left;}" << endl;
    ofs << ".alignright {text-align:right;}" << endl;
    ofs << ".figure {width:680px;height:600px;}" << endl;
    ofs << ".header {color:#ffffff;padding:1px;height:20px;background:#000000;}" << endl;
    ofs << ".section_title {color:#ffffff;font-size:20px;padding:5px;text-align:left;background:#663355; margin-top:10px;}" << endl;
    ofs << ".section_table {width:100%;}" << endl;
    ofs << ".subsection_title {font-size:16px;padding:5px;margin-top:10px;text-align:left;color:#663355}" << endl;
    ofs << "#container {text-align:center;padding:3px 3px 3px 10px;font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}" << endl;
    ofs << ".menu_item {text-align:left;padding-top:5px;font-size:18px;}" << endl;
    ofs << ".highlight {text-align:left;padding-top:30px;padding-bottom:30px;font-size:20px;line-height:35px;}" << endl;
    ofs << "#helper {text-align:left;border:1px dotted #fafafa;color:#777777;font-size:12px;}" << endl;
    ofs << "#footer {text-align:left;padding:15px;color:#ffffff;font-size:10px;background:#663355;font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}" << endl;
    ofs << ".kmer_table {text-align:center;font-size:8px;padding:2px;}" << endl;
    ofs << ".kmer_table td{text-align:center;font-size:8px;padding:0px;color:#ffffff}" << endl;
    ofs << ".sub_section_tips {color:#999999;font-size:10px;padding-left:5px;padding-bottom:3px;}" << endl;
    ofs << "</style>" << endl;
}

void HtmlReporter::printJS(ofstream& ofs){
    ofs << "<script src='https://opengene.org/plotly-1.2.0.min.js'></script>" << endl;
    ofs << "\n<script type='text/javascript'>" << endl;
    ofs << "    window.Plotly || document.write('<script src=\"https://cdn.plot.ly/plotly-1.2.0.min.js\"><\\/script>')" << endl;
    ofs << "</script>" << endl;
    ofs << "\n<script type=\"text/javascript\">" << endl;
    ofs << "    function showOrHide(divname) {" << endl;
    ofs << "        div = document.getElementById(divname);" << endl;
    ofs << "        if(div.style.display == 'none')" << endl;
    ofs << "            div.style.display = 'block';" << endl;
    ofs << "        else" << endl;
    ofs << "            div.style.display = 'none';" << endl;
    ofs << "    }" << endl;
    ofs << "</script>" << endl;
}

const string HtmlReporter::getCurrentSystemTime()
{
  auto tt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  struct tm* ptm = localtime(&tt);
  char date[60] = {0};
  sprintf(date, "%d-%02d-%02d      %02d:%02d:%02d",
    (int)ptm->tm_year + 1900,(int)ptm->tm_mon + 1,(int)ptm->tm_mday,
    (int)ptm->tm_hour,(int)ptm->tm_min,(int)ptm->tm_sec);
  return std::string(date);
}

void HtmlReporter::printFooter(ofstream& ofs){
    ofs << "\n</div>" << endl;
    ofs << "<div id='footer'> ";
    ofs << "<p>"<<command<<"</p>";
    ofs << "fastplong " << FASTPLONG_VER << ", at " << getCurrentSystemTime() << " </div>";
    ofs << "</body></html>";
}
