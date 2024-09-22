#include "processor.h"
#include "seprocessor.h"

Processor::Processor(Options* opt){
    mOptions = opt;
}


Processor::~Processor(){
}

bool Processor::process() {
    SingleEndProcessor p(mOptions);
    p.process();

    return true;
}