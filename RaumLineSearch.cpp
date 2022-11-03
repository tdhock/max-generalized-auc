#include "Rcpp.h"
#include "aumLineSearch.h"
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
int aumLineSearch(const DataFrame df) {
    // extract columnds from dataframe
    NumericVector fpDiff = df["fp.diff"];
    NumericVector fnDiff = df["fn.diff"];
    NumericVector intercept = df["intercept"];
    NumericVector slope = df["slope"];
    int lineCount = df.size();
    
    // build lines
    vector<Line> lines;
    lines.reserve(lineCount);
    for (int i = 0; i < lineCount; i++) {
        Line line = Line { .intercept = intercept[i], .slope = slope[i] };
        lines.push_back(line);
    }
    
    double FP[10] = {0};
    double FN[10] = {0};
    double M[10] = {0};
    NumericVector stepSizeVec;
    NumericVector aumVec;
    
    lineSearch(
        &lines[0],
        lineCount,
        &fpDiff[0],
        &fnDiff[0],
        FP, FN, M,
        20,
        &stepSizeVec[0],
        &aumVec[0]
    );
    
    return 0;
}
