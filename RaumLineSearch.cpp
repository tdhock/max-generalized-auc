#include "Rcpp.h"
#include "aumLineSearch.h"
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List aumLineSearch(const DataFrame df, const double initialAum) {
    // extract columnds from dataframe
    NumericVector fpDiff = df["fp.diff"];
    NumericVector fnDiff = df["fn.diff"];
    NumericVector intercept = df["intercept"];
    NumericVector slope = df["slope"];
    int lineCount = df.nrow();
    cout << "Running line intersection with " << lineCount << " lines." << endl;

    // build lines
    vector<Line> lines;
    lines.reserve(lineCount);
    for (int i = 0; i < lineCount; i++) {
        Line line = Line { .intercept = intercept[i], .slope = slope[i] };
        lines.push_back(line);
    }

    double FP[50] = {0};
    double FN[50] = {0};
    double M[50] = {0};
    NumericVector stepSizeVec(50);
    NumericVector aumVec(50);

    lineSearch(
        &lines[0],
        lineCount,
        &fpDiff[0],
        &fnDiff[0],
        initialAum,
        20,
        FP, FN, M,
        &stepSizeVec[0],
        &aumVec[0]
    );
    
    return List::create(Named("aum", aumVec), Named("step.size", stepSizeVec));
}
