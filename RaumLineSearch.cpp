#include "Rcpp.h"
#include "aumLineSearch.h"
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List aumLineSearch(const DataFrame df, const double initialAum, int maxIterations) {
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

    // stack allocated, not returned
    double FP[maxIterations];
    double FN[maxIterations];
    double M[maxIterations];
    // heap allocated so we can return it to R
    NumericVector stepSizeVec(maxIterations);
    NumericVector aumVec(maxIterations);

    lineSearch(
        &lines[0],
        lineCount,
        &fpDiff[0],
        &fnDiff[0],
        initialAum,
        maxIterations,
        FP, FN, M,
        &stepSizeVec[0],
        &aumVec[0]
    );
    
    return List::create(Named("aum", aumVec), Named("step.size", stepSizeVec));
}
