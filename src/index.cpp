#include <Rcpp.h>

using namespace Rcpp;


//  check if all four classes exist
// [[Rcpp::export]]
bool isInteraction(RawVector x1, RawVector x2) {
    size_t n = x1.size();
    bool oneone = false;
    bool zerozero = false;
    bool zeroone = false;
    bool onezero = false;
    bool all = false;
    // setup variables and find the mean
    for (size_t i = 0; i < n; i++) {
        if (x1[i] == 0x01 && x2[i] == 0x01)  {
            zerozero = true;
        } else if (x1[i] == 0x03 && x2[i] == 0x03)  {
            oneone = true;
        } else if (x1[i] == 0x01 && x2[i] == 0x03)  {
             zeroone = true;
        } else if (x1[i] == 0x03 && x2[i] == 0x01)  {
            onezero = true;
        }
        if (zerozero && oneone && zeroone && onezero) {
            all = true;
            break;
        }
    } 
    return all;
}


// make index
// [[Rcpp::export]]
DataFrame interactionIndices(RawMatrix x, CharacterVector chr, NumericVector pos, 
                             double minDist) {
    size_t n = x.ncol();
    bool interaction;
    std::vector<int> inx1;
    std::vector<int> inx2;
    for (size_t i = 0; i < n-1; i ++) {
        for (size_t j = i+1; j < n; j ++) {
            if ((chr[i] != chr[j]) || 
                ((pos[j] - pos[i]) >= minDist)) {
                interaction = isInteraction(x(_, i), x(_, j));
                if (interaction) {
                    inx1.push_back(i + 1);
                    inx2.push_back(j + 1);
                }
            }
        }
    }
    return DataFrame::create(Named("inx1") = inx1,
                             Named("inx2") = inx2);
}


// chack index
// [[Rcpp::export]]
DataFrame testIndices(RawMatrix x, IntegerVector inx1, IntegerVector inx2, 
                      CharacterVector chr, NumericVector pos, double minDist) {
    size_t n = inx1.size();
    bool interaction;
    std::vector<int> inx1new;
    std::vector<int> inx2new;
    size_t i, j;
    for (size_t k = 0; k < n; k ++) {
        i = inx1[k] - 1;
        j = inx2[k] - 1;
        if ((chr[i] != chr[j]) || 
            ((pos[j] - pos[i]) >= minDist)) {
            interaction = isInteraction(x(_, i), x(_, j));
            if (interaction) {
                inx1new.push_back(i + 1);
                inx2new.push_back(j + 1);
            }
        }
    }
    return DataFrame::create(Named("inx1") = inx1new,
                             Named("inx2") = inx2new);
} 
