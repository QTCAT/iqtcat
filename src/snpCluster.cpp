#include <Rcpp.h>

using namespace Rcpp;


//  2 x 2 correlation 
// [[Rcpp::export]]
double cor_x(RawVector x1, RawVector x2, RawVector y1, RawVector y2) {
     size_t n = x1.size();
    IntegerVector X(n), Y(n);
    double ex(0), ey(0), xt(0), yt(0), sxx(0), syy(0), sxy(0);
    // setup variables and find the mean
    for (size_t i = 0; i < n; i ++) {
        if (x2[i] != 0x01)  {
            if (x1[i] == 0x03)  {
                if (x2[i] == 0x03) {
                    X[i] = 4;
                } else if (x2[i] == 0x02) {
                    X[i] = 2;
                }
            } else if (x1[i] == 0x02)  {
                if (x2[i] == 0x03) {
                    X[i] = 2;
                } else if (x2[i] == 0x02) {
                    X[i] = 1;
                }
            }
        } 
        if (y2[i] != 0x01)  {
            if (y1[i] == 0x03)  {
                if (y2[i] == 0x03) {
                    Y[i] = 4;
                } else if (y2[i] == 0x02) {
                    Y[i] = 2;
                }
            } else if (y1[i] == 0x02)  {
                if (y2[i] == 0x03) {
                    Y[i] = 2;
                } else if (y2[i] == 0x02) {
                    Y[i] = 1;
                }
            }
        }
        ex += X[i];
        ey += Y[i];
    }
    ex /= n;
    ey /= n;
    // correlation coefficent
    for (size_t i = 0; i < n; i ++) {
        xt = X[i] - ex;
        yt = Y[i] - ey;
        sxx += xt * xt;
        syy += yt * yt;
        sxy += xt * yt;
    }
    double cor = sxy / (sqrt(sxx * syy) + 1e-16);
    return cor;
}


//  2 x 2 correlation distance
// [[Rcpp::export]]
double corDist_x(RawVector x1, RawVector x2, RawVector y1, RawVector y2) {
    double dcor = 1 - std::fabs(cor_x(x1, x2, y1, y2));
    return dcor;
}


// distance from all marker
// [[Rcpp::export]]
NumericVector corDists_x(RawMatrix x, IntegerVector inx1, IntegerVector inx2) {
    size_t nXact = inx1.size();
    NumericVector dist(nXact * (nXact - 1) / 2);
    size_t i1, i2, j1, j2;
    size_t count = 0;
    for (size_t i = 0; i < nXact - 1; i ++) {
         i1 = inx1[i] - 1;
         i2 = inx2[i] - 1;
        for (size_t j = i + 1; j < nXact; j ++) {
            j1 = inx1[j] - 1;
            j2 = inx2[j] - 1;
            dist[count] = corDist_x(x(_, i1), x(_, i2), 
                                    x(_, j1), x(_, j2));
            count ++;
        }
    }
   return dist;
}


// cluster with zero distance
// [[Rcpp::export]]
List corPreIdenticals_x(RawMatrix x, IntegerVector inx1, IntegerVector inx2, 
                        const int step) {
    size_t nXact = inx1.size();
    std::map< int, std::vector<int> > preclust;
    size_t k = 1;
    std::vector<int> premedoInx;
    size_t i1, i2, j0, j1, j2;
    // medoids equally distributed over the search space
    if (nXact >= (step * 2)) {
        bool smallDist = false;
        size_t startStep = (int)step * .5;
        premedoInx.push_back(startStep);
        double medodist = 0;
        for (size_t i = startStep + step; i < nXact; i += step) {
            i1 = inx1[i] - 1;
            i2 = inx2[i] - 1;
            for(size_t j = startStep; j < i; j += step) {
                j1 = inx1[j] - 1;
                j2 = inx2[j] - 1;
                medodist = corDist_x(x(_, i1), x(_, i2), 
                                     x(_, j1), x(_, j2));
                if (medodist < .05) {
                    smallDist = true;
                    break;
                }
            }
            if (smallDist) {
                smallDist = false;
                continue;
            }
            premedoInx.push_back(i);
        }
        k = premedoInx.size();
    }
    // cluster all variables to the nearest medoid
    if (k > 1) {
        NumericVector predist(k);
        for (size_t i = 0; i < nXact; i ++) {
            i1 = inx1[i] - 1;
            i2 = inx2[i] - 1;
            for (size_t j = 0; j < k; j ++) {
                j0 = premedoInx[j];
                j1 = inx1[j0] - 1;
                j2 = inx2[j0] - 1;
                predist[j] = corDist_x(x(_, i1), x(_, i2), 
                                       x(_, j1), x(_, j2));
            }
            preclust[which_min(predist)].push_back(i + 1);
        }
    } else {
        for (size_t i = 0; i < nXact; i ++) {
            preclust[0].push_back(i + 1);
        }
    }
    return wrap(preclust);
}


// find all zero clusters in between each of the pre-clusters
// [[Rcpp::export]]
List corIdenticals_x(RawMatrix x, IntegerVector inx1, IntegerVector inx2, 
                   IntegerVector clustIdx) {
    size_t nCluster = clustIdx.size();
    double dist = 0;
    size_t clustCount = 1;
    std::vector<int> grp1, grp2, medoInx1, medoInx2;
    size_t grpMedo;
    IntegerVector clusters(nCluster);
    size_t i0, i1, i2, j0, j1, j2;
    for (size_t i = 0; i < nCluster; i ++) {
        i0 = clustIdx[i] - 1;
        if (clusters[i] == 0) {
            clusters[i] = clustCount;
            i1 = inx1[i0] - 1;
            i2 = inx2[i0] - 1;
            grp1.push_back(i1);
            grp2.push_back(i2);
            for (size_t j = i + 1; j < nCluster; j ++) {
                j0 = clustIdx[j] - 1;
                if (clusters[j] == 0) {
                    j1 = inx1[j0] - 1;
                    j2 = inx2[j0] - 1;
                    dist = corDist_x(x(_, i1), x(_, i2), 
                                     x(_, j1), x(_, j2));
                    if (dist <= 1e-7) {
                        clusters[j] = clustCount;
                        grp1.push_back(j1);
                        grp2.push_back(j2);
                    }
                }
            }
            grpMedo = grp1.size() / 2;
            medoInx1.push_back(grp1[grpMedo] + 1);
            medoInx2.push_back(grp2[grpMedo] + 1);
            clustCount ++;
            grp1.clear();
            grp2.clear();
        }
    }
    DataFrame medoids = DataFrame::create(Named("inx1") = medoInx1,
                                          Named("inx2") = medoInx2);
    return List::create(clusters, medoids);
}


// join data to one object
// [[Rcpp::export]]
List joinCorIdenticals_x(int n, List preclust, List ClustMedo) {
     IntegerVector clusters(n);
     std::vector<int> medoInx1, medoInx2;
     List SubClustMedo, medo;
     IntegerVector clustIdx, clust, medo1, medo2;
     size_t count = 0, j0 = 0, j1 = 0, j2 = 0;
     for (size_t i = 0; i < preclust.size(); i ++) {
         clustIdx = preclust[i];
         SubClustMedo = ClustMedo[i];
         clust = SubClustMedo[0];
         medo = SubClustMedo[1];
         medo1 = medo[0];
         medo2 = medo[1];
         for (size_t j = 0; j < clustIdx.size(); j ++) {
             j0 = clustIdx[j] - 1;
             clusters[j0] = clust[j] + count;
         }
         for (size_t j = 0; j < medo1.size(); j ++) {
             medoInx1.push_back(medo1[j]);
             medoInx2.push_back(medo2[j]);
         }
         count = medoInx1.size();
     }
     DataFrame medoids = DataFrame::create(Named("inx1") = medoInx1,
                                           Named("inx2") = medoInx2);
     return List::create(clusters, medoids);
}


// clustering of IDB by correlation dictamce with CLARANS
// [[Rcpp::export]]
List corClarans_x(RawMatrix x, IntegerVector inx1, IntegerVector inx2, 
                 const int k, const int maxNeigbours) {
    RNGScope scope;
    Environment base("package:base");
    Function sample_int = base["sample.int"];
    size_t nXact = inx1.size();
    NumericMatrix medoDist(nXact, k);
    IntegerVector medoInx1(k);
    IntegerVector medoInx2(k);
    double dist = 0;
    IntegerVector clusters(nXact);
    NumericVector pdist(nXact, 1.0);
    size_t i0, i1, i2, j1, j2, l0;
    // starting medoid and clusters
    for (size_t i = 0; i < k; i ++) {
        i0 = as<size_t>(sample_int(nXact, 1, false, pdist)) - 1;
        i1 = inx1[i0] - 1;
        i2 = inx2[i0] - 1;
        medoInx1[i] = i1;
        medoInx2[i] = i2;
        for (size_t j = 0; j < nXact; j ++) {
            j1 = inx1[j] - 1;
            j2 = inx2[j] - 1;
            dist = corDist_x(x(_, i1), x(_, i2), x(_, j1), x(_, j2));
            medoDist(j, i) = dist;
        }
    }
    size_t clust = 0;
    double costs = 0;
    double objective = 0;
    for (size_t i= 0; i < nXact; i ++) {
        // cluster membership of objects
        clust = which_min(medoDist(i, _));
        clusters[i] = clust;
        costs += medoDist(i, clust);
    }
    objective = costs / nXact;
    IntegerVector ranVar(maxNeigbours);
    double minDist = 0;
    NumericVector objectDist(nXact);
    NumericVector objectMinDist(nXact);
    double object_cost = 0;
    // iteration of inner clarans loop
    size_t i = 0;
    while (i < maxNeigbours) {
        if ( i == 0) {
            ranVar = floor(runif(maxNeigbours, 0, nXact - 1e-15));
        }
        // replacing medoids with objects if costs are smaller
        i0 = ranVar[i];
        i1 = inx1[i0] - 1;
        i2 = inx2[i0] - 1;
        l0 = clusters[i0];
        for (size_t j = 0; j < nXact; j ++) {
            j1 = inx1[j] - 1;
            j2 = inx2[j] - 1;
            dist = corDist_x(x(_, i1), x(_, i2), x(_, j1), x(_, j2));
            objectDist[j] = dist;
            minDist = dist;
            for (size_t l = 0; l < k; l ++) {
                if (l != l0) {
                    minDist = std::min(minDist, medoDist(j, l));
                }
            }
            objectMinDist[j] = minDist;
        }
        object_cost = sum(objectMinDist);
        if (object_cost < costs) {
            // new medoids
            medoDist(_, l0) = objectDist;
            medoInx1[l0] = i1 + 1;
            medoInx2[l0] = i2 + 1;
            //  new costs
            costs = object_cost;
            objective = costs / nXact;
            // new membership to cluster
            for (size_t j = 0; j < nXact; j ++) {
                clusters[j] = which_min(medoDist(j, _));
            }
            // restarting loop 
            i = 0;
        } else {
            i ++;
        }
    }
    // Result as list
    DataFrame medoids = DataFrame::create(Named("inx1") = medoInx1,
                                          Named("inx2") = medoInx2);
    return List::create(clusters + 1, medoids, objective);
} // corClarans_x
