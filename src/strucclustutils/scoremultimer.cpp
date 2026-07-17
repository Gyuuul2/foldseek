#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "LocalParameters.h"
#include "StructureUtil.h"
#include "Coordinate16.h"
#include "MultimerUtil.h"
#include "set"
#include "unordered_set"
#include "tmalign/basic_fun.h"
#include "FileUtil.h"
#include "LDDT.h"
#include <map>
#include <algorithm>

#ifdef OPENMP
#include <omp.h>
#endif
#define INTERFACE_THRESHOLD 8

// carrying chainToChainAlignments from the same query and target complex
struct SearchResult {
    SearchResult() {}
    // alnVec is moved in: for a very large complex it holds up to qChains*dbChains coordinate-heavy
    // entries, so a deep copy here would double the read-path peak memory. The caller (currAlns) is
    // cleared before its next use, so the move is safe and the contents are byte-identical.
    SearchResult(std::vector<unsigned int> &qChainKeys, unsigned int qResidueLen, std::vector<unsigned int> &dbChainKeys, unsigned int dbResidueLen, std::vector<ChainToChainAln> &alnVec) : qChainKeys(qChainKeys), dbChainKeys(dbChainKeys), qResidueLen(qResidueLen), dbResidueLen(dbResidueLen), alnVec(std::move(alnVec)) {}

    std::vector<unsigned int> qChainKeys;
    std::vector<unsigned int> dbChainKeys;
    unsigned int qResidueLen;
    unsigned int dbResidueLen;
    std::vector<ChainToChainAln> alnVec;

    void resetDbComplex(std::vector<unsigned int> &chainKeys, unsigned int residueLen) {
        dbChainKeys = chainKeys;
        dbResidueLen = residueLen;
    }

    void standardize(int MonomerIncludeMode) {
        if (dbResidueLen == 0)
            alnVec.clear();

        if (MonomerIncludeMode == SKIP_MONOMERS && dbChainKeys.size() < MULTIPLE_CHAINED_COMPLEX)
            alnVec.clear();

        if (alnVec.empty())
            return;

        double length = alnVec.size();
        for (size_t i = 0; i < SIZE_OF_SUPERPOSITION_VECTOR; i++) {
            double mean = 0.0;
            double var = 0.0;
            for (auto &aln : alnVec) {
                mean += aln.superposition[i] / length;
            }
            for (auto &aln : alnVec) {
                var += std::pow(aln.superposition[i] - mean, 2) / length;
            }
            double sd = std::sqrt(var);
            double cv = (abs(mean) > TOO_SMALL_MEAN) ? sd / std::abs(mean) : sd;
            for (auto &aln : alnVec) {
                aln.superposition[i] = cv < TOO_SMALL_CV ? FILTERED_OUT : (aln.superposition[i] - mean) / sd;
            }
        }
    }

    void clear() {
        alnVec.clear();
        qChainKeys.clear();
        dbChainKeys.clear();
    }
};

// compute complex tm score
// carrying final output lines
struct Assignment {
    Assignment() {}
    Assignment(unsigned int qLength, unsigned int dbLength): qResidueLength(qLength), dbResidueLength(dbLength), matches(0) {}
    unsigned int qResidueLength;
    unsigned int dbResidueLength;
    unsigned int matches;
    std::vector<float> qCaXVec;
    std::vector<float> qCaYVec;
    std::vector<float> qCaZVec;
    std::vector<float> dbCaXVec;
    std::vector<float> dbCaYVec;
    std::vector<float> dbCaZVec;
    std::vector<unsigned int> matchLenVec; 
    double qTmScore;
    double dbTmScore;
    std::string tString;
    std::string uString;
    std::string backtrace;
    resultToWrite_t assignmentResult;
    TMaligner::TMscoreResult tmResult;
    std::vector<resultToWriteWithKey_t> chainToChainResults;
    resultToWrite_t filterResult;
    unsigned int assignmentId;

    // Matched CA coordinates are appended by ComplexScorer::appendClusterAln (it re-decodes them from
    // the CA DB), since ChainToChainAln no longer carries the coordinate vectors.

    void reset() {
        matches = 0;
        qCaXVec.clear();
        qCaYVec.clear();
        qCaZVec.clear();
        dbCaXVec.clear();
        dbCaYVec.clear();
        dbCaZVec.clear();
        uString.clear();
        tString.clear();
        chainToChainResults.clear();
        assignmentResult.clear();
        filterResult.clear();
    }

    bool getTmScore(TMaligner &tmAligner) {
        // for safety
        if (matches == 0) {
            return false;
        }
        backtrace = std::string(matches, 'M');
        unsigned int normLen = std::min(qResidueLength, dbResidueLength);
        tmAligner.initQuery(&qCaXVec[0], &qCaYVec[0], &qCaZVec[0], NULL, matches);
        tmResult = tmAligner.computeTMscore(&dbCaXVec[0], &dbCaYVec[0], &dbCaZVec[0],matches,0,0, backtrace,normLen);
        qTmScore = tmResult.tmscore * normLen / qResidueLength;
        dbTmScore = tmResult.tmscore * normLen / dbResidueLength;
        return true;
    }

    void updateResultToWriteLines() {
        char sep = ',';
        char tab = '\t';
        tString.append(std::to_string(tmResult.t[0]) + sep + std::to_string(tmResult.t[1]) + sep + std::to_string(tmResult.t[2]));
        uString.append(std::to_string(tmResult.u[0][0]) + sep + std::to_string(tmResult.u[0][1]) + sep + std::to_string(tmResult.u[0][2]) + sep);
        uString.append(std::to_string(tmResult.u[1][0]) + sep + std::to_string(tmResult.u[1][1]) + sep + std::to_string(tmResult.u[1][2]) + sep);
        uString.append(std::to_string(tmResult.u[2][0]) + sep + std::to_string(tmResult.u[2][1]) + sep + std::to_string(tmResult.u[2][2]));
        assignmentResult = tab + std::to_string(qTmScore) + tab + std::to_string(dbTmScore) + tab + uString + tab + tString;
        uString.clear();
        tString.clear();
    }

    void getChainToChainResult(const unsigned int qChainKey, std::string &chainToChainResult) const {
        for (auto & chainResultWithKey: chainToChainResults) {
            if (chainResultWithKey.first != qChainKey) {
                continue;
            }
            chainToChainResult = chainResultWithKey.second;
            return;
        }
    }
};

struct NeighborsWithDist {
    NeighborsWithDist(unsigned int neighbor, float dist) : neighbor(neighbor), dist(dist) {}
    unsigned int neighbor;
    float dist;
};

bool compareChainToChainAlnByDbComplexId(const ChainToChainAln &first, const ChainToChainAln &second) {
    if (first.dbChain.complexId < second.dbChain.complexId) {
        return true;
    }
    if (first.dbChain.complexId > second.dbChain.complexId) {
        return false;
    }

    if (first.qChain.chainKey < second.qChain.chainKey) {
        return true;
    }
    if (first.qChain.chainKey > second.qChain.chainKey) {
        return false;
    }

    if (first.dbChain.chainKey < second.dbChain.chainKey) {
        return true;
    }
    if (first.dbChain.chainKey > second.dbChain.chainKey) {
        return false;
    }
    return false;
}

bool compareChainToChainAln(const ChainToChainAln &first, const ChainToChainAln &second) {
    if (first.qChain.chainKey < second.qChain.chainKey) {
        return true;
    }
    if (first.qChain.chainKey > second.qChain.chainKey) {
        return false;
    }

    if (first.dbChain.chainKey < second.dbChain.chainKey) {
        return true;
    }
    if (first.dbChain.chainKey > second.dbChain.chainKey) {
        return false;
    }
    return false;
}

bool compareAssignment(const Assignment &first, const Assignment &second) {
    if (first.qTmScore > second.qTmScore)
        return true;
    if (first.qTmScore < second.qTmScore)
        return false;
    if (first.dbTmScore > second.dbTmScore)
        return true;
    if (first.dbTmScore < second.dbTmScore)
        return false;
    return false;
}

bool compareNeighborWithDist(const NeighborsWithDist &first, const NeighborsWithDist &second) {
    if (first.dist < second.dist)
        return true;
    if (first.dist > second.dist)
        return false;
    return false;
}

void getResult(std::string &result, std::string &currentResult, const Assignment &assignment) {
    // currentResult.append("\t" + assignment.assignmentResult + "\t" + std::to_string(assignment.assignmentId));
    // if (!assignment.filterResult.empty()) {
    //     currentResult.append("\t" + assignment.filterResult);
    // }
    currentResult.append("\t" + assignment.assignmentResult);
    if (!assignment.filterResult.empty()) {
        currentResult.append("\t" + assignment.filterResult);
    }
    currentResult.append(  + "\t" + std::to_string(assignment.assignmentId));
    currentResult.append("\n");
    result.append(currentResult);
}

// Below this many alignments the brute-force O(n) neighbor scan is cheaper than building a tree.
// (Memory is O(n) on both paths; this only trades tree-build overhead vs scan cost.)
#define KD_TREE_MIN_POINTS 64

// Adaptive neighbor search for DBSCAN. For large alignment counts we build a k-d tree over the
// 12-D superposition vectors so an eps-radius neighbor query costs ~O(log n) instead of O(n),
// with only O(n) memory (no O(n^2) distance matrix). Descent uses a small margin so a true
// neighbor (by the exact ChainToChainAln::getDistance) is never pruned away; the exact distance
// test decides membership, so the returned set is identical to the brute-force scan.
class SuperpositionKdTree {
public:
    void build(const std::vector<ChainToChainAln> &alns) {
        alns_ = &alns;
        size_t n = alns.size();
        point_.assign(n, 0);
        splitAxis_.assign(n, 0);
        left_.assign(n, -1);
        right_.assign(n, -1);
        order_.resize(n);
        for (size_t i = 0; i < n; i++) {
            order_[i] = (unsigned int) i;
        }
        nodeCount_ = 0;
        root_ = buildRec(0, n, 0);
    }

    // Appends every j (j != queryIdx) with getDistance(queryIdx, j) < eps. Order is unspecified;
    // the caller sorts to match the brute-force ordering.
    void rangeQuery(unsigned int queryIdx, float eps, std::vector<unsigned int> &out) const {
        queryRec(root_, queryIdx, eps, out);
    }

private:
    const std::vector<ChainToChainAln> *alns_;
    std::vector<unsigned int> point_;   // node id -> alignment index
    std::vector<int> splitAxis_;
    std::vector<int> left_;
    std::vector<int> right_;
    std::vector<unsigned int> order_;   // working index buffer (partitioned during build)
    int root_;
    int nodeCount_;

    int buildRec(size_t lo, size_t hi, int depth) {
        if (lo >= hi) {
            return -1;
        }
        const std::vector<ChainToChainAln> &alns = *alns_;
        int axis = depth % (int) SIZE_OF_SUPERPOSITION_VECTOR;
        size_t mid = lo + (hi - lo) / 2;
        std::nth_element(order_.begin() + lo, order_.begin() + mid, order_.begin() + hi,
            [&alns, axis](unsigned int a, unsigned int b) {
                return alns[a].superposition[axis] < alns[b].superposition[axis];
            });
        int node = nodeCount_++;
        point_[node] = order_[mid];
        splitAxis_[node] = axis;
        left_[node] = buildRec(lo, mid, depth + 1);
        right_[node] = buildRec(mid + 1, hi, depth + 1);
        return node;
    }

    void queryRec(int node, unsigned int q, float eps, std::vector<unsigned int> &out) const {
        if (node < 0) {
            return;
        }
        const std::vector<ChainToChainAln> &alns = *alns_;
        unsigned int p = point_[node];
        if (p != q && alns[q].getDistance(alns[p]) < eps) {
            out.emplace_back(p);
        }
        int axis = splitAxis_[node];
        double diff = (double) alns[q].superposition[axis] - (double) alns[p].superposition[axis];
        // conservative: also descend the far side unless it is provably beyond eps
        double reach = (double) eps + ((double) eps * 1e-3 + 1e-4);
        if (diff <= 0.0) {
            queryRec(left_[node], q, eps, out);
            if (-diff < reach) {
                queryRec(right_[node], q, eps, out);
            }
        } else {
            queryRec(right_[node], q, eps, out);
            if (diff < reach) {
                queryRec(left_[node], q, eps, out);
            }
        }
    }
};

class DBSCANCluster {
public:
    DBSCANCluster(SearchResult &searchResult, std::set<cluster_t> &finalClusters, float minCov) : searchResult(searchResult), finalClusters(finalClusters) {
        useKdTree = false;
        cLabel = 0;
        minimumClusterSize = std::ceil((float) searchResult.qChainKeys.size() * minCov);
        maximumClusterSize = std::min(searchResult.qChainKeys.size(), searchResult.dbChainKeys.size());
        maximumClusterNum = searchResult.alnVec.size() / maximumClusterSize;
        prevMaxClusterSize = 0;
        maxDist = FLT_MIN;
        minDist = FLT_MAX;
        learningRate = LEARNING_RATE;
    }

    bool getAlnClusters() {
        // if Query or Target is a Monomer Complex.
        if (std::min(searchResult.qChainKeys.size(), searchResult.dbChainKeys.size()) < MULTIPLE_CHAINED_COMPLEX)
            return earlyStopForMonomers();

        // rbh filter
        filterAlnsByRBH();
        initEpsRange();
        // For large alignment counts, a k-d tree makes neighbor queries ~O(log n) with O(n) memory;
        // small counts stay on the brute-force scan. alnVec is fixed from here (only labels change).
        useKdTree = searchResult.alnVec.size() > KD_TREE_MIN_POINTS;
        if (useKdTree) {
            kdTree.build(searchResult.alnVec);
        }
        // To skip DBSCAN clustering when alignments are few enough.
        if (searchResult.alnVec.size() <= maximumClusterSize)
            return checkClusteringNecessity();

        return runDBSCAN();
    }

private:
    SearchResult &searchResult;
    float eps;
    float maxDist;
    float minDist;
    float learningRate;
    unsigned int cLabel;
    unsigned int maximumClusterNum;
    unsigned int prevMaxClusterSize;
    unsigned int currMaxClusterSize;
    unsigned int maximumClusterSize;
    unsigned int minimumClusterSize;
    std::vector<unsigned int> neighbors;
    std::vector<unsigned int> neighborsOfCurrNeighbor;
    // Membership set of alignment indices already pulled into the current cluster, kept as a reused
    // epoch-stamped array instead of a hash set: foundEpoch[i] == curEpoch means "already added". Bumping
    // curEpoch per cluster is an O(1) clear, so there is no per-cluster hashing or reallocation.
    std::vector<size_t> foundEpoch;
    size_t curEpoch;
    std::vector<NeighborsWithDist> neighborsWithDist;
    // "chain already seen" sets, kept as reused epoch-stamped arrays indexed by local chain id instead
    // of hash sets of chain keys: qFoundEpoch[localId] == chainEpoch means seen. Bumping chainEpoch is
    // an O(1) clear. Sized to the rosters in filterAlnsByRBH (the first user of a query/target pair).
    std::vector<size_t> qFoundEpoch;
    std::vector<size_t> dbFoundEpoch;
    size_t chainEpoch;
    std::vector<cluster_t> currClusters;
    std::set<cluster_t> &finalClusters;
    // Best tmScore per query/target chain, indexed by the chain's local id (position in its roster),
    // reused across pairs -- replaces the old chainKey-keyed std::map.
    std::vector<float> qBestScore;
    std::vector<float> dbBestScore;
    bool useKdTree;
    SuperpositionKdTree kdTree;

    bool earlyStopForMonomers() {
        if (minimumClusterSize >= MULTIPLE_CHAINED_COMPLEX)
            return finishDBSCAN();

        getSingleChainedCluster();
        return finishDBSCAN();
    }

    void getSingleChainedCluster() {
        finalClusters.clear();
        for (unsigned int alnIdx = 0; alnIdx < searchResult.alnVec.size(); alnIdx++ ) {
            neighbors = {alnIdx};
            finalClusters.insert(neighbors);
        }
    }

    bool runDBSCAN() {
        unsigned int neighborIdx;
        unsigned int neighborAlnIdx;
        foundEpoch.assign(searchResult.alnVec.size(), 0);
        curEpoch = 0;
        while (eps < maxDist) {
            initializeAlnLabels();
            for (size_t centerAlnIdx = 0; centerAlnIdx < searchResult.alnVec.size(); centerAlnIdx++) {
                ChainToChainAln &centerAln = searchResult.alnVec[centerAlnIdx];
                if (centerAln.label != 0)
                    continue;

                getNeighbors(centerAlnIdx, neighbors);
                if (neighbors.size() < MIN_PTS)
                    continue;

                centerAln.label = ++cLabel;
                curEpoch++;
                for (auto initNeighbor : neighbors) {
                    foundEpoch[initNeighbor] = curEpoch;
                }
                neighborIdx = 0;
                while (neighborIdx < neighbors.size()) {
                    neighborAlnIdx = neighbors[neighborIdx++];
                    if (centerAlnIdx == neighborAlnIdx)
                        continue;

                    ChainToChainAln &neighborAln = searchResult.alnVec[neighborAlnIdx];
                    neighborAln.label = cLabel;
                    getNeighbors(neighborAlnIdx, neighborsOfCurrNeighbor);
                    if (neighborsOfCurrNeighbor.size() < MIN_PTS)
                        continue;

                    for (auto neighbor : neighborsOfCurrNeighbor) {
                        if (foundEpoch[neighbor] != curEpoch) {
                            foundEpoch[neighbor] = curEpoch;
                            neighbors.emplace_back(neighbor);
                        }
                    }
                }
                if (neighbors.size() > maximumClusterSize || checkChainRedundancy())
                    getNearestNeighbors(centerAlnIdx);

                // too small cluster
                if (neighbors.size() < currMaxClusterSize)
                    continue;

                // new Biggest cluster
                if (neighbors.size() > currMaxClusterSize) {
                    currMaxClusterSize = neighbors.size();
                    currClusters.clear();
                }
                SORT_SERIAL(neighbors.begin(), neighbors.end());
                currClusters.emplace_back(neighbors);
            }

            if (!finalClusters.empty() && currClusters.empty())
                return finishDBSCAN();

            if (currMaxClusterSize < prevMaxClusterSize)
                return finishDBSCAN();

            if (currMaxClusterSize > prevMaxClusterSize) {
                finalClusters.clear();
                prevMaxClusterSize = currMaxClusterSize;
            }

            if (currMaxClusterSize >= minimumClusterSize)
                finalClusters.insert(currClusters.begin(), currClusters.end());

            if (currMaxClusterSize == maximumClusterSize && finalClusters.size() == maximumClusterNum)
                return finishDBSCAN();

            eps += learningRate;
        }

        if (minimumClusterSize < MULTIPLE_CHAINED_COMPLEX && prevMaxClusterSize < MULTIPLE_CHAINED_COMPLEX)
            getSingleChainedCluster();

        return finishDBSCAN();
    }

    // Lazy DBSCAN: no O(n^2) distance matrix is stored -- materializing N(N-1)/2 floats OOMs for very
    // large complexes (e.g. a 1356-mer would need terabytes). This single pass only determines the eps
    // range (min/max pairwise distance); individual distances are recomputed on demand in
    // getNeighbors/getNearestNeighbors via ChainToChainAln::getDistance (a cheap 12-dim vector distance).
    void initEpsRange() {
        for (size_t i = 0; i < searchResult.alnVec.size(); i++) {
            const ChainToChainAln &prevAln = searchResult.alnVec[i];
            for (size_t j = i + 1; j < searchResult.alnVec.size(); j++) {
                float dist = prevAln.getDistance(searchResult.alnVec[j]);
                maxDist = std::max(maxDist, dist);
                minDist = std::min(minDist, dist);
            }
        }
        eps = minDist;
    }

    void getNeighbors(size_t centerIdx, std::vector<unsigned int> &neighborVec) const {
        neighborVec.clear();
        neighborVec.emplace_back(centerIdx);
        if (useKdTree) {
            // tree returns the same set as the scan below; sort to reproduce ascending-index order
            kdTree.rangeQuery(centerIdx, eps, neighborVec);
            std::sort(neighborVec.begin() + 1, neighborVec.end());
            return;
        }
        const ChainToChainAln &centerAln = searchResult.alnVec[centerIdx];
        for (size_t neighborIdx = 0; neighborIdx < searchResult.alnVec.size(); neighborIdx++) {

            if (neighborIdx == centerIdx)
                continue;

            if (centerAln.getDistance(searchResult.alnVec[neighborIdx]) >= eps)
                continue;

            neighborVec.emplace_back(neighborIdx);
        }
    }

    void initializeAlnLabels() {
        for (auto &aln : searchResult.alnVec) {
            aln.label = INITIALIZED_LABEL;
        }
        cLabel = INITIALIZED_LABEL;
        currMaxClusterSize = 0;
        currClusters.clear();
    }

    bool checkChainRedundancy() {
        chainEpoch++;
        for (auto neighborIdx : neighbors) {
            const ChainToChainAln &aln = searchResult.alnVec[neighborIdx];
            if (qFoundEpoch[aln.qLocalId] == chainEpoch)
                return true;
            qFoundEpoch[aln.qLocalId] = chainEpoch;

            if (dbFoundEpoch[aln.dbLocalId] == chainEpoch)
                return true;
            dbFoundEpoch[aln.dbLocalId] = chainEpoch;
        }
        return false;
    }

    bool checkClusteringNecessity() {
        // Too few alns => do nothing and finish it
        if (searchResult.alnVec.size() < minimumClusterSize)
            return finishDBSCAN();
        // All alns as a cluster
        for (size_t alnIdx=0; alnIdx<searchResult.alnVec.size(); alnIdx++) {
            neighbors.emplace_back(alnIdx);
        }
        // Redundant chains => DBSCAN clustering
        if (checkChainRedundancy()) {
            neighbors.clear();
            return runDBSCAN();
        }
        // Already good => finish it without clustering
        finalClusters.insert(neighbors);
        return finishDBSCAN();
    }

    bool finishDBSCAN() {
        initializeAlnLabels();
        neighbors.clear();
        neighborsOfCurrNeighbor.clear();
        neighborsWithDist.clear();
        return !finalClusters.empty();
    }

    void filterAlnsByRBH() {
        // best tmScore per chain, indexed by dense local id (position in the roster) instead of a map
        qBestScore.assign(searchResult.qChainKeys.size(), FLT_MIN);
        dbBestScore.assign(searchResult.dbChainKeys.size(), FLT_MIN);
        // size (and reset) the "chain seen" epoch arrays for this query/target pair -- reused by
        // checkChainRedundancy/getNearestNeighbors, which only bump chainEpoch afterwards.
        qFoundEpoch.assign(searchResult.qChainKeys.size(), 0);
        dbFoundEpoch.assign(searchResult.dbChainKeys.size(), 0);
        chainEpoch = 0;
        for (auto &aln: searchResult.alnVec) {
            qBestScore[aln.qLocalId] = std::max(aln.tmScore, qBestScore[aln.qLocalId]);
            dbBestScore[aln.dbLocalId] = std::max(aln.tmScore, dbBestScore[aln.dbLocalId]);
        }
        // Single-pass in-place compaction (erase-remove idiom): drop alns below the RBH margin by
        // moving the kept ones forward, then truncate once. Equivalent to the previous per-element
        // vector::erase but O(n) moves instead of O(n^2) -- the old erase shifted the whole tail of a
        // coordinate-heavy vector on every drop, which is catastrophic for very large complexes. Kept
        // set and order are identical, so the output is byte-identical. qFoundNum/dbFoundNum count the
        // distinct kept query/target chains (what the old qFound/dbFound set sizes measured).
        chainEpoch++;
        size_t qFoundNum = 0;
        size_t dbFoundNum = 0;
        size_t writeIdx = 0;
        for (size_t readIdx = 0; readIdx < searchResult.alnVec.size(); readIdx++) {
            ChainToChainAln &aln = searchResult.alnVec[readIdx];
            if (aln.tmScore < std::max(qBestScore[aln.qLocalId], dbBestScore[aln.dbLocalId]) * TM_SCORE_MARGIN) {
                continue;
            }
            if (qFoundEpoch[aln.qLocalId] != chainEpoch) {
                qFoundEpoch[aln.qLocalId] = chainEpoch;
                qFoundNum++;
            }
            if (dbFoundEpoch[aln.dbLocalId] != chainEpoch) {
                dbFoundEpoch[aln.dbLocalId] = chainEpoch;
                dbFoundNum++;
            }
            if (writeIdx != readIdx) {
                searchResult.alnVec[writeIdx] = std::move(searchResult.alnVec[readIdx]);
            }
            writeIdx++;
        }
        searchResult.alnVec.resize(writeIdx);

        if (std::min(qFoundNum, dbFoundNum) < minimumClusterSize)
            searchResult.alnVec.clear();
    }

    void getNearestNeighbors(unsigned int centerIdx) {
        chainEpoch++;
        neighborsWithDist.clear();
        neighborsWithDist.emplace_back(centerIdx, 0.0);
        for (auto neighborIdx: neighbors) {
            if (neighborIdx == centerIdx)
                continue;
            neighborsWithDist.emplace_back(neighborIdx, searchResult.alnVec[centerIdx].getDistance(searchResult.alnVec[neighborIdx]));
        }
        SORT_SERIAL(neighborsWithDist.begin(), neighborsWithDist.end(), compareNeighborWithDist);
        neighbors.clear();
        for (auto neighborWithDist : neighborsWithDist) {
            const ChainToChainAln &aln = searchResult.alnVec[neighborWithDist.neighbor];
            if (qFoundEpoch[aln.qLocalId] == chainEpoch)
                break;
            qFoundEpoch[aln.qLocalId] = chainEpoch;
            if (dbFoundEpoch[aln.dbLocalId] == chainEpoch)
                break;
            dbFoundEpoch[aln.dbLocalId] = chainEpoch;
            neighbors.emplace_back(neighborWithDist.neighbor);
        }
    }
};

// One parsed alignment line, tagged with its target complex. The whole query complex's lines are
// streamed into a reused vector (mmseqs style: no std::map, no random access), then sorted so that a
// target complex's chains are contiguous. `line` is the raw (compressed) alignment record copied out
// of the DB, because the DBReader data pointer is only valid until the next getData on another key.
struct AlnLine {
    AlnLine(unsigned int dbComplexId, unsigned int qChainKey, unsigned int dbChainKey, unsigned int ordinal, unsigned int lineOffset)
        : dbComplexId(dbComplexId), qChainKey(qChainKey), dbChainKey(dbChainKey), ordinal(ordinal), lineOffset(lineOffset) {}
    unsigned int dbComplexId;
    unsigned int qChainKey;
    unsigned int dbChainKey;
    unsigned int ordinal;      // position in reading order; tiebreaker so the first row of a duplicate
                               // (qChainKey,dbChainKey) is the one the old std::map<...>.insert kept.
    unsigned int lineOffset;   // start of the null-terminated alignment row in ComplexScorer::lineArena
};

// Sort by (target complex, query chain, target chain, reading order). The primary key groups a target
// complex's lines into one contiguous run in ascending complex-id order, which keeps the later
// (unstable) assignment sort byte-identical to the previous full loop; the secondary keys make the
// order within a group deterministic. The ordinal tiebreaker guarantees that among duplicate
// (qChainKey,dbChainKey) rows the reading-first one comes first, so scoreDbComplexGroup can keep just
// it -- matching the old map, which ignored later rows for an already-seen (qChainKey,dbChainKey).
bool compareAlnLine(const AlnLine &first, const AlnLine &second) {
    if (first.dbComplexId != second.dbComplexId) {
        return first.dbComplexId < second.dbComplexId;
    }
    if (first.qChainKey != second.qChainKey) {
        return first.qChainKey < second.qChainKey;
    }
    if (first.dbChainKey != second.dbChainKey) {
        return first.dbChainKey < second.dbChainKey;
    }
    return first.ordinal < second.ordinal;
}

class ComplexScorer {
public:
    ComplexScorer(IndexReader *q3diDbr, IndexReader *t3diDbr, DBReader<unsigned int> &alnDbr, DBReader<unsigned int> *qCaDbr, DBReader<unsigned int> *tCaDbr, DBReader<unsigned int> *dbLookupReader, unsigned int thread_idx, float minAssignedChainsRatio, int monomerIncludeMode)
        : q3diDbr(q3diDbr), t3diDbr(t3diDbr), alnDbr(alnDbr), qCaDbr(qCaDbr), tCaDbr(tCaDbr), dbLookupReader(dbLookupReader), thread_idx(thread_idx), minAssignedChainsRatio(minAssignedChainsRatio), monomerIncludeMode(monomerIncludeMode)  {
        maxChainLen = std::max(q3diDbr->sequenceReader->getMaxSeqLen()+1, t3diDbr->sequenceReader->getMaxSeqLen()+1);
        maxResLen = maxChainLen * 2;
        tmAligner = new TMaligner(maxResLen, false, true, false);
    }
    ~ComplexScorer() {
        delete tmAligner;
    }

    // Stream every alignment line of the whole query complex into a reused vector (mmseqs style: one
    // getData per query chain, parse the record keys in place, no std::map, no random access), tagging
    // each with its target complex, then sort so each target complex is one contiguous run. Target
    // chains absent from the target sequence DB (e.g. cluster-search DBs keep only representatives) are
    // skipped here -- exactly the chains the old code kept in the map but never consumed. Sets qResLen.
    void buildAlnLines(std::vector<unsigned int> &qChainKeys, std::vector<AlnLine> &alnLines) {
        alnLines.clear();
        lineArena.clear();
        qResLen = getQueryResidueLength(qChainKeys);
        if (qResLen == 0) {
            return;
        }

        for (auto qChainKey: qChainKeys) {
            unsigned int qDbKey = alnDbr.getId(qChainKey);
            if (qDbKey == NOT_AVAILABLE_CHAIN_KEY) {
                continue;
            }
            char *data = alnDbr.getData(qDbKey, thread_idx);
            size_t dataSize = alnDbr.getDataSize();
            if (*data == '\0') {
                continue;
            }
            char dbKeyBuffer[256];
            char lineBuffer[1024];
            while (*data != '\0') {
                Util::parseKey(data, dbKeyBuffer);
                const auto dbChainKey = static_cast<unsigned int>(strtoul(dbKeyBuffer, NULL, 10));
                size_t lookupId = dbLookupReader->getLookupIdByKey(dbChainKey);
                // keep only target chains that resolve to a complex AND are present in the target
                // sequence DB (== the chains dbComplexIdToChainKeys buckets, i.e. the ones the old loop
                // actually processed); anything else never contributed to the output.
                if (lookupId != SIZE_MAX && t3diDbr->sequenceReader->getId(dbChainKey) != NOT_AVAILABLE_CHAIN_KEY) {
                    unsigned int dbComplexId = dbLookupReader->getLookupFileNumber(lookupId);
                    Util::getLine(data, dataSize, lineBuffer, 1024);
                    unsigned int lineOffset = (unsigned int) lineArena.size();
                    lineArena.insert(lineArena.end(), lineBuffer, lineBuffer + strlen(lineBuffer) + 1);
                    alnLines.emplace_back(dbComplexId, qChainKey, dbChainKey, (unsigned int)alnLines.size(), lineOffset);
                }
                data = Util::skipLine(data);
            } // while end
        } // for end
        SORT_SERIAL(alnLines.begin(), alnLines.end(), compareAlnLine);
    }

    // Score one target complex: the contiguous run alnLines[lo, hi) all share dbComplexId and is sorted
    // by (qChainKey, dbChainKey). Builds currAlns and fills searchResult, or leaves it empty (skip).
    void scoreDbComplexGroup(unsigned int qComplexId, unsigned int dbComplexId, std::vector<unsigned int> &qChainKeys, std::vector<unsigned int> &dbChainKeys, std::vector<AlnLine> &alnLines, size_t lo, size_t hi, SearchResult &searchResult) {
        unsigned int dbResLen = getDbResidueLength(dbChainKeys);
        if (dbResLen == 0) {
            return;
        }
        // Abandon rule (mirrors the old per-query getQueryAlnResult early-return): the pair is scored
        // only if EVERY query chain of the complex has at least one alignment to this target complex.
        // The run is sorted by qChainKey, so the number of distinct query chains present is the number
        // of qChainKey boundaries; if it does not cover the whole roster, skip.
        size_t distinctQ = 0;
        for (size_t i = lo; i < hi; i++) {
            if (i == lo || alnLines[i].qChainKey != alnLines[i - 1].qChainKey) {
                distinctQ++;
            }
        }
        if (distinctQ != qChainKeys.size()) {
            return;
        }

        Coordinate16 qCoords;
        Coordinate16 tCoords;
        currAlns.clear();
        // reset the per-target-complex CA decode cache
        tCaOffset.assign(dbChainKeys.size(), ~0u);
        tCaArena.clear();
        bool hasBacktrace = false;
        // for each query chain (its lines are contiguous in the sorted run)
        size_t i = lo;
        while (i < hi) {
            unsigned int qChainKey = alnLines[i].qChainKey;
            unsigned int qLocalId = (unsigned int)(std::lower_bound(qChainKeys.begin(), qChainKeys.end(), qChainKey) - qChainKeys.begin());
            Chain qChain(qComplexId, qChainKey);
            bool queryInited = false;
            bool havePrevDb = false;
            unsigned int prevDbKey = 0;
            while (i < hi && alnLines[i].qChainKey == qChainKey) {
                unsigned int dbChainKey = alnLines[i].dbChainKey;
                // Keep only the first row per (qChainKey, dbChainKey); duplicates are adjacent after the
                // ordinal-tiebroken sort. The old std::map<{q,db},line> likewise ignored later rows for an
                // already-seen (q,db) pair, so keeping just the first preserves byte-identical behaviour.
                if (havePrevDb && dbChainKey == prevDbKey) {
                    i++;
                    continue;
                }
                havePrevDb = true;
                prevDbKey = dbChainKey;
                Matcher::result_t dbAlnResult = Matcher::parseAlignmentRecord(&lineArena[alnLines[i].lineOffset]);
                // init the aligner once per query chain (qLen is the query chain length, identical for
                // all its lines), before skipping any empty-backtrace lines -- matches the old order.
                if (!queryInited) {
                    size_t qDbId = qCaDbr->getId(qChainKey);
                    char *qCaData = qCaDbr->getData(qDbId, thread_idx);
                    size_t qCaLength = qCaDbr->getEntryLen(qDbId);
                    unsigned int qLen = dbAlnResult.qLen;
                    float *queryCaData = qCoords.read(qCaData, qLen, qCaLength);
                    tmAligner->initQuery(queryCaData, &queryCaData[qLen], &queryCaData[qLen * 2], NULL, qLen);
                    queryInited = true;
                }
                if (dbAlnResult.backtrace.empty()) {
                    i++;
                    continue;
                }
                hasBacktrace = true;
                unsigned int dbLocalId = (unsigned int)(std::lower_bound(dbChainKeys.begin(), dbChainKeys.end(), dbChainKey) - dbChainKeys.begin());
                unsigned int dbLen = dbAlnResult.dbLen;
                // decode this target chain's CA only the first time it is seen in this target complex;
                // copy into the owned arena (Coordinate16 reuses its scratch buffer). Same floats as a
                // fresh decode, so byte-identical.
                if (tCaOffset[dbLocalId] == ~0u) {
                    size_t tCaId = tCaDbr->getId(dbChainKey);
                    float *decoded = tCoords.read(tCaDbr->getData(tCaId, thread_idx), dbLen, tCaDbr->getEntryLen(tCaId));
                    tCaOffset[dbLocalId] = (unsigned int) tCaArena.size();
                    tCaArena.insert(tCaArena.end(), decoded, decoded + 3 * (size_t)dbLen);
                }
                float *targetCaData = &tCaArena[tCaOffset[dbLocalId]];
                Chain dbChain(dbComplexId, dbChainKey);
                TMaligner::TMscoreResult tmResult = tmAligner->computeTMscore(targetCaData, &targetCaData[dbLen], &targetCaData[dbLen * 2], dbLen, dbAlnResult.qStartPos, dbAlnResult.dbStartPos, Matcher::uncompressAlignment(dbAlnResult.backtrace), dbAlnResult.qLen);
                currAlns.emplace_back(qChain, dbChain, dbAlnResult, tmResult);
                currAlns.back().qLocalId = qLocalId;
                currAlns.back().dbLocalId = dbLocalId;
                currAlns.back().lineIdx = (unsigned int)i;
                i++;
            }
        } // while end

        // When no alignment is found.
        if (currAlns.empty()) {
            return;
        }

        // When alignments have no backtrace
        if (!hasBacktrace) {
            Debug(Debug::ERROR) << "Backtraces are required. Please run search with '-a' option.\n";
            EXIT(EXIT_FAILURE);
        }

        SORT_SERIAL(currAlns.begin(), currAlns.end(), compareChainToChainAln);
        searchResult = SearchResult(qChainKeys, qResLen, dbChainKeys, dbResLen, currAlns);
        searchResult.standardize(monomerIncludeMode);
    }

    void getAssignments(SearchResult &searchResult, std::vector<Assignment> &assignments, const std::vector<AlnLine> &alnLines) {
        if (maxResLen < maxChainLen * std::min(searchResult.qChainKeys.size(),  searchResult.dbChainKeys.size())) {
            delete tmAligner;
            maxResLen = std::max(searchResult.qChainKeys.size(), searchResult.dbChainKeys.size()) * maxChainLen;
            tmAligner = new TMaligner(maxResLen, false, true, false);
        }
        finalClusters.clear();
        DBSCANCluster dbscanCluster(searchResult, finalClusters, minAssignedChainsRatio);
        if (!dbscanCluster.getAlnClusters()) {
            finalClusters.clear();
            return;
        }
        // Reused per-thread buffers for re-decoding the CA coordinates of final-cluster alns.
        Coordinate16 qCoords;
        Coordinate16 tCoords;
        for (auto &cluster : finalClusters) {
            Assignment assignment(searchResult.qResidueLen, searchResult.dbResidueLen);
            // Pre-size the accumulators to their exact final size (each aln contributes `matches` CA
            // coords to all six vectors and one entry to the per-chain vectors), so appendClusterAln
            // never reallocates. reserve() only touches capacity, so the result is byte-identical.
            unsigned int clusterMatches = 0;
            for (auto alnIdx : cluster) {
                clusterMatches += searchResult.alnVec[alnIdx].matches;
            }
            assignment.qCaXVec.reserve(clusterMatches);
            assignment.qCaYVec.reserve(clusterMatches);
            assignment.qCaZVec.reserve(clusterMatches);
            assignment.dbCaXVec.reserve(clusterMatches);
            assignment.dbCaYVec.reserve(clusterMatches);
            assignment.dbCaZVec.reserve(clusterMatches);
            assignment.matchLenVec.reserve(cluster.size());
            assignment.chainToChainResults.reserve(cluster.size());
            for (auto alnIdx : cluster) {
                appendClusterAln(assignment, searchResult.alnVec[alnIdx], alnLines, qCoords, tCoords);
            }
            // matches==0 ? skip
            if (!assignment.getTmScore(*tmAligner)) {
                continue;
            }
            assignment.updateResultToWriteLines();
            assignments.emplace_back(std::move(assignment));
        }
        finalClusters.clear();
    }

private:
    IndexReader *q3diDbr;
    IndexReader *t3diDbr;
    DBReader<unsigned int> &alnDbr;
    DBReader<unsigned int> *qCaDbr;
    DBReader<unsigned int> *tCaDbr;
    DBReader<unsigned int> *dbLookupReader;
    const unsigned int thread_idx;
    const float minAssignedChainsRatio;
    const int monomerIncludeMode;
    TMaligner *tmAligner;
    unsigned int maxChainLen;
    unsigned int maxResLen;
    std::vector<ChainToChainAln> currAlns;
    std::set<cluster_t> finalClusters;
    unsigned int qResLen;
    // Contiguous storage for the query complex's alignment row strings (null-terminated, back to back),
    // reused across complexes. AlnLine holds an offset into this instead of a per-row std::string --
    // avoids millions of small allocations and keeps the rows cache-local.
    std::vector<char> lineArena;
    // Per target-complex cache of decoded target CA coordinates (x|y|z, 3*dbLen floats per chain), so a
    // target chain that aligns to many query chains is decoded once per complex instead of once per
    // alignment. tCaOffset[dbLocalId] is the chain's offset in tCaArena, or ~0u if not yet decoded.
    std::vector<float> tCaArena;
    std::vector<unsigned int> tCaOffset;

    // Re-decode one final-cluster alignment's matched CA coordinates and append them (plus its match
    // count / per-chain result) to the Assignment. ChainToChainAln deliberately does not store the CA
    // vectors (memory), so they are reconstructed here for the few alns that reach a cluster: parse the
    // stored alignment line to recover qStart/dbStart/len/backtrace, re-read the chain CA from the same
    // CA DB, and replay the exact backtrace walk the ChainToChainAln ctor used -- the appended floats
    // are therefore bit-identical to the previously in-struct copy, and the per-cluster order is
    // preserved, so the output is byte-identical.
    void appendClusterAln(Assignment &assignment, const ChainToChainAln &aln, const std::vector<AlnLine> &alnLines, Coordinate16 &qCoords, Coordinate16 &tCoords) {
        Matcher::result_t r = Matcher::parseAlignmentRecord(&lineArena[alnLines[aln.lineIdx].lineOffset]);
        size_t qDbId = qCaDbr->getId(aln.qChain.chainKey);
        float *qCa = qCoords.read(qCaDbr->getData(qDbId, thread_idx), r.qLen, qCaDbr->getEntryLen(qDbId));
        size_t tDbId = tCaDbr->getId(aln.dbChain.chainKey);
        float *tCa = tCoords.read(tCaDbr->getData(tDbId, thread_idx), r.dbLen, tCaDbr->getEntryLen(tDbId));
        unsigned int qPos = r.qStartPos;
        unsigned int dbPos = r.dbStartPos;
        for (char cigar : r.backtrace) {
            if (cigar == 'M') {
                assignment.qCaXVec.emplace_back(qCa[qPos]);
                assignment.qCaYVec.emplace_back(qCa[r.qLen + qPos]);
                assignment.qCaZVec.emplace_back(qCa[r.qLen * 2 + qPos]);
                assignment.dbCaXVec.emplace_back(tCa[dbPos]);
                assignment.dbCaYVec.emplace_back(tCa[r.dbLen + dbPos]);
                assignment.dbCaZVec.emplace_back(tCa[r.dbLen * 2 + dbPos]);
                qPos++;
                dbPos++;
            } else if (cigar == 'I') {
                qPos++;
            } else if (cigar == 'D') {
                dbPos++;
            }
        }
        assignment.matches += aln.matches;
        assignment.matchLenVec.emplace_back(aln.matches);
        // Re-serialize the per-chain output line only now (final-cluster alns only) instead of storing
        // it in every ChainToChainAln. resultToBuffer(parseAlignmentRecord(line)) reproduces exactly the
        // string the old in-struct resultToWrite held, so the output is byte-identical.
        char buffer[4096];
        size_t len = Matcher::resultToBuffer(buffer, r, true, true, false);
        assignment.chainToChainResults.emplace_back(aln.qChain.chainKey, std::string(buffer, len - 1));
    }

    unsigned int getQueryResidueLength(std::vector<unsigned int> &qChainKeys) const {
        unsigned int qResidueLen = 0;
        size_t qDbId;
        for (auto qChainKey: qChainKeys) {
            qDbId = q3diDbr->sequenceReader->getId(qChainKey);
            // Not accessible
            if (qDbId == NOT_AVAILABLE_CHAIN_KEY)
                return 0;

            qResidueLen += q3diDbr->sequenceReader->getSeqLen(qDbId);
        }
        return qResidueLen;
    }

    unsigned int getDbResidueLength(const std::vector<unsigned int> &dbChainKeys) const {
        unsigned int dbResidueLen = 0;
        size_t tDbId;
        for (auto dbChainKey: dbChainKeys) {
            tDbId = t3diDbr->sequenceReader->getId(dbChainKey);
            // Not accessible
            if (tDbId == NOT_AVAILABLE_CHAIN_KEY)
                return 0;

            dbResidueLen += t3diDbr->sequenceReader->getSeqLen(tDbId);
        }
        return dbResidueLen;
    }
};

class ComplexFilter {
public:
    ComplexFilter(
        const std::vector<unsigned int> &qChainKeys,
        IndexReader *qDbr,
        DBReader<unsigned int> *qStructDbr,
        LocalParameters &par,
        unsigned int thread_idx
    )
        : qChainKeys_(qChainKeys),
          qDbr_(qDbr),
          qStructDbr_(qStructDbr),
          par_(par),
          thread_idx(thread_idx)
    {
        qInterfaceVec_.resize(qChainKeys_.size());
    }

    void computeInterfaceRegion()
    {
        float threshold = INTERFACE_THRESHOLD;
        float d2 = threshold * threshold;
        Coordinate16 coords, coords2;

        for (size_t chainIdx = 0; chainIdx < qChainKeys_.size(); chainIdx++) {
            unsigned int chainKey = qChainKeys_[chainIdx];
            qChainKeyTochainIdx_[chainKey] = chainIdx;

            unsigned int chainDbId = qDbr_->sequenceReader->getId(chainKey);
            char *cadata = qStructDbr_->getData(chainDbId, thread_idx);
            size_t caLength = qStructDbr_->getEntryLen(chainDbId);
            size_t chainLen = qDbr_->sequenceReader->getSeqLen(chainDbId);
            float *chainData = coords.read(cadata, chainLen, caLength);

            for (size_t chainIdx2 = 0; chainIdx2 < qChainKeys_.size(); chainIdx2++) {
                if (chainIdx == chainIdx2) continue;

                unsigned int chainKey2 = qChainKeys_[chainIdx2];
                unsigned int chainDbId2 = qDbr_->sequenceReader->getId(chainKey2);
                char *cadata2 = qStructDbr_->getData(chainDbId2, thread_idx);
                size_t caLength2 = qStructDbr_->getEntryLen(chainDbId2);
                size_t chainLen2 = qDbr_->sequenceReader->getSeqLen(chainDbId2);
                float *chainData2 = coords2.read(cadata2, chainLen2, caLength2);

                for (size_t res1 = 0; res1 < chainLen; res1++) {
                    for (size_t res2 = 0; res2 < chainLen2; res2++) {
                        float dist = BasicFunction::dist(
                            chainData[res1], chainData[chainLen + res1], chainData[2 * chainLen + res1],
                            chainData2[res2], chainData2[chainLen2 + res2], chainData2[2 * chainLen2 + res2]
                        );
                        if (dist < d2) {
                            qInterfaceVec_[chainIdx].push_back(res1);
                            break;
                        }
                    }
                }
            }
        }
    }
    void filterAssignment(
        unsigned int assignmentId,
        Assignment &assignment,
        std::map<unsigned int, std::pair<Assignment, unsigned int>> &tCompBestAssignment,
        DBReader<unsigned int> &dbLookupReader,
        std::vector<std::vector<unsigned int>> &dbComplexIdToChainKeys
    ) 
    {
        unsigned int tComplexId;
        unsigned int qalnlen = 0;
        unsigned int talnlen = 0;
        unsigned int adjustAlnLen;
        assignment.assignmentId = assignmentId;
        assignment.filterResult.clear();
        size_t alnChainNum = assignment.chainToChainResults.size();

        std::vector<int> qChainLengths(alnChainNum);
        std::vector<int> dbChainLengths(alnChainNum);
        std::vector<std::string> backtraceVec(alnChainNum);
        std::vector<int> qStartPosVec(alnChainNum);
        std::vector<int> dbStartPosVec(alnChainNum);
        std::vector<float> qTMscores(alnChainNum);
        std::vector<float> dbTMscores(alnChainNum);
        std::vector<size_t> qidxToqChainKey(alnChainNum);
        std::string uString, tString;

        if (alnChainNum < (size_t)par_.minAlignedChains)
            return;

        unsigned int dbChainNum = 0;

        for (size_t i = 0; i < alnChainNum; i++) {

            const unsigned int &qKey = assignment.chainToChainResults[i].first;
            std::string resultToWrite = assignment.chainToChainResults[i].second+"\t"+assignment.assignmentResult+"\t"+std::to_string(assignment.assignmentId);
            const char* data = resultToWrite.c_str();

            Matcher::result_t res;
            ComplexDataHandler retComplex = parseScoreComplexResult(data, res);

            unsigned int tChainKey = res.dbKey;
            size_t lookupId = dbLookupReader.getLookupIdByKey(tChainKey);
            if (lookupId == SIZE_MAX) {
                return;
            }
            tComplexId = dbLookupReader.getLookupFileNumber(lookupId);
            if (tComplexId >= dbComplexIdToChainKeys.size()) {
                return;
            }

            qalnlen += abs(res.qEndPos - res.qStartPos) + 1;
            talnlen += abs(res.dbEndPos - res.dbStartPos) + 1;

            qChainLengths[i] = res.qLen;
            dbChainLengths[i] = res.dbLen;
            backtraceVec[i] = res.backtrace;
            qStartPosVec[i] = res.qStartPos;
            dbStartPosVec[i] = res.dbStartPos;

            uString = retComplex.uString;
            tString = retComplex.tString;

            dbChainNum = dbComplexIdToChainKeys[tComplexId].size();
            qidxToqChainKey[i] = qKey;
        }

        // check if multimer tm matches the threshold parameter
        if (par_.covMode == Parameters::COV_MODE_BIDIRECTIONAL && (assignment.qTmScore < par_.filtMultTmThr || assignment.dbTmScore < par_.filtMultTmThr)) {
            return;
        } else if (par_.covMode == Parameters::COV_MODE_TARGET && assignment.dbTmScore < par_.filtMultTmThr) {
            return;
        } else if (par_.covMode == Parameters::COV_MODE_QUERY && assignment.qTmScore < par_.filtMultTmThr) {
            return;
        }

        // check if multimer coverage matches the threshold parameter
        float qCov = static_cast<float>(qalnlen)/static_cast<float>(assignment.qResidueLength);
        float tCov = static_cast<float>(talnlen)/static_cast<float>(assignment.dbResidueLength);
        if (par_.covMode == Parameters::COV_MODE_BIDIRECTIONAL) {
            adjustAlnLen = (qCov+tCov)/2;
        } else if (par_.covMode == Parameters::COV_MODE_TARGET) {
            adjustAlnLen = tCov;
        } else if (par_.covMode == Parameters::COV_MODE_QUERY) {
            adjustAlnLen = qCov;
        }
        if(! Util::hasCoverage(par_.covThr, par_.covMode, qCov, tCov)){
            return;
        }
        // chain by chain tmscore
        unsigned int matchLen = 0;
        for (unsigned int qidx = 0; qidx < alnChainNum; qidx++) {
            unsigned int qmatchLen = assignment.matchLenVec[qidx];
            Coordinates tmt(qmatchLen), tchain(qmatchLen);
            for (unsigned int residx = matchLen; residx< qmatchLen + matchLen; residx++) {
                tchain.x[residx - matchLen] = assignment.dbCaXVec[residx];
                tchain.y[residx - matchLen] = assignment.dbCaYVec[residx];
                tchain.z[residx - matchLen] = assignment.dbCaZVec[residx];
            }
            // fill t and u
            float t[3];
            float u[3][3];
            std::string tmp;
            int ti = 0;
            const int tlen = static_cast<int>(tString.size());
            for (int k=0; k<tlen; k++) {
                if (k ==tlen-1) {
                    t[ti] = std::stof(tmp);
                } else if (tString[k] == ',') {
                    t[ti] = std::stof(tmp);
                    tmp.clear();
                    ti++;
                } else {
                    tmp.push_back(tString[k]);
                }
            }
            std::string tmp2;
            int ui = 0;
            int uj = 0;
            const int ulen = static_cast<int>(uString.size());
            for (int k=0; k < ulen; k++) {
                if (k==ulen-1) {
                    u[ui][uj] = std::stof(tmp2);
                } else if (uString[k] == ',') {
                    u[ui][uj] = std::stof(tmp2);
                    tmp2.clear();
                    uj++;
                } else {
                    tmp2.push_back(uString[k]);
                }
                if (uj == 3) {
                    ui++;
                    uj = 0;
                }
            }
            // based on t and u, calculate chain tm
            BasicFunction::do_rotation(tchain, tmt, qmatchLen, t, u);
            float d0 = 1.24*(cbrt(dbChainLengths[qidx]-15)) -1.8;
            float d02 = d0*d0;
            float tmScore = 0;
            for (unsigned int ci=0; ci<qmatchLen; ci++) {
                float xa_x = assignment.qCaXVec[matchLen + ci];
                float xa_y = assignment.qCaYVec[matchLen + ci];
                float xa_z = assignment.qCaZVec[matchLen + ci];
                float ya_x = tmt.x[ci];
                float ya_y = tmt.y[ci];
                float ya_z = tmt.z[ci];
                float di = BasicFunction::dist(xa_x, xa_y, xa_z, ya_x, ya_y, ya_z);
                float oneDividedDist = 1/(1+di/d02);
                tmScore += oneDividedDist;
            }
            float qTmScore = tmScore / qChainLengths[qidx];
            float dbTmScore = tmScore / dbChainLengths[qidx];
            qTMscores[qidx] = qTmScore;
            dbTMscores[qidx] = dbTmScore;
            matchLen+=qmatchLen;
        }
        if (par_.filtChainTmThr > 0) {
            // check chain-tm-threshold 
            // if cov-mode 0, every chain should be aligned, every tm scores should be higher than the threshold
            int chainpassNum = 0;
            if (par_.covMode == Parameters::COV_MODE_BIDIRECTIONAL) {
                if (dbChainNum != qChainKeys_.size() || dbChainNum != alnChainNum) {
                    return;
                }
                for (size_t qidx = 0; qidx < alnChainNum; qidx++) {
                    if(qTMscores[qidx] < par_.filtChainTmThr) {
                        return;
                    }
                    if(dbTMscores[qidx] < par_.filtChainTmThr) {
                        return;
                    }
                }
            // if cov-mode 1 or 2, min-aligned-chains should have tm higher than the threshold.
            } else if (par_.covMode == Parameters::COV_MODE_TARGET) {
                for (size_t qidx = 0; qidx < alnChainNum; qidx++) {
                    if (dbTMscores[qidx]>= par_.filtChainTmThr) {
                    chainpassNum++; 
                    }
                }
                if(chainpassNum < par_.minAlignedChains) {
                    return;
                }
            } else if (par_.covMode == Parameters::COV_MODE_QUERY) {
                for (size_t qidx = 0; qidx < alnChainNum; qidx++) {
                    if (qTMscores[qidx]>= par_.filtChainTmThr) {
                    chainpassNum++; 
                    }
                }
                if(chainpassNum < par_.minAlignedChains) {
                    return;
                }
            }
        }
        
        // interface-lddt, only check if aligned chains > 1
        // if the interface lddt parameter is set and aligned chain num ==1, then this assignment doesn't pass
        // if the interface lddt parameter isn't set, interfacelddt is printed out as 0 if aligned chain num ==1
        // if aligned chain num > 1, always calculate}
        float interfaceLddt = 0;
        if (alnChainNum == 1 && par_.filtInterfaceLddtThr > 0) {
            return;
        } else if (alnChainNum > 1) {
            std::vector<int> qcaToResidue(assignment.qResidueLength, -1);
            std::vector<int> dbcaToResidue(assignment.dbResidueLength, -1);
            unsigned int numIncrease = 0;
            unsigned int alnChainSum = 0;
            for (size_t qidx = 0; qidx < alnChainNum; qidx++) {
                std::string backtraceString = Matcher::uncompressAlignment(backtraceVec[qidx]);
                unsigned int qStart = qStartPosVec[qidx];
                unsigned int dbStart = dbStartPosVec[qidx];
                unsigned int qPos = 0;
                unsigned int dbPos = 0;
                for (size_t btPos = 0; btPos < backtraceString.size(); btPos++) {
                    if (backtraceString[btPos] == 'M') {
                        qcaToResidue[alnChainSum + qStart + qPos] = numIncrease;
                        dbcaToResidue[alnChainSum + dbStart + dbPos] = numIncrease++;
                        qPos++;
                        dbPos++;
                    } else if (backtraceString[btPos] == 'I') {
                        qPos++;
                    } else if (backtraceString[btPos] == 'D') {
                        dbPos++;
                    }
                }
                alnChainSum += assignment.matchLenVec[qidx];
            }
            std::vector<float> qIntVecX, qIntVecY, qIntVecZ, dbIntVecX, dbIntVecY, dbIntVecZ;
            alnChainSum = 0;
            unsigned int wholeIntLen = 0;
            for (size_t qidx = 0; qidx < alnChainNum; qidx++) {
                std::string backtraceString = Matcher::uncompressAlignment(backtraceVec[qidx]);
                std::vector<int> qchainIsMatch(qChainLengths[qidx], -1);
                std::vector<bool> qchainIsVisited(qChainLengths[qidx], 0);
                unsigned int qStart = qStartPosVec[qidx];
                unsigned int dbStart = dbStartPosVec[qidx];
                unsigned int qPos = 0;
                unsigned int dbPos = 0;
                for (size_t btPos = 0; btPos < backtraceString.size(); btPos++) {
                    if (backtraceString[btPos] == 'M') {
                        qchainIsMatch[qStart + qPos] = dbStart + dbPos;
                        qPos++;
                        dbPos++;
                    } else if (backtraceString[btPos] == 'I') {
                        qPos++;
                    } else if (backtraceString[btPos] == 'D') {
                        dbPos++;
                    }
                }
                std::vector<unsigned int>& qChainInterfaceVec =  qInterfaceVec_[qChainKeyTochainIdx_[qidxToqChainKey[qidx]]];
                wholeIntLen += qChainInterfaceVec.size();
                for (size_t resIdIdx = 0; resIdIdx< qChainInterfaceVec.size(); resIdIdx++) {
                    if(qchainIsMatch[qChainInterfaceVec[resIdIdx]] > -1 && qchainIsVisited[qChainInterfaceVec[resIdIdx]] == 0) {
                        qIntVecX.push_back(assignment.qCaXVec[qcaToResidue[alnChainSum + qChainInterfaceVec[resIdIdx]]]);
                        qIntVecY.push_back(assignment.qCaYVec[qcaToResidue[alnChainSum + qChainInterfaceVec[resIdIdx]]]);
                        qIntVecZ.push_back(assignment.qCaZVec[qcaToResidue[alnChainSum + qChainInterfaceVec[resIdIdx]]]);
                        dbIntVecX.push_back(assignment.dbCaXVec[dbcaToResidue[alnChainSum + qchainIsMatch[qChainInterfaceVec[resIdIdx]]]]);
                        dbIntVecY.push_back(assignment.dbCaYVec[dbcaToResidue[alnChainSum + qchainIsMatch[qChainInterfaceVec[resIdIdx]]]]);
                        dbIntVecZ.push_back(assignment.dbCaZVec[dbcaToResidue[alnChainSum + qchainIsMatch[qChainInterfaceVec[resIdIdx]]]]);
                        qchainIsVisited[qChainInterfaceVec[resIdIdx]] = 1;
                    }
                }
                alnChainSum += assignment.matchLenVec[qidx];
            }
            unsigned int intAlnLen = dbIntVecX.size();
            if (intAlnLen > 0) {
                std::string alnbt(intAlnLen, 'M');
                LDDTCalculator lddtcalculator(intAlnLen+1, intAlnLen+1);
                lddtcalculator.initQuery(intAlnLen, &qIntVecX[0], &qIntVecY[0], &qIntVecZ[0]);
                LDDTCalculator::LDDTScoreResult lddtres = lddtcalculator.computeLDDTScore(intAlnLen, 0, 0, alnbt, &dbIntVecX[0], &dbIntVecY[0], &dbIntVecZ[0]);
                interfaceLddt = lddtres.avgLddtScore * lddtres.scoreLength / wholeIntLen;

                if(interfaceLddt < par_.filtInterfaceLddtThr) {
                    return;
                }
            } else if (par_.filtInterfaceLddtThr > 0) {
                return;
            }
        }
        // write down qcov, tcov, qchaintms, tchaintms, interface lddt if everything passed
        std::string result;
        result.append(SSTR(qCov));
        result.append("\t");
        result.append(SSTR(tCov));
        result.append("\t");
        for (unsigned int i = 0; i < qTMscores.size(); i++) {
            result.append(SSTR(qTMscores[i]));
            if(i < qTMscores.size() - 1) {
                result.append(",");
            }
        }
        result.append("\t");
        for (unsigned int i = 0; i < dbTMscores.size(); i++) {
            result.append(SSTR(dbTMscores[i]));
            if(i < dbTMscores.size() - 1) {
                result.append(",");
            }
        }
        result.append("\t");
        result.append(SSTR(interfaceLddt));
        assignment.filterResult = result;

        auto it = tCompBestAssignment.find(tComplexId);
        if (it == tCompBestAssignment.end() || adjustAlnLen > it->second.second)
            tCompBestAssignment[tComplexId] = {assignment, adjustAlnLen};
    }

private:
    const std::vector<unsigned int> &qChainKeys_;
    IndexReader *qDbr_;
    DBReader<unsigned int> *qStructDbr_;
    LocalParameters &par_;

    std::vector<std::vector<unsigned int>> qInterfaceVec_;
    std::map<unsigned int, unsigned int> qChainKeyTochainIdx_;
    unsigned int thread_idx;
};


// Use the DBReader lookup table as the source of complex membership instead of reparsing .lookup.
static void buildComplexIdToChainKeys(
        DBReader<unsigned int> &lookupReader,
        IndexReader *seqDbr,
        std::vector<std::vector<unsigned int>> &complexIdToChainKeys,
        std::vector<unsigned int> *complexIds = NULL
) {
    DBReader<unsigned int>::LookupEntry *lookup = lookupReader.getLookup();
    const size_t lookupSize = lookupReader.getLookupSize();
    unsigned int maxComplexId = 0;
    for (size_t i = 0; i < lookupSize; i++) {
        maxComplexId = std::max(maxComplexId, lookup[i].fileNumber);
    }
    complexIdToChainKeys.assign(maxComplexId + 1, std::vector<unsigned int>());
    std::vector<bool> visited;
    if (complexIds != NULL) {
        complexIds->clear();
        visited.assign(maxComplexId + 1, false);
    }
    for (size_t i = 0; i < lookupSize; i++) {
        const unsigned int chainKey = lookup[i].id;
        if (seqDbr->sequenceReader->getId(chainKey) == NOT_AVAILABLE_CHAIN_KEY) {
            continue;
        }
        const unsigned int complexId = lookup[i].fileNumber;
        if (complexIds != NULL && visited[complexId] == false) {
            complexIds->emplace_back(complexId);
            visited[complexId] = true;
        }
        complexIdToChainKeys[complexId].emplace_back(chainKey);
    }
    if (complexIds != NULL) {
        SORT_SERIAL(complexIds->begin(), complexIds->end());
    }
}

int scoremultimer(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    int dbType = alnDbr.getDbtype();
    uint16_t extended = DBReader<unsigned int>::getExtendedDbtype(dbType);
    bool needSrc = false;
    if (extended & Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC) {
        needSrc = true;
        dbType = DBReader<unsigned int>::setExtendedDbtype(dbType, Parameters::DBTYPE_EXTENDED_INDEX_NEED_SRC);
    }
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed, dbType);
    resultWriter.open();

    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);

    std::string t3DiDbrName =  StructureUtil::getIndexWithSuffix(par.db2, "_ss");
    bool is3DiIdx = Parameters::isEqualDbtype(FileUtil::parseDbType(t3DiDbrName.c_str()), Parameters::DBTYPE_INDEX_DB);
    IndexReader* t3DiDbr = new IndexReader(
        is3DiIdx ? t3DiDbrName : par.db2,
        par.threads,
        needSrc ? IndexReader::SRC_SEQUENCES : IndexReader::SEQUENCES,
        touch ? IndexReader::PRELOAD_INDEX : 0,
        DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
        needSrc ? "_seq_ss" : "_ss"
    );
    DBReader<unsigned int>* tCaDbr = new DBReader<unsigned int>(
        needSrc? (par.db2 + "_seq_ca").c_str() : (par.db2 + "_ca").c_str(), 
        needSrc? (par.db2 + "_seq_ca.index").c_str() : (par.db2 + "_ca.index").c_str(),
        par.threads,
        DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA
    );
    tCaDbr->open(DBReader<unsigned int>::NOSORT);

    IndexReader* q3DiDbr = NULL;
    DBReader<unsigned int> *qCaDbr = NULL;
    bool sameDB = false;
    if (par.db1 == par.db2) {
        sameDB = true;
        q3DiDbr = t3DiDbr;
        qCaDbr = tCaDbr;
    } else {
        q3DiDbr = new IndexReader(
                StructureUtil::getIndexWithSuffix(par.db1, "_ss"),
                par.threads, IndexReader::SEQUENCES,
                touch ? IndexReader::PRELOAD_INDEX : 0,
                DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA
        );
        qCaDbr = new DBReader<unsigned int>((par.db1 + "_ca").c_str(), (par.db1 + "_ca.index").c_str(),
        par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        qCaDbr->open(DBReader<unsigned int>::NOSORT);
    }

    const float minAssignedChainsRatio = par.minAssignedChainsThreshold > MAX_ASSIGNED_CHAIN_RATIO ? MAX_ASSIGNED_CHAIN_RATIO : par.minAssignedChainsThreshold;
    int monomerIncludeMode = par.monomerIncludeMode;
    std::vector<unsigned int> qComplexIndices;
    std::vector<std::vector<unsigned int>> qComplexIdToChainKeys;
    std::vector<std::vector<unsigned int>> dbComplexIdToChainKeysStorage;

    DBReader<unsigned int> qLookupReader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_LOOKUP);
    qLookupReader.open(DBReader<unsigned int>::NOSORT);
    buildComplexIdToChainKeys(qLookupReader, q3DiDbr, qComplexIdToChainKeys, &qComplexIndices);

    DBReader<unsigned int> *dbLookupReader = &qLookupReader;
    std::vector<std::vector<unsigned int>> &dbComplexIdToChainKeys = sameDB ? qComplexIdToChainKeys : dbComplexIdToChainKeysStorage;
    if (!sameDB) {
        dbLookupReader = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_LOOKUP);
        dbLookupReader->open(DBReader<unsigned int>::NOSORT);
        buildComplexIdToChainKeys(*dbLookupReader, t3DiDbr, dbComplexIdToChainKeysStorage);
    }
    Debug::Progress progress(qComplexIndices.size());

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
        Coordinate16 qcoords;
        Coordinate16 tcoords;
        Matcher::result_t res;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        resultToWrite_t resultToWrite;
        resultToWrite_t currentResultToWrite;
        std::vector<AlnLine> alnLines;
        SearchResult searchResult;
        std::vector<Assignment> assignments;
        std::map<unsigned int, std::pair<Assignment, unsigned int>> tCompBestAssignment;
        ComplexScorer complexScorer(q3DiDbr, t3DiDbr, alnDbr, qCaDbr, tCaDbr, dbLookupReader, thread_idx, minAssignedChainsRatio, monomerIncludeMode);
#pragma omp for schedule(dynamic, 1)
        // for each q complex
        for (size_t qCompIdx = 0; qCompIdx < qComplexIndices.size(); qCompIdx++) {
            unsigned int qComplexId = qComplexIndices[qCompIdx];
            std::vector<unsigned int> &qChainKeys = qComplexIdToChainKeys[qComplexId];
            if (monomerIncludeMode == SKIP_MONOMERS && qChainKeys.size() < MULTIPLE_CHAINED_COMPLEX) {
                progress.updateProgress();
                continue;
            }
            // read the whole query complex's alignments once, streamed into a vector grouped (sorted)
            // by target complex. The hit target complexes fall out as the distinct dbComplexId runs.
            complexScorer.buildAlnLines(qChainKeys, alnLines);
            if (alnLines.empty()) {
                progress.updateProgress();
                continue;
            }
            // each contiguous run of equal dbComplexId is one target complex, in ascending order
            size_t runStart = 0;
            while (runStart < alnLines.size()) {
                unsigned int dbComplexId = alnLines[runStart].dbComplexId;
                size_t runEnd = runStart;
                while (runEnd < alnLines.size() && alnLines[runEnd].dbComplexId == dbComplexId) {
                    runEnd++;
                }
                std::vector<unsigned int> &dbChainKeys = dbComplexIdToChainKeys[dbComplexId];
                if (!dbChainKeys.empty()) {
                    complexScorer.scoreDbComplexGroup(qComplexId, dbComplexId, qChainKeys, dbChainKeys, alnLines, runStart, runEnd, searchResult);
                    if (!searchResult.alnVec.empty()) {
                        complexScorer.getAssignments(searchResult, assignments, alnLines);
                        searchResult.clear();
                    }
                }
                runStart = runEnd;
            }
            SORT_SERIAL(assignments.begin(), assignments.end(), compareAssignment);
            // Filter when multimercluster OR filtering paramters are set by user.
            if (par.filtMultTmThr + par.filtChainTmThr + par.filtInterfaceLddtThr > 0 ){
                // for each query chain key
                ComplexFilter filter(qChainKeys, q3DiDbr, qCaDbr, par, thread_idx);
                filter.computeInterfaceRegion();
                // for each assignment, filter
                for (unsigned int assignmentId = 0; assignmentId < assignments.size(); assignmentId++){
                    filter.filterAssignment(assignmentId, assignments[assignmentId], tCompBestAssignment, *dbLookupReader, dbComplexIdToChainKeys);
                }
                for (size_t qChainKeyIdx = 0; qChainKeyIdx < qChainKeys.size(); qChainKeyIdx++) {
                    resultToWrite.clear();
                    unsigned int & qKey = qChainKeys[qChainKeyIdx];
                    for(const auto &pair : tCompBestAssignment) {
                        const Assignment &assignment = pair.second.first;
                        currentResultToWrite.clear();
                        assignment.getChainToChainResult(qKey, currentResultToWrite);
                        if (currentResultToWrite.empty()) {
                            continue;
                        }
                        getResult(resultToWrite, currentResultToWrite, assignment);
                    }   
                    resultWriter.writeData(resultToWrite.c_str(), resultToWrite.length(), qKey, thread_idx);
                }
            } else {
                // for each assignment
                for (unsigned int assignmentId = 0; assignmentId < assignments.size(); assignmentId++){
                    Assignment &assignment = assignments[assignmentId];
                    assignment.assignmentId = assignmentId;
                }
                for (size_t qChainKeyIdx = 0; qChainKeyIdx < qChainKeys.size(); qChainKeyIdx++) {
                    resultToWrite.clear();
                    unsigned int & qKey = qChainKeys[qChainKeyIdx];
                    for (unsigned int assignmentId = 0; assignmentId < assignments.size(); assignmentId++){
                        Assignment &assignment = assignments[assignmentId];
                        currentResultToWrite.clear();
                        assignment.getChainToChainResult(qKey, currentResultToWrite);
                        if (currentResultToWrite.empty()) {
                            continue;
                        }
                        getResult(resultToWrite, currentResultToWrite, assignment);
                    }   
                    resultWriter.writeData(resultToWrite.c_str(), resultToWrite.length(), qKey, thread_idx);
                }
            }

            alnLines.clear();
            assignments.clear();
            currentResultToWrite.clear();
            resultToWrite.clear();
            tCompBestAssignment.clear();
            progress.updateProgress();
        }
    }
    alnDbr.close();
    delete t3DiDbr;
    tCaDbr->close();
    delete tCaDbr;
    if (!sameDB) {
        delete q3DiDbr;
        qCaDbr->close();
        delete qCaDbr;
    }
    qLookupReader.close();
    if (!sameDB) {
        dbLookupReader->close();
        delete dbLookupReader;
    }
    resultWriter.close(false);
    return EXIT_SUCCESS;
}
