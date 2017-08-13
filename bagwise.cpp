#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <random>
#include <cmath>
#include <limits>
#include <algorithm>
#include <utility>
#include <map>
#include <cmath>

#include "GaussianProcess.h"

 
struct Node {
    std::string word;
    std::map<int, int> child;
    Node(std::string w) : word(w), child() {
    }
};
 
class BKTree {
    std::vector<Node> tree;

    int min(int a, int b, int c);
    int editDistance(const std::string& a, const std::string& b);
    void add(Node& root, const std::string& str);
    std::vector <std::string> getSimilarWords(Node& root, const std::string& s, const int d);
public:
    void add(const std::string& str);
    std::vector <std::string> getSimilarWords(const std::string& s, const int d);
};
// stores the root Node
 
 
int BKTree::min(int a, int b, int c) {
    return std::min(a, std::min(b, c));
}
 
int BKTree::editDistance(const std::string& a, const std::string& b) {
    int m = a.length(), n = b.length();
    int dp[m+1][n+1];
 
    // filling base cases
    for (int i=0; i<=m; i++)
        dp[i][0] = i;
    for (int j=0; j<=n; j++)
        dp[0][j] = j;
 
    // populating matrix using dp-approach
    for (int i=1; i<=m; i++) {
        for (int j=1; j<=n; j++) {
            if (a[i-1] != b[j-1]) {
                dp[i][j] = min( 1 + dp[i-1][j],  // deletion
                                1 + dp[i][j-1],  // insertion
                                1 + dp[i-1][j-1] // replacement
                              );
            } else {
                dp[i][j] = dp[i-1][j-1];
            }
        }
    }
    return dp[m][n];
}
 
// adds curr Node to the tree
void BKTree::add(const std::string& str) {
    if(tree.size() == 0) {
        tree.push_back(Node(str));
    } else {
        add(tree[0], str);
    }
}
        
void BKTree::add(Node& root, const std::string& str) {
    int dist = editDistance(root.word, str);

    if (root.child.count(dist) == 0) {
        root.child[dist] = tree.size();
        tree.push_back(Node(str));
    } else {
        add(tree[root.child[dist]], str);
    }
}
 
std::vector <std::string> BKTree::getSimilarWords(const std::string& s, const int d) {
    if(tree.size() == 0) return std::vector<std::string>();
    else return getSimilarWords(tree[0], s, d);

}

std::vector <std::string> BKTree::getSimilarWords(Node& root, const std::string& s, const int d) {
    std::vector <std::string> ret;

    int dist = editDistance(root.word, s);
 
    if (dist <= d) ret.push_back(root.word);
 
    int start = dist - d;
    if (start < 0) start = 1;
 
    for (;start < dist+d; start++) {
        if(root.child.count(start) > 0) {
            std::vector <std::string> tmp =
             getSimilarWords(tree[root.child[start]], s, d);
            for (auto i : tmp)
                ret.push_back(i);
        }
    }
    return ret;
}

double CDF(double x) {
    return 0.5 + 0.5 * erf(x * M_SQRT1_2);
}

std::string random_triplet(int i) {
    char dna[] = {'A', 'T', 'G', 'C'};
    std::string triplet = "";
    triplet += dna[(i>>2)&3];
    triplet += dna[(i>>0)&3];
    triplet += dna[(i>>4)&3];
    return triplet;
}
std::map<char, std::vector<char>> nucleotide;
std::map<char, std::array<double, 4>> whatisgood;
std::map<std::string, char> codon;
std::map<char, std::string> r_codon;

std::vector<std::string> amino2base(std::string amino) {
    std::vector<std::string> ret = {""};
    for(auto c : amino) {
        std::vector<std::string> base;
        for(auto it=codon.begin(); it!=codon.end(); ++it) {
            if(it->second == c) {
                base.push_back(it->first);
            }
        }
        std::vector<std::string> ret2;
        for(auto& s : ret) {
            for(auto & b : base) {
                ret2.push_back(s + b);
            }
        }
        ret = std::move(ret2);
    }
    return ret;
}

int main(int argc, char** argv) {
    // iupac http://bioinformatics.org/sms/iupac.html
    nucleotide['A'] = {'A'};
    nucleotide['C'] = {'C'};
    nucleotide['G'] = {'G'};
    nucleotide['T'] = {'T'};
    nucleotide['R'] = {'A', 'G'};
    nucleotide['Y'] = {'C', 'T'};
    nucleotide['S'] = {'G', 'C'};
    nucleotide['W'] = {'A', 'T'};
    nucleotide['K'] = {'G', 'T'};
    nucleotide['M'] = {'A', 'C'};
    nucleotide['B'] = {'C', 'G', 'T'};
    nucleotide['D'] = {'A', 'G', 'T'};
    nucleotide['H'] = {'A', 'C', 'T'};
    nucleotide['V'] = {'A', 'C', 'G'};
    nucleotide['N'] = {'A', 'C', 'G', 'T'};

    whatisgood['A'] = {1, 0, 0, 0};
    whatisgood['C'] = {0, 1, 0, 0};
    whatisgood['G'] = {0, 0, 1, 0};
    whatisgood['T'] = {0, 0, 0, 1};
    whatisgood['R'] = {0.5, 0, 0.5, 0};
    whatisgood['Y'] = {0, 0.5, 0, 0.5};
    whatisgood['S'] = {0, 0.5, 0.5, 0};
    whatisgood['W'] = {0.5, 0, 0, 0.5};
    whatisgood['K'] = {0, 0, 0.5, 0.5};
    whatisgood['M'] = {0.5, 0.5, 0, 0};
    whatisgood['B'] = {0, 1/3.0, 1/3.0, 1/3.0};
    whatisgood['D'] = {1/3.0, 0, 1/3.0, 1/3.0};
    whatisgood['H'] = {1/3.0, 1/3.0, 0, 1/3.0};
    whatisgood['V'] = {1/3.0, 1/3.0, 1/3.0, 0};
    whatisgood['N'] = {0.25, 0.25, 0.25, 0.25};

    codon["TTT"] = 'F';
    codon["TTC"] = 'F';
    codon["TTA"] = 'L';
    codon["TTG"] = 'L';
    codon["CTT"] = 'L';
    codon["CTC"] = 'L';
    codon["CTA"] = 'L';
    codon["CTG"] = 'L';
    codon["ATT"] = 'I';
    codon["ATC"] = 'I';
    codon["ATA"] = 'I';
    codon["ATG"] = 'M';
    codon["GTT"] = 'V';
    codon["GTC"] = 'V';
    codon["GTA"] = 'V';
    codon["GTG"] = 'V';

    codon["TCT"] = 'S';
    codon["TCC"] = 'S';
    codon["TCA"] = 'S';
    codon["TCG"] = 'S';
    codon["CCT"] = 'P';
    codon["CCC"] = 'P';
    codon["CCA"] = 'P';
    codon["CCG"] = 'P';
    codon["ACT"] = 'T';
    codon["ACC"] = 'T';
    codon["ACA"] = 'T';
    codon["ACG"] = 'T';
    codon["GCT"] = 'A';
    codon["GCC"] = 'A';
    codon["GCA"] = 'A';
    codon["GCG"] = 'A';

    codon["TAT"] = 'Y';
    codon["TAC"] = 'Y';
    codon["TAA"] = '-';
    codon["TAG"] = '-';
    codon["CAT"] = 'H';
    codon["CAC"] = 'H';
    codon["CAA"] = 'Q';
    codon["CAG"] = 'Q';
    codon["AAT"] = 'N';
    codon["AAC"] = 'N';
    codon["AAA"] = 'K';
    codon["AAG"] = 'K';
    codon["GAT"] = 'D';
    codon["GAC"] = 'D';
    codon["GAA"] = 'E';
    codon["GAG"] = 'E';

    codon["TGT"] = 'C';
    codon["TGC"] = 'C';
    codon["TGA"] = '-';
    codon["TGG"] = 'W';
    codon["CGT"] = 'R';
    codon["CGC"] = 'R';
    codon["CGA"] = 'R';
    codon["CGG"] = 'R';
    codon["AGT"] = 'S';
    codon["AGC"] = 'S';
    codon["AGA"] = 'R';
    codon["AGG"] = 'R';
    codon["GGT"] = 'G';
    codon["GGC"] = 'G';
    codon["GGA"] = 'G';
    codon["GGG"] = 'G';

    // load train data
    std::vector<std::vector<double>> X;
    std::vector<double> t;

    std::ifstream ifs("experiment3-train.csv");
    std::string line;
    std::getline(ifs, line); // ignore the first line
    while(std::getline(ifs, line)) {
        std::stringstream ss(line);
        std::string elem;
        int row = 0;
        std::vector<double> row_x;
        double row_t;
        while(std::getline(ss, elem, ',')) {
            std::stringstream ss2(elem);
            row++;
            if(row==1) {
                std::string row_name;
                ss2 >> row_name;
            } else {
                double e;
                ss2 >> e;
                if(row <= 21) row_x.push_back(e);
                else row_t = e;
            }
        }
        X.push_back(row_x);
        t.push_back(row_t);
    }
    ifs.close();

    GaussianProcess<std::vector<double>> GP(X, t);

    std::vector<std::pair<double, std::string>> pred;
    ifs.open("experiment3-pred.csv");
    std::getline(ifs, line); // ignore the first line
    while(std::getline(ifs, line)) {
        std::stringstream ss(line);
        std::string elem;
        int row = 0;
        std::vector<double> row_x;
        std::string row_name;
        while(std::getline(ss, elem, ',')) {
            std::stringstream ss2(elem);
            row++;
            if(row==1) {
                ss2 >> row_name;
            } else {
                double e;
                ss2 >> e;
                if(row <= 21) row_x.push_back(e);
                else break;
            }
        }
        double t_pred = GP.predict(row_x).mu;
        pred.push_back(std::make_pair(t_pred, row_name));
    }

    // show result
    std::sort(pred.begin(), pred.end(), std::greater<std::pair<double, std::string>>());


    BKTree bkt;
    for(int i=0; i<500; ++i) {
        bkt.add(pred[i].second);
    }

    auto match = bkt.getSimilarWords("GAYF", 1);

    std::vector<std::pair<std::string, double>> candidate;
    for (auto x : match) {
         std::cout << x << std::endl;
        auto base_list = amino2base(x);
        double pt = 1.0 / base_list.size();
        for(auto s : base_list) {
            candidate.push_back({s, pt});
        }
    }
    std::vector<std::array<double, 4>> probability_list(12);
    for(auto cand : candidate) {
        for(int i=0; i<12; ++i) {
            if(cand.first[i]=='A') probability_list[i][0] += cand.second;
            if(cand.first[i]=='C') probability_list[i][1] += cand.second;
            if(cand.first[i]=='G') probability_list[i][2] += cand.second;
            if(cand.first[i]=='T') probability_list[i][3] += cand.second;
        }
    }

    for(auto probability : probability_list) {
        std::cout << "( ";
        for(auto& p : probability) {
            p /= match.size();
            std::cout << p << " ";
        }
        std::cout << ")" << std::endl;
        double min_dist = 1e10;
        char best;
        for(auto it=whatisgood.begin(); it!=whatisgood.end(); it++) {
            double dist = 0;
            for(int i=0; i<4; ++i) {
                dist += std::pow(probability[i]-it->second[i], 2.0);
            }
            if(dist < min_dist) {
                min_dist = dist;
                best = it->first;
            }
        }
        std::cout << best << std::endl;
    }
}
