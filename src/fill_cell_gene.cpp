// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <fstream>
#include <string>
#include <unordered_map>
#include <iterator>
#include <algorithm>
#include <unordered_set>
using namespace std;
using namespace Rcpp;
using namespace arma;
// Intersection of gene vectors
vector<string> gene_intersect(vector<string> a, vector<string> b) {
  vector<string> out;
  sort(a.begin(), a.end());
  sort(b.begin(), b.end());
  set_intersection(a.begin(), a.end(), b.begin(), b.end(), back_inserter(out));
  return out;
}

//[[Rcpp::export]]
List fill_cell_gene(const char* fn, List genes, 
                    int est_ncells, int est_ngenes,
                    std::vector<std::string> whitelist,
                    bool display_progress = true,
                    int progress_unit = 5e6) {
  ifstream infile(fn);
  string bc, umi, ec_str, cts;
  string pbar = "", pumi = "";
  int ec, n, i = 0; // i to keep track of # of iterations
  vector<string> gs, gl;
  unordered_map<string, unordered_map<string, double>> cell_gene;
  cell_gene.reserve(est_ncells);
  int wl_size = whitelist.size();
  // Convert whitelist into unordered_set to speed up lookup
  unordered_set<string> wl(whitelist.begin(), whitelist.end());
  if (display_progress) {
    Rcout << "Reading data" << endl;
  }
  while (infile >> bc >> umi >> ec_str >> cts) {
    if (i % 1000 == 0) {
      checkUserInterrupt();
    }
    // If barcode is not in whitelist, skip to the next barcode
    // Won't skip if whitelist has length 1, a place holder for NA
    if (wl_size > 1 && wl.find(bc) == wl.end()) {
      continue;
    }
    ec = stoi(ec_str);
    if (bc == pbar) {
      // Same barcode
      if (umi == pumi) {
        // Same umi, get intersection of gene list
        gl = as<vector<string>>(genes[ec]);
        gs = gene_intersect(gs, gl);
      } else {
        // New UMI, process the previous gene list
        n = gs.size();
        for (int j = 0; j < n; j++) {
          cell_gene[bc][gs[j]] += 1.0/(double)n;
        }
        pumi = umi;
        gs = as<vector<string>>(genes[ec]);
      }
    } else {
      // Previous gene list
      n = gs.size();
      for (int j = 0; j < n; j++) {
        cell_gene[pbar][gs[j]] += 1.0/(double)n;
      }
      pumi = umi; pbar = bc;
      gs = as<vector<string>>(genes[ec]);
    }
    // Some sense of progress
    if (display_progress) {
      if (i % progress_unit == 0 && i > 0) {
        Rcout << "Read " << i/1e6 << " million lines" << endl;
      }
    }
    i++;
  }
  // Remember the last gene
  n = gs.size();
  for (int j = 0; j < n; j++) {
    cell_gene[pbar][gs[j]] += 1.0/(double)n;
  }
  
  // Convert the unordered map into a sparse matrix
  // I'm using stl vectors here since they grow more nicely than arma::vec
  vector<string> barcodes, geneIDs;
  barcodes.reserve(est_ncells); geneIDs.reserve(est_ngenes);
  vector<double> values;
  values.reserve(i); // Now i is the number of line in the output file
  vector<int> rowind, colind;
  rowind.reserve(i); colind.reserve(i);
  unordered_map<string, int> rowind_map;
  unordered_map<string, double> g; // for individual genes
  i = 0; // Here keep track of number of iteration in case there's interruption
  int gene_row = 0; // keep track of how many entries
  string gn; // gene name
  if (display_progress) {
    Rcout << "Constructing sparse matrix" << endl;
  }
  Progress p(cell_gene.size(), display_progress);
  // Consider parallelizing
  for (auto el : cell_gene) {
    if (i % 1000 == 0) {
      checkUserInterrupt();
    }
    barcodes.push_back(el.first);
    // Construct the row index, iterate through each gene for this barcode
    g = el.second;
    for (auto k: g) {
      gn = k.first;
      if (rowind_map.find(gn) == rowind_map.end()) {
        // If not found
        rowind_map[gn] = gene_row;
        geneIDs.push_back(gn);
        gene_row++;
      }
      rowind.push_back(rowind_map[gn]);
      colind.push_back(i);
      values.push_back(k.second); // The UMI count
    }
    p.increment();
    i++;
  }
  
  // Convert to sparse matrix
  // Convert to arma types
  vec val_use(values);
  uvec rind_use = conv_to<uvec>::from(rowind);
  uvec colind_use = conv_to<uvec>::from(colind);
  umat locations = join_rows(rind_use, colind_use);
  locations = locations.t();
  sp_mat res_mat(locations, val_use, geneIDs.size(), barcodes.size());
  
  return List::create(_["matrix"] = res_mat,
                      _["barcodes"] = barcodes,
                      _["genes"] = geneIDs);
}
