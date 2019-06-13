// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppParallel)]]
#define BOOST_DISABLE_ASSERTS
#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <progress.hpp>
#include <progress_bar.hpp>
#include <fstream>
#include <string>
#include <unordered_map>
#include <iterator>
#include <algorithm>
#include <unordered_set>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/join.hpp>
#include <RcppParallel.h>
using namespace std;
using namespace Rcpp;
using namespace arma;
using namespace boost;
using namespace tbb;

typedef tbb::concurrent_unordered_map<int, vector<string> > cumivs;
// Tokenize the individual ECs
vector<string> tokenize(string s) {
  tokenizer<> tok(s);
  vector<string> out;
  for (tokenizer<>::iterator p = tok.begin(); p != tok.end(); ++p){
    out.push_back(*p);
  }
  return out;
}

// Get unique elements of vector
void gene_unique(vector<string>& g) {
  sort(g.begin(), g.end());
  g.erase(unique(g.begin(), g.end()), g.end());
}

// Get EC index and corresponding ECs
cumivs matrix_ec(int est_ngenes, std::string kallisto_out_path, int ncores = 0, 
                 bool verbose = true) {
  // Read in the matrix.ec file
  const string fn = kallisto_out_path + "/matrix.ec";
  ifstream infile(fn);
  if (infile.fail()) {
    stop("The file matrix.ec does not exist in kallisto_out_path.");
  }
  string ecind_str, ecs;
  int ct = 0;
  vector<string> ec2g_tmp;
  ec2g_tmp.reserve(15 * est_ngenes);
  if (verbose) Rcout << "Reading matrix.ec" << endl;
  while (infile >> ecind_str >> ecs) {
    if (ct % 1000 == 0) {
      checkUserInterrupt();
    }
    ec2g_tmp.push_back(ecs);
    ct++;
  }
  // Process ec2g_tmp in parallel
  if (verbose) Rcout << "Processing ECs" << endl;
  cumivs ec_vec;
  Progress p(ec2g_tmp.size(), verbose);
#ifdef _OPENMP
  if (ncores > 0)
    omp_set_num_threads(ncores);
#endif
#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < ec2g_tmp.size(); i++) {
    if (!Progress::check_abort()) {
      p.increment();
      ec_vec[i] = tokenize(ec2g_tmp[i]);
    }
  }
  return ec_vec;
}

// Get genes each EC maps to
cumivs EC2geneC(Rcpp::DataFrame tr2g, cumivs ec_vec, int ncores, bool verbose) {
  vector<string> g = tr2g["gene"];
  cumivs ec2g;
  if (verbose) Rcout << "Matching genes to ECs" << endl;
  Progress p(ec_vec.size(), verbose);
#ifdef _OPENMP
  if (ncores > 0)
    omp_set_num_threads(ncores);
#endif
#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < ec_vec.size(); i++) {
    if (!Progress::check_abort()) {
      p.increment();
      vector<string> genes;
      for (auto gi : ec_vec[i]) { 
        genes.push_back(g[stoi(gi)]);
      }
      gene_unique(genes);
      ec2g[i] = genes;
    }
  }
  return ec2g;
}

// This function is exported to R so users can query genes corresponding to
// ECs
//[[Rcpp::export]]
Rcpp::List EC2gene_export(Rcpp::DataFrame tr2g, std::string kallisto_out_path,
               int ncores = 0, bool verbose = true) {
  tbb::concurrent_unordered_map<int, vector<string> > ec_vec, ec2g;
  ec_vec = matrix_ec(tr2g.nrow(), kallisto_out_path, ncores, verbose);
  ec2g = EC2geneC(tr2g, ec_vec, ncores, verbose);
  return List::create(_["ec_vec"] = wrap(ec_vec),
                      _["ec2g"] = wrap(ec2g));
}

// Intersection of gene vectors
vector<string> gene_intersect(vector<string> a, vector<string> b) {
  vector<string> out;
  sort(a.begin(), a.end());
  sort(b.begin(), b.end());
  set_intersection(a.begin(), a.end(), b.begin(), b.end(), back_inserter(out));
  return out;
}

// Construct the row and column indices
void fill_spmat(unordered_map<string, unordered_map<string, double> >& cell_gene,
                vector<string>& barcodes, vector<string>& geneIDs,
                vector<double>& values, vector<size_t>& rowind,
                vector<size_t>& colind, bool verbose) {
  size_t i = 0, gene_row = 0; string gn;
  unordered_map<string, size_t> rowind_map;
  Progress p(cell_gene.size(), verbose);
  // Consider parallelizing
  for (auto el : cell_gene) {
    if (i % 1000 == 0) {
      checkUserInterrupt();
    }
    barcodes.push_back(el.first);
    // Construct the row index, iterate through each gene for this barcode
    for (auto k: el.second) {
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
}

// Convert to the actual sparse matrix
sp_mat conv_spmat(vector<double> values, vector<size_t> rowind,
                  vector<size_t> colind, size_t nrow, size_t ncol) {
  // Convert to sparse matrix
  // Convert to arma types
  vec val_use(values);
  uvec rind_use = conv_to<uvec>::from(rowind);
  uvec colind_use = conv_to<uvec>::from(colind);
  umat locations = join_rows(rind_use, colind_use);
  locations = locations.t();
  sp_mat res_mat(locations, val_use, nrow, ncol);
  return res_mat;
}

//[[Rcpp::export]]
List fill_cell_gene(std::string fn, std::string kallisto_out_path, 
                    DataFrame tr2g, int est_ncells, int est_ngenes,
                    std::vector<std::string> whitelist,
                    bool gene_count = true, bool tcc = true,
                    bool single_gene = true,
                    int ncores = 0, bool verbose = true,
                    int progress_unit = 5e6) {
  if (!(gene_count || tcc)) {
    stop("At least one of gene_count and tcc must be true.");
  }
  // Here assume that the path is normalized in R
  ifstream infile(kallisto_out_path + "/" + fn);
  if (infile.fail()) {
    stop("File in fn does not exist in kallisto_out_path.");
  }
  string bc, umi, ec_str, cts, ec_use;
  string pbar = "", pumi = "", pbf = "";
  int ec, n, i = 0; // i to keep track of # of iterations
  vector<string> gs, gl, trs;
  unordered_map<string, unordered_map<string, double>> cell_gene, cell_ec;
  if (gene_count) cell_gene.reserve(est_ncells);
  if (tcc) cell_ec.reserve(est_ncells);
  int wl_size = whitelist.size();
  // Convert whitelist into unordered_set to speed up lookup
  unordered_set<string> wl(whitelist.begin(), whitelist.end());
  // Get genes and ecs
  cumivs ec_vec, ec2g;
  ec_vec = matrix_ec(est_ngenes, kallisto_out_path, ncores, verbose);
  if (gene_count) {
    ec2g = EC2geneC(tr2g, ec_vec, ncores, verbose);
  }
  if (verbose) {
    Rcout << "Reading data" << endl;
  }
  while (infile >> bc >> umi >> ec_str >> cts) {
    if (i % 1000 == 0) {
      checkUserInterrupt();
    }
    ec = stoi(ec_str);
    if (bc == pbar) {
      // Same barcode
      if (umi == pumi) {
        // For TCC matrix, get intersection of transcript list
        if (tcc) trs = gene_intersect(trs, ec_vec[ec]);
        // Same umi, get intersection of gene list
        if (gene_count) gs = gene_intersect(gs, ec2g[ec]);
      } else {
        if (tcc) {
          if (trs.size() > 0) {
            ec_use = boost::algorithm::join(trs, ",");
            cell_ec[pbar][ec_use] += 1;
          }
          trs = ec_vec[ec];
        }
        if (gene_count) {
          // New UMI, process the previous gene list
          n = gs.size();
          // Single gene mode: skip if n > 1
          if (!(single_gene && n > 1)) {
            for (int j = 0; j < n; j++) {
              cell_gene[pbar][gs[j]] += 1.0/(double)n;
            }
          }
          gs = ec2g[ec];
        }
        pumi = umi;
      }
    } else {
      // If barcode is not in whitelist, skip to the next barcode
      // Won't skip if whitelist has length 1, a place holder for NA
      if ((pbf == bc) || (wl_size > 1 && wl.find(bc) == wl.end())) {
        pbf = bc; // Avoid checking whitelist again if the invalid barcode does
        continue; // reappear
      }
      if (tcc) {
        // For TCC matrix, previous transcript list
        if (trs.size() > 0) {
          ec_use = boost::algorithm::join(trs, ",");
          cell_ec[pbar][ec_use] += 1;
        }
        trs = ec_vec[ec];
      }
      if (gene_count) {
        // Previous gene list
        n = gs.size();
        if (!(single_gene && n > 1)) {
          for (int j = 0; j < n; j++) {
            cell_gene[pbar][gs[j]] += 1.0/(double)n;
          }
        }
        gs = ec2g[ec];
      }
      pumi = umi; pbar = bc;
    }
    // Some sense of progress
    if (verbose) {
      if (i % progress_unit == 0 && i > 0) {
        Rcout << "Read " << i/1e6 << " million reads" << endl;
      }
    }
    i++;
  }
  if (tcc) {
    // Remember the last EC
    if (trs.size() > 0) {
      ec_use = boost::algorithm::join(trs, ",");
      cell_ec[pbar][ec_use] += 1;
    }
  }
  if (gene_count) {
    // Remember the last gene
    n = gs.size();
    if (!(single_gene && n > 1)) {
      for (int j = 0; j < n; j++) {
        cell_gene[pbar][gs[j]] += 1.0/(double)n;
      }
    }
  }
  
  // Convert the unordered map into a sparse matrix
  // I'm using stl vectors here since they grow more nicely than arma::vec
  vector<string> barcodes_gc, barcodes_tcc, geneIDs, ec_inds;
  vector<double> values_gc, values_tcc;
  vector<size_t> rowind_gc, colind_gc, rowind_tcc, colind_tcc;
  if (gene_count) {
    barcodes_gc.reserve(est_ncells); geneIDs.reserve(est_ngenes);
    values_gc.reserve(i); rowind_gc.reserve(i); colind_gc.reserve(i);
    if (verbose) {
      Rcout << "Constructing gene count matrix" << endl;
    }
    fill_spmat(cell_gene, barcodes_gc, geneIDs, values_gc, rowind_gc, colind_gc,
               verbose);
  }
  if (tcc) {
    barcodes_tcc.reserve(est_ncells); ec_inds.reserve(est_ncells);
    values_tcc.reserve(i); rowind_tcc.reserve(i); colind_tcc.reserve(i);
    if (verbose) {
      Rcout << "Constructing TCC matrix" << endl;
    }
    fill_spmat(cell_ec, barcodes_tcc, ec_inds, values_tcc, rowind_tcc, colind_tcc,
               verbose);
  }
  
  sp_mat gc_mat, tcc_mat;
  gc_mat = conv_spmat(values_gc, rowind_gc, colind_gc, 
                      geneIDs.size(), barcodes_gc.size());
  tcc_mat = conv_spmat(values_tcc, rowind_tcc, colind_tcc,
                       ec_inds.size(), barcodes_tcc.size());
  
  // Output
  List out;
  if (gene_count && tcc) {
    out = List::create(_["gene_count"] = List::create(_["matrix"] = gc_mat,
                                            _["barcodes"] = barcodes_gc,
                                            _["genes"] = geneIDs),
                       _["TCC"] = List::create(_["matrix"] = tcc_mat,
                                            _["barcodes"] = barcodes_tcc,
                                            _["EC_index"] = ec_inds));
  } else if (gene_count) {
    out = List::create(_["matrix"] = gc_mat,
                       _["barcodes"] = barcodes_gc,
                       _["genes"] = geneIDs);
  } else {
    out = List::create(_["matrix"] = tcc_mat,
                       _["barcodes"] = barcodes_tcc,
                       _["EC_index"] = ec_inds);
  }
  return out;
}
