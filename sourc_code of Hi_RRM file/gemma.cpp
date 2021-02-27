/*
    Genome-wide Efficient Mixed Model Association (GEMMA)
    Copyright © 2011-2017, Xiang Zhou
    Copyright © 2017, Peter Carbonetto
    Copyright © 2017-2018, Pjotr Prins

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <cmath>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <sys/stat.h>
#ifdef OPENBLAS
extern "C" {
  // these functions are defined in cblas.h - but if we include that we
  // conflicts with other BLAS includes (GSL)
  int openblas_get_num_threads(void);
  int openblas_get_parallel(void);
  char* openblas_get_config(void);
  char* openblas_get_corename(void);
}
#else
#pragma message "Not compiling with OPENBLAS"
#endif

#include "gsl/gsl_blas.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_version.h"

#include "bslmm.h"
#include "bslmmdap.h"
#include <csignal> // for gsl_error_handler
#include "gemma.h"
#include "gemma_io.h"
#include "lapack.h"
#include "ldr.h"
#include "lm.h"
#include "lmm.h"
#include "mathfunc.h"
#include "mvlmm.h"
#include "prdt.h"
#include "varcov.h"
#include "vc.h"
#include "debug.h"
#include "version.h"

using namespace std;

GEMMA::GEMMA(void) : version(GEMMA_VERSION), date(GEMMA_DATE), year(GEMMA_YEAR) {}

void gemma_gsl_error_handler (const char * reason,
                              const char * file,
                              int line, int gsl_errno) {
  cerr << "GSL ERROR: " << reason << " in " << file
       << " at line " << line << " errno " << gsl_errno <<endl;
  std::raise(SIGINT); // keep the stack trace for gdb
}

#if defined(OPENBLAS) && !defined(OPENBLAS_LEGACY)
#include <openblas_config.h>
#endif

void GEMMA::PrintHeader(void) {

  cout <<
    "This GEMMAX software is developed based on GEMMA" <<  endl;
  return;
}

void GEMMA::PrintLicense(void) {
  cout <<
    "This GEMMAX software is developed based on GEMMA" <<  endl;
  cout << endl;

  return;
}

void GEMMA::PrintHelp(size_t option) {

  cout << "This GEMMAX software is developed based on GEMMA" << endl;
  return;
}

// OPTIONS
// -------
// gk:      21-22
// gs:      25-26
// gq:      27-28
// eigen:   31-32
// lmm:     1-5
// bslmm:   11-15
// predict: 41-43
// lm:      51
// vc:      61
// ci:      66-67
// calccor: 71
// gw:      72

void GEMMA::Assign(int argc, char **argv, PARAM &cPar) {
  string str;

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-bfile") == 0 || strcmp(argv[i], "--bfile") == 0 ||
        strcmp(argv[i], "-b") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.file_bfile = str;
    } else if (strcmp(argv[i], "-mbfile") == 0 ||
               strcmp(argv[i], "--mbfile") == 0 ||
               strcmp(argv[i], "-mb") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.file_mbfile = str;
    } else if (strcmp(argv[i], "-silence") == 0 || strcmp(argv[i], "--quiet") == 0) {
      debug_set_quiet_mode(true);
    } else if (strcmp(argv[i], "-g") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.file_geno = str;
    } else if (strcmp(argv[i], "-mg") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.file_mgeno = str;
    } else if (strcmp(argv[i], "-p") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.file_pheno = str;
    } else if (strcmp(argv[i], "-k") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.file_kin = str;
    } else if (strcmp(argv[i], "-n") == 0) { // set pheno column (list/range)
      (cPar.p_column).clear();
      while (argv[i + 1] != NULL && argv[i + 1][0] != '-') {
        ++i;
        str.clear();
        str.assign(argv[i]);
        (cPar.p_column).push_back(atoi(str.c_str()));
      }
    }  else if (strcmp(argv[i], "-outdir") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.path_out = str;
    } else if (strcmp(argv[i], "-o") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.file_out = str;
    } else if (strcmp(argv[i], "-miss") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.miss_level = atof(str.c_str());
    } else if (strcmp(argv[i], "-maf") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      if (cPar.maf_level != -1) {
        cPar.maf_level = atof(str.c_str());
      }
    } else if (strcmp(argv[i], "-hwe") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.hwe_level = atof(str.c_str());
    } else if (strcmp(argv[i], "-r2") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.r2_level = atof(str.c_str());
    } else if (strcmp(argv[i], "-gk") == 0) {
      if (cPar.a_mode != 0) {
        cPar.error = true;
        cout << "error! only one of -gk -gs -eigen -vc -lm -lmm -bslmm "
                "-predict -calccor options is allowed."
             << endl;
        break;
      }
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        cPar.a_mode = M_KIN;
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.a_mode = 20 + atoi(str.c_str());
    } else if (strcmp(argv[i], "-gs") == 0) {
      if (cPar.a_mode != 0) {
        cPar.error = true;
        cout << "error! only one of -gk -gs -eigen -vc -lm -lmm -bslmm "
                "-predict -calccor options is allowed."
             << endl;
        break;
      }
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        cPar.a_mode = 25;
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.a_mode = 24 + atoi(str.c_str());
    } else if (strcmp(argv[i], "-gq") == 0) {
      if (cPar.a_mode != 0) {
        cPar.error = true;
        cout << "error! only one of -gk -gs -eigen -vc -lm -lmm -bslmm "
                "-predict -calccor options is allowed."
             << endl;
        break;
      }
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        cPar.a_mode = 27;
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.a_mode = 26 + atoi(str.c_str());
    } else if (strcmp(argv[i], "-gw") == 0) {
      if (cPar.a_mode != 0) {
        cPar.error = true;
        cout << "error! only one of -gk -gs -eigen -vc -lm -lmm -bslmm "
                "-predict -calccor options is allowed."
             << endl;
        break;
      }
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        cPar.a_mode = 72;
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.a_mode = 71 + atoi(str.c_str());
    } else if (strcmp(argv[i], "-sample") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.ni_subsample = atoi(str.c_str());
    } else if (strcmp(argv[i], "-eigen") == 0) {
      if (cPar.a_mode != 0) {
        cPar.error = true;
        cout << "error! only one of -gk -gs -eigen -vc -lm -lmm -bslmm "
                "-predict -calccor options is allowed."
             << endl;
        break;
      }
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        cPar.a_mode = 31;
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.a_mode = 30 + atoi(str.c_str());
    } else if (strcmp(argv[i], "-calccor") == 0) {
      if (cPar.a_mode != 0) {
        cPar.error = true;
        cout << "error! only one of -gk -gs -eigen -vc -lm -lmm -bslmm "
                "-predict -calccor options is allowed."
             << endl;
        break;
      }
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        cPar.a_mode = 71;
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.a_mode = 70 + atoi(str.c_str());
    } else if (strcmp(argv[i], "-fa") == 0 || strcmp(argv[i], "-lmm") == 0|| strcmp(argv[i], "-JiontGEMMAX") == 0|| strcmp(argv[i], "-SeparateGEMMAX") == 0) {
      if (cPar.a_mode != 0) {
        cPar.error = true;
        cout << "error! only one of -gk -gs -eigen -vc -lm -lmm -bslmm "
                "-predict -calccor options is allowed."
             << endl;
        break;
      }
      if (strcmp(argv[i], "-JiontGEMMAX") == 0){
        cPar.a_mode = 4;
        continue;
      }
      if (strcmp(argv[i], "-SeparateGEMMAX") == 0){
        cPar.a_mode = 1;
        continue;
      }
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        cPar.a_mode = 1;
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.a_mode = atoi(str.c_str());
    } else if (strcmp(argv[i], "-lmin") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.l_min = atof(str.c_str());
    } else if (strcmp(argv[i], "-lmax") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.l_max = atof(str.c_str());
    } else if (strcmp(argv[i], "-region") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.n_region = atoi(str.c_str());
    } else if (strcmp(argv[i], "-pnr") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.p_nr = atof(str.c_str());
    } else if (strcmp(argv[i], "-emi") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.em_iter = atoi(str.c_str());
    } else if (strcmp(argv[i], "-nri") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.nr_iter = atoi(str.c_str());
    } else if (strcmp(argv[i], "-nind") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.ni_max = atoi(str.c_str()); // for testing purposes
      enforce(cPar.ni_max > 0);
    } else if (strcmp(argv[i], "-issue") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      auto issue = atoi(str.c_str()); // for testing purposes
      enforce(issue > 0);
      debug_set_issue(issue);
    } else if (strcmp(argv[i], "-emp") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.em_prec = atof(str.c_str());
    } else if (strcmp(argv[i], "-nrp") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.nr_prec = atof(str.c_str());
    } else if (strcmp(argv[i], "-crt") == 0) {
      cPar.crt = 1;
    } else if (strcmp(argv[i], "-hmin") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.h_min = atof(str.c_str());
    } else if (strcmp(argv[i], "-hmax") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.h_max = atof(str.c_str());
    } else if (strcmp(argv[i], "-rmin") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.rho_min = atof(str.c_str());
    } else if (strcmp(argv[i], "-rmax") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.rho_max = atof(str.c_str());
    } else if (strcmp(argv[i], "-pmin") == 0) {
      if (argv[i + 1] == NULL) {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.logp_min = atof(str.c_str()) * log(10.0);
    } else if (strcmp(argv[i], "-pmax") == 0) {
      if (argv[i + 1] == NULL) {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.logp_max = atof(str.c_str()) * log(10.0);
    } else if (strcmp(argv[i], "-smin") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.s_min = atoi(str.c_str());
    } else if (strcmp(argv[i], "-smax") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.s_max = atoi(str.c_str());
    } else if (strcmp(argv[i], "-gmean") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.geo_mean = atof(str.c_str());
    } else if (strcmp(argv[i], "-hscale") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.h_scale = atof(str.c_str());
    } else if (strcmp(argv[i], "-rscale") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.rho_scale = atof(str.c_str());
    } else if (strcmp(argv[i], "-pscale") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.logp_scale = atof(str.c_str()) * log(10.0);
    } else if (strcmp(argv[i], "-w") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.w_step = atoi(str.c_str());
    } else if (strcmp(argv[i], "-s") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.s_step = atoi(str.c_str());
    } else if (strcmp(argv[i], "-rpace") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.r_pace = atoi(str.c_str());
    } else if (strcmp(argv[i], "-wpace") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.w_pace = atoi(str.c_str());
    } else if (strcmp(argv[i], "-seed") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.randseed = atol(str.c_str());
    } else if (strcmp(argv[i], "-mh") == 0) {
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.n_mh = atoi(str.c_str());
    } else if (strcmp(argv[i], "-predict") == 0) {
      if (cPar.a_mode != 0) {
        cPar.error = true;
        cout << "error! only one of -gk -gs -eigen -vc -lm -lmm -bslmm "
                "-predict -calccor options is allowed."
             << endl;
        break;
      }
      if (argv[i + 1] == NULL || argv[i + 1][0] == '-') {
        cPar.a_mode = 41;
        continue;
      }
      ++i;
      str.clear();
      str.assign(argv[i]);
      cPar.a_mode = 40 + atoi(str.c_str());
    }  else {
      cout << "error! unrecognized option: " << argv[i] << endl;
      cPar.error = true;
      continue;
    }
  }

  // Change prediction mode to 43 if the epm file is not provided.
  if (cPar.a_mode == 41 && cPar.file_epm.empty()) {
    cPar.a_mode = 43;
  }

  return;
}

void GEMMA::BatchRun(PARAM &cPar) {
  clock_t time_begin, time_start;
  time_begin = clock();

  if (is_check_mode()) enable_segfpe(); // fast NaN checking by default

  // Read Files.
  cout << "Reading Files ... " << endl;
  cPar.ReadFiles();
  if (cPar.error == true) {
    cout << "error! fail to read files. " << endl;
    return;
  }
  cPar.CheckData();
  if (cPar.error == true) {
    cout << "error! fail to check data. " << endl;
    return;
  }

  // Prediction for bslmm
  if (cPar.a_mode == 41 || cPar.a_mode == 42) {
    gsl_vector *y_prdt;

    y_prdt = gsl_vector_safe_alloc(cPar.ni_total - cPar.ni_test);

    // set to zero
    gsl_vector_set_zero(y_prdt);

    PRDT cPRDT;
    cPRDT.CopyFromParam(cPar);

    // add breeding value if needed
    if (!cPar.file_kin.empty() && !cPar.file_ebv.empty()) {
      cout << "Adding Breeding Values ... " << endl;

      gsl_matrix *G = gsl_matrix_safe_alloc(cPar.ni_total, cPar.ni_total);
      gsl_vector *u_hat = gsl_vector_safe_alloc(cPar.ni_test);

      // read kinship matrix and set u_hat
      vector<int> indicator_all;
      size_t c_bv = 0;
      for (size_t i = 0; i < cPar.indicator_idv.size(); i++) {
        indicator_all.push_back(1);
        if (cPar.indicator_bv[i] == 1) {
          gsl_vector_set(u_hat, c_bv, cPar.vec_bv[i]);
          c_bv++;
        }
      }

      ReadFile_kin(cPar.file_kin, indicator_all, cPar.mapID2num, cPar.k_mode,
                   cPar.error, G);
      if (cPar.error == true) {
        cout << "error! fail to read kinship/relatedness file. " << endl;
        return;
      }

      // read u
      cPRDT.AddBV(G, u_hat, y_prdt);

      gsl_matrix_safe_free(G);
      gsl_vector_safe_free(u_hat);
    }

    // add beta
    if (!cPar.file_bfile.empty()) {
      cPRDT.AnalyzePlink(y_prdt);
    } else {
      cPRDT.AnalyzeBimbam(y_prdt);
    }

    // add mu
    gsl_vector_add_constant(y_prdt, cPar.pheno_mean);

    // convert y to probability if needed
    if (cPar.a_mode == 42) {
      double d;
      for (size_t i = 0; i < y_prdt->size; i++) {
        d = gsl_vector_get(y_prdt, i);
        d = gsl_cdf_gaussian_P(d, 1.0);
        gsl_vector_set(y_prdt, i, d);
      }
    }

    cPRDT.CopyToParam(cPar);

    cPRDT.WriteFiles(y_prdt);

    gsl_vector_safe_free(y_prdt);
  }

  // Prediction with kinship matrix only; for one or more phenotypes
  if (cPar.a_mode == 43) {
    // first, use individuals with full phenotypes to obtain estimates of Vg and
    // Ve
    gsl_matrix *Y = gsl_matrix_safe_alloc(cPar.ni_test, cPar.n_ph);
    gsl_matrix *W = gsl_matrix_safe_alloc(Y->size1, cPar.n_cvt);
    gsl_matrix *G = gsl_matrix_safe_alloc(Y->size1, Y->size1);
    gsl_matrix *U = gsl_matrix_safe_alloc(Y->size1, Y->size1);
    gsl_matrix *UtW = gsl_matrix_safe_alloc(Y->size1, W->size2);
    gsl_matrix *UtY = gsl_matrix_safe_alloc(Y->size1, Y->size2);
    gsl_vector *eval = gsl_vector_safe_alloc(Y->size1);

    gsl_matrix *Y_full = gsl_matrix_safe_alloc(cPar.ni_cvt, cPar.n_ph);
    gsl_matrix *W_full = gsl_matrix_safe_alloc(Y_full->size1, cPar.n_cvt);

    // set covariates matrix W and phenotype matrix Y
    // an intercept should be included in W,
    cPar.CopyCvtPhen(W, Y, 0);
    cPar.CopyCvtPhen(W_full, Y_full, 1);

    gsl_matrix *Y_hat = gsl_matrix_safe_alloc(Y_full->size1, cPar.n_ph);
    gsl_matrix *G_full = gsl_matrix_safe_alloc(Y_full->size1, Y_full->size1);
    gsl_matrix *H_full = gsl_matrix_safe_alloc(Y_full->size1 * Y_hat->size2,
                                          Y_full->size1 * Y_hat->size2);

    // read relatedness matrix G, and matrix G_full
    ReadFile_kin(cPar.file_kin, cPar.indicator_idv, cPar.mapID2num, cPar.k_mode,
                 cPar.error, G);
    if (cPar.error == true) {
      cout << "error! fail to read kinship/relatedness file. " << endl;
      return;
    }
    // This is not so elegant. Reads twice to select on idv and then cvt
    ReadFile_kin(cPar.file_kin, cPar.indicator_cvt, cPar.mapID2num, cPar.k_mode,
                 cPar.error, G_full);
    if (cPar.error == true) {
      cout << "error! fail to read kinship/relatedness file. " << endl;
      return;
    }

    // center matrix G
    CenterMatrix(G);
    CenterMatrix(G_full);
    validate_K(G);

    // eigen-decomposition and calculate trace_G
    cout << "Start Eigen-Decomposition..." << endl;
    time_start = clock();
    cPar.trace_G = EigenDecomp_Zeroed(G, U, eval, 0);
    // write(eval,"eval zeroed");
    cPar.time_eigen = (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

    // calculate UtW and Uty
    CalcUtX(U, W, UtW);
    CalcUtX(U, Y, UtY);

    // calculate variance component and beta estimates
    // and then obtain predicted values
    if (cPar.n_ph == 1) {
      gsl_vector *beta = gsl_vector_safe_alloc(W->size2);
      gsl_vector *se_beta = gsl_vector_safe_alloc(W->size2);

      double lambda, logl, vg, ve;
      gsl_vector_view UtY_col = gsl_matrix_column(UtY, 0);

      // obtain estimates
      CalcLambda('R', eval, UtW, &UtY_col.vector, cPar.l_min, cPar.l_max,
                 cPar.n_region, lambda, logl);
      CalcLmmVgVeBeta(eval, UtW, &UtY_col.vector, lambda, vg, ve, beta,
                      se_beta);

      cout << "REMLE estimate for vg in the null model = " << vg << endl;
      cout << "REMLE estimate for ve in the null model = " << ve << endl;
      cPar.vg_remle_null = vg;
      cPar.ve_remle_null = ve;

      // obtain Y_hat from fixed effects
      gsl_vector_view Yhat_col = gsl_matrix_column(Y_hat, 0);
      gsl_blas_dgemv(CblasNoTrans, 1.0, W_full, beta, 0.0, &Yhat_col.vector);

      // obtain H
      gsl_matrix_set_identity(H_full);
      gsl_matrix_scale(H_full, ve);
      gsl_matrix_scale(G_full, vg);
      gsl_matrix_add(H_full, G_full);

      // free matrices
      gsl_vector_safe_free(beta);
      gsl_vector_safe_free(se_beta);
    } else {
      gsl_matrix *Vg = gsl_matrix_safe_alloc(cPar.n_ph, cPar.n_ph);
      gsl_matrix *Ve = gsl_matrix_safe_alloc(cPar.n_ph, cPar.n_ph);
      gsl_matrix *B = gsl_matrix_safe_alloc(cPar.n_ph, W->size2);
      gsl_matrix *se_B = gsl_matrix_safe_alloc(cPar.n_ph, W->size2);

      // obtain estimates
      CalcMvLmmVgVeBeta(eval, UtW, UtY, cPar.em_iter, cPar.nr_iter,
                        cPar.em_prec, cPar.nr_prec, cPar.l_min, cPar.l_max,
                        cPar.n_region, Vg, Ve, B, se_B);

      cout << "REMLE estimate for Vg in the null model: " << endl;
      for (size_t i = 0; i < Vg->size1; i++) {
        for (size_t j = 0; j <= i; j++) {
          cout << tab(j) << gsl_matrix_get(Vg, i, j);
        }
        cout << endl;
      }
      cout << "REMLE estimate for Ve in the null model: " << endl;
      for (size_t i = 0; i < Ve->size1; i++) {
        for (size_t j = 0; j <= i; j++) {
          cout << tab(j) << gsl_matrix_get(Ve, i, j);
        }
        cout << endl;
      }
      cPar.Vg_remle_null.clear();
      cPar.Ve_remle_null.clear();
      for (size_t i = 0; i < Vg->size1; i++) {
        for (size_t j = i; j < Vg->size2; j++) {
          cPar.Vg_remle_null.push_back(gsl_matrix_get(Vg, i, j));
          cPar.Ve_remle_null.push_back(gsl_matrix_get(Ve, i, j));
        }
      }

      // obtain Y_hat from fixed effects
      gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, W_full, B, 0.0, Y_hat);

      // obtain H
      KroneckerSym(G_full, Vg, H_full);
      for (size_t i = 0; i < G_full->size1; i++) {
        gsl_matrix_view H_sub = gsl_matrix_submatrix(
            H_full, i * Ve->size1, i * Ve->size2, Ve->size1, Ve->size2);
        gsl_matrix_add(&H_sub.matrix, Ve);
      }

      // free matrices
      gsl_matrix_safe_free(Vg);
      gsl_matrix_safe_free(Ve);
      gsl_matrix_safe_free(B);
      gsl_matrix_safe_free(se_B);
    }

    PRDT cPRDT;

    cPRDT.CopyFromParam(cPar);

    cout << "Predicting Missing Phentypes ... " << endl;
    time_start = clock();
    cPRDT.MvnormPrdt(Y_hat, H_full, Y_full);
    cPar.time_opt = (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

    cPRDT.WriteFiles(Y_full);

    gsl_matrix_safe_free(Y);
    gsl_matrix_safe_free(W);
    gsl_matrix_safe_free(G);
    gsl_matrix_safe_free(U);
    gsl_matrix_safe_free(UtW);
    gsl_matrix_safe_free(UtY);
    gsl_vector_safe_free(eval);

    gsl_matrix_safe_free(Y_full);
    gsl_matrix_safe_free(Y_hat);
    gsl_matrix_safe_free(W_full);
    gsl_matrix_safe_free(G_full);
    gsl_matrix_safe_free(H_full);
  }

  // Generate Kinship matrix (optionally using LOCO)
  if (cPar.a_mode == M_KIN || cPar.a_mode == M_KIN2) {
    cout << "Calculating Relatedness Matrix ... " << endl;

    gsl_matrix *G = gsl_matrix_safe_alloc(cPar.ni_total, cPar.ni_total);
    enforce_msg(G, "allocate G"); // just to be sure

    time_start = clock();

    cPar.CalcKin(G);

    cPar.time_G = (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);
    if (cPar.error == true) {
      cout << "error! fail to calculate relatedness matrix. " << endl;
      return;
    }

    // Now we have the Kinship matrix test it
    validate_K(G);

    if (cPar.a_mode == M_KIN) {
      cPar.WriteMatrix(G, "cXX");
    } else {
      cPar.WriteMatrix(G, "sXX");
    }

    gsl_matrix_safe_free(G);
  }

  // Compute the LDSC weights (not implemented yet)
  if (cPar.a_mode == 72) {
    cout << "Calculating Weights ... " << endl;

    VARCOV cVarcov;
    cVarcov.CopyFromParam(cPar);

    if (is_check_mode()) disable_segfpe(); // disable fast NaN checking for now

    if (!cPar.file_bfile.empty()) {
      cVarcov.AnalyzePlink();
    } else {
      cVarcov.AnalyzeBimbam();
    }

    cVarcov.CopyToParam(cPar);
  }

  // Compute the S matrix (and its variance), that is used for
  // variance component estimation using summary statistics.
  if (cPar.a_mode == 25 || cPar.a_mode == 26) {
    cout << "Calculating the S Matrix ... " << endl;

    gsl_matrix *S = gsl_matrix_safe_alloc(cPar.n_vc * 2, cPar.n_vc);
    gsl_vector *ns = gsl_vector_safe_alloc(cPar.n_vc + 1);
    gsl_matrix_set_zero(S);
    gsl_vector_set_zero(ns);

    gsl_matrix_view S_mat = gsl_matrix_submatrix(S, 0, 0, cPar.n_vc, cPar.n_vc);
    gsl_matrix_view Svar_mat =
        gsl_matrix_submatrix(S, cPar.n_vc, 0, cPar.n_vc, cPar.n_vc);
    gsl_vector_view ns_vec = gsl_vector_subvector(ns, 0, cPar.n_vc);

    gsl_matrix *K = gsl_matrix_safe_alloc(cPar.ni_test, cPar.n_vc * cPar.ni_test);
    gsl_matrix *A = gsl_matrix_safe_alloc(cPar.ni_test, cPar.n_vc * cPar.ni_test);
    gsl_matrix_set_zero(K);
    gsl_matrix_set_zero(A);

    gsl_vector *y = gsl_vector_safe_alloc(cPar.ni_test);
    gsl_matrix *W = gsl_matrix_safe_alloc(cPar.ni_test, cPar.n_cvt);

    cPar.CopyCvtPhen(W, y, 0);

    set<string> setSnps_beta;
    map<string, double> mapRS2wA, mapRS2wK;

    cPar.ObtainWeight(setSnps_beta, mapRS2wK);

    time_start = clock();
    cPar.CalcS(mapRS2wA, mapRS2wK, W, A, K, &S_mat.matrix, &Svar_mat.matrix,
               &ns_vec.vector);
    cPar.time_G = (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);
    if (cPar.error == true) {
      cout << "error! fail to calculate the S matrix. " << endl;
      return;
    }

    gsl_vector_set(ns, cPar.n_vc, cPar.ni_test);

    cPar.WriteMatrix(S, "S");
    cPar.WriteVector(ns, "size");
    cPar.WriteVar("snps");

    gsl_matrix_safe_free(S);
    gsl_vector_safe_free(ns);

    gsl_matrix_safe_free(A);
    gsl_matrix_safe_free(K);

    gsl_vector_safe_free(y);
    gsl_matrix_safe_free(K);
  }

  // Compute the q vector, that is used for variance component estimation using
  // summary statistics
  if (cPar.a_mode == 27 || cPar.a_mode == 28) {
    gsl_matrix *Vq = gsl_matrix_safe_alloc(cPar.n_vc, cPar.n_vc);
    gsl_vector *q = gsl_vector_safe_alloc(cPar.n_vc);
    gsl_vector *s = gsl_vector_safe_alloc(cPar.n_vc + 1);
    gsl_vector_set_zero(q);
    gsl_vector_set_zero(s);

    gsl_vector_view s_vec = gsl_vector_subvector(s, 0, cPar.n_vc);

    vector<size_t> vec_cat, vec_ni;
    vector<double> vec_weight, vec_z2;
    map<string, double> mapRS2weight;
    mapRS2weight.clear();

    time_start = clock();
    ReadFile_beta(cPar.file_beta, cPar.mapRS2cat, mapRS2weight, vec_cat, vec_ni,
                  vec_weight, vec_z2, cPar.ni_total, cPar.ns_total,
                  cPar.ns_test);
    cout << "## number of total individuals = " << cPar.ni_total << endl;
    cout << "## number of total SNPs/var = " << cPar.ns_total << endl;
    cout << "## number of analyzed SNPs/var = " << cPar.ns_test << endl;
    cout << "## number of variance components = " << cPar.n_vc << endl;
    cout << "Calculating the q vector ... " << endl;
    Calcq(cPar.n_block, vec_cat, vec_ni, vec_weight, vec_z2, Vq, q,
          &s_vec.vector);
    cPar.time_G = (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

    if (cPar.error == true) {
      cout << "error! fail to calculate the q vector. " << endl;
      return;
    }

    gsl_vector_set(s, cPar.n_vc, cPar.ni_total);

    cPar.WriteMatrix(Vq, "Vq");
    cPar.WriteVector(q, "q");
    cPar.WriteVector(s, "size");
    gsl_matrix_safe_free(Vq);
    gsl_vector_safe_free(q);
    gsl_vector_safe_free(s);
  }

  // Calculate SNP covariance.
  if (cPar.a_mode == 71) {
    VARCOV cVarcov;
    cVarcov.CopyFromParam(cPar);

    if (is_check_mode()) disable_segfpe(); // fast NaN checking for now

    if (!cPar.file_bfile.empty()) {
      cVarcov.AnalyzePlink();
    } else {
      cVarcov.AnalyzeBimbam();
    }

    cVarcov.CopyToParam(cPar);
  }

  // LM.
  if (cPar.a_mode == 51 || cPar.a_mode == 52 || cPar.a_mode == 53 ||
      cPar.a_mode == 54) { // Fit LM
    gsl_matrix *Y = gsl_matrix_safe_alloc(cPar.ni_test, cPar.n_ph);
    gsl_matrix *W = gsl_matrix_safe_alloc(Y->size1, cPar.n_cvt);

    // set covariates matrix W and phenotype matrix Y
    // an intercept should be included in W,
    cPar.CopyCvtPhen(W, Y, 0);

    // Fit LM or mvLM
    if (cPar.n_ph == 1) {
      LM cLm;
      cLm.CopyFromParam(cPar);

      gsl_vector_view Y_col = gsl_matrix_column(Y, 0);

      // if (is_check_mode()) disable_segfpe(); // disable fast NaN checking for now

      if (!cPar.file_gene.empty()) {
        cLm.AnalyzeGene(W,
                        &Y_col.vector); // y is the predictor, not the phenotype
      } else if (!cPar.file_bfile.empty()) {
        cLm.AnalyzePlink(W, &Y_col.vector);
      } else {
        cLm.AnalyzeBimbam(W, &Y_col.vector);
      }

      cLm.WriteFiles();
      cLm.CopyToParam(cPar);
    }
    // release all matrices and vectors
    gsl_matrix_safe_free(Y);
    gsl_matrix_safe_free(W);
  }

  // VC estimation with one or multiple kinship matrices
  // REML approach only
  // if file_kin or file_ku/kd is provided, then a_mode is changed to 5 already,
  // in param.cpp
  // for one phenotype only;
  if (cPar.a_mode == 61 || cPar.a_mode == 62 || cPar.a_mode == 63) {
    if (!cPar.file_beta.empty()) {
      // need to obtain a common set of SNPs between beta file and the genotype
      // file; these are saved in mapRS2wA and mapRS2wK
      // normalize the weight in mapRS2wK to have an average of one; each
      // element of mapRS2wA is 1
      // update indicator_snps, so that the numbers are in accordance with
      // mapRS2wK
      set<string> setSnps_beta;
      ReadFile_snps_header(cPar.file_beta, setSnps_beta);

      map<string, double> mapRS2wA, mapRS2wK;
      cPar.ObtainWeight(setSnps_beta, mapRS2wK);

      cPar.UpdateSNP(mapRS2wK);

      // Setup matrices and vectors.
      gsl_matrix *S = gsl_matrix_safe_alloc(cPar.n_vc * 2, cPar.n_vc);
      gsl_matrix *Vq = gsl_matrix_safe_alloc(cPar.n_vc, cPar.n_vc);
      gsl_vector *q = gsl_vector_safe_alloc(cPar.n_vc);
      gsl_vector *s = gsl_vector_safe_alloc(cPar.n_vc + 1);

      gsl_matrix *K = gsl_matrix_safe_alloc(cPar.ni_test, cPar.n_vc * cPar.ni_test);
      gsl_matrix *A = gsl_matrix_safe_alloc(cPar.ni_test, cPar.n_vc * cPar.ni_test);

      gsl_vector *y = gsl_vector_safe_alloc(cPar.ni_test);
      gsl_matrix *W = gsl_matrix_safe_alloc(cPar.ni_test, cPar.n_cvt);

      gsl_matrix_set_zero(K);
      gsl_matrix_set_zero(A);

      gsl_matrix_set_zero(S);
      gsl_matrix_set_zero(Vq);
      gsl_vector_set_zero(q);
      gsl_vector_set_zero(s);

      cPar.CopyCvtPhen(W, y, 0);

      gsl_matrix_view S_mat =
          gsl_matrix_submatrix(S, 0, 0, cPar.n_vc, cPar.n_vc);
      gsl_matrix_view Svar_mat =
          gsl_matrix_submatrix(S, cPar.n_vc, 0, cPar.n_vc, cPar.n_vc);
      gsl_vector_view s_vec = gsl_vector_subvector(s, 0, cPar.n_vc);

      vector<size_t> vec_cat, vec_ni;
      vector<double> vec_weight, vec_z2;

      // read beta, based on the mapRS2wK
      ReadFile_beta(cPar.file_beta, cPar.mapRS2cat, mapRS2wK, vec_cat, vec_ni,
                    vec_weight, vec_z2, cPar.ni_study, cPar.ns_study,
                    cPar.ns_test);

      cout << "Study Panel: " << endl;
      cout << "## number of total individuals = " << cPar.ni_study << endl;
      cout << "## number of total SNPs/var = " << cPar.ns_study << endl;
      cout << "## number of analyzed SNPs/var = " << cPar.ns_test << endl;
      cout << "## number of variance components = " << cPar.n_vc << endl;

      // compute q
      Calcq(cPar.n_block, vec_cat, vec_ni, vec_weight, vec_z2, Vq, q,
            &s_vec.vector);

      // compute S
      time_start = clock();
      cPar.CalcS(mapRS2wA, mapRS2wK, W, A, K, &S_mat.matrix, &Svar_mat.matrix,
                 &s_vec.vector);
      cPar.time_G += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);
      if (cPar.error == true) {
        cout << "error! fail to calculate the S matrix. " << endl;
        return;
      }

      // compute vc estimates
      CalcVCss(Vq, &S_mat.matrix, &Svar_mat.matrix, q, &s_vec.vector,
               cPar.ni_study, cPar.v_pve, cPar.v_se_pve, cPar.pve_total,
               cPar.se_pve_total, cPar.v_sigma2, cPar.v_se_sigma2,
               cPar.v_enrich, cPar.v_se_enrich);

      assert(!has_nan(cPar.v_se_pve));

      // if LDSC weights, then compute the weights and run the above steps again
      if (cPar.a_mode == 62) {
        // compute the weights and normalize the weights for A
        cPar.UpdateWeight(1, mapRS2wK, cPar.ni_study, &s_vec.vector, mapRS2wA);

        // read beta file again, and update weigths vector
        ReadFile_beta(cPar.file_beta, cPar.mapRS2cat, mapRS2wA, vec_cat, vec_ni,
                      vec_weight, vec_z2, cPar.ni_study, cPar.ns_total,
                      cPar.ns_test);

        // compute q
        Calcq(cPar.n_block, vec_cat, vec_ni, vec_weight, vec_z2, Vq, q,
              &s_vec.vector);

        // compute S
        time_start = clock();
        cPar.CalcS(mapRS2wA, mapRS2wK, W, A, K, &S_mat.matrix, &Svar_mat.matrix,
                   &s_vec.vector);
        cPar.time_G += (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);
        if (cPar.error == true) {
          cout << "error! fail to calculate the S matrix. " << endl;
          return;
        }

        // compute vc estimates
        CalcVCss(Vq, &S_mat.matrix, &Svar_mat.matrix, q, &s_vec.vector,
                 cPar.ni_study, cPar.v_pve, cPar.v_se_pve, cPar.pve_total,
                 cPar.se_pve_total, cPar.v_sigma2, cPar.v_se_sigma2,
                 cPar.v_enrich, cPar.v_se_enrich);
        assert(!has_nan(cPar.v_se_pve));
      }


      gsl_vector_set(s, cPar.n_vc, cPar.ni_test);

      cPar.WriteMatrix(S, "S");
      cPar.WriteMatrix(Vq, "Vq");
      cPar.WriteVector(q, "q");
      cPar.WriteVector(s, "size");

      gsl_matrix_safe_free(S);
      gsl_matrix_safe_free(Vq);
      gsl_vector_safe_free(q);
      gsl_vector_safe_free(s);

      gsl_matrix_safe_free(A);
      gsl_matrix_safe_free(K);
      gsl_vector_safe_free(y);
      gsl_matrix_safe_free(W);
    } else if (!cPar.file_study.empty() || !cPar.file_mstudy.empty()) {
      if (!cPar.file_study.empty()) {
        string sfile = cPar.file_study + ".size.txt";
        CountFileLines(sfile, cPar.n_vc);
      } else {
        string file_name;
        igzstream infile(cPar.file_mstudy.c_str(), igzstream::in);
        if (!infile) {
          cout << "error! fail to open mstudy file: " << cPar.file_study
               << endl;
          return;
        }

        safeGetline(infile, file_name);

        infile.clear();
        infile.close();

        string sfile = file_name + ".size.txt";
        CountFileLines(sfile, cPar.n_vc);
      }

      cPar.n_vc = cPar.n_vc - 1;

      gsl_matrix *S = gsl_matrix_safe_alloc(2 * cPar.n_vc, cPar.n_vc);
      gsl_matrix *Vq = gsl_matrix_safe_alloc(cPar.n_vc, cPar.n_vc);
      // gsl_matrix *V=gsl_matrix_safe_alloc (cPar.n_vc+1,
      // (cPar.n_vc*(cPar.n_vc+1))/2*(cPar.n_vc+1) );
      // gsl_matrix *Vslope=gsl_matrix_safe_alloc (n_lines+1,
      // (n_lines*(n_lines+1))/2*(n_lines+1) );
      gsl_vector *q = gsl_vector_safe_alloc(cPar.n_vc);
      gsl_vector *s_study = gsl_vector_safe_alloc(cPar.n_vc);
      gsl_vector *s_ref = gsl_vector_safe_alloc(cPar.n_vc);
      gsl_vector *s = gsl_vector_safe_alloc(cPar.n_vc + 1);

      gsl_matrix_set_zero(S);
      gsl_matrix_view S_mat =
          gsl_matrix_submatrix(S, 0, 0, cPar.n_vc, cPar.n_vc);
      gsl_matrix_view Svar_mat =
          gsl_matrix_submatrix(S, cPar.n_vc, 0, cPar.n_vc, cPar.n_vc);

      gsl_matrix_set_zero(Vq);
      // gsl_matrix_set_zero(V);
      // gsl_matrix_set_zero(Vslope);
      gsl_vector_set_zero(q);
      gsl_vector_set_zero(s_study);
      gsl_vector_set_zero(s_ref);

      if (!cPar.file_study.empty()) {
        ReadFile_study(cPar.file_study, Vq, q, s_study, cPar.ni_study);
      } else {
        ReadFile_mstudy(cPar.file_mstudy, Vq, q, s_study, cPar.ni_study);
      }

      if (!cPar.file_ref.empty()) {
        ReadFile_ref(cPar.file_ref, &S_mat.matrix, &Svar_mat.matrix, s_ref,
                     cPar.ni_ref);
      } else {
        ReadFile_mref(cPar.file_mref, &S_mat.matrix, &Svar_mat.matrix, s_ref,
                      cPar.ni_ref);
      }

      cout << "## number of variance components = " << cPar.n_vc << endl;
      cout << "## number of individuals in the sample = " << cPar.ni_study
           << endl;
      cout << "## number of individuals in the reference = " << cPar.ni_ref
           << endl;

      CalcVCss(Vq, &S_mat.matrix, &Svar_mat.matrix, q, s_study, cPar.ni_study,
               cPar.v_pve, cPar.v_se_pve, cPar.pve_total, cPar.se_pve_total,
               cPar.v_sigma2, cPar.v_se_sigma2, cPar.v_enrich,
               cPar.v_se_enrich);
      assert(!has_nan(cPar.v_se_pve));

      gsl_vector_view s_sub = gsl_vector_subvector(s, 0, cPar.n_vc);
      gsl_vector_safe_memcpy(&s_sub.vector, s_ref);
      gsl_vector_set(s, cPar.n_vc, cPar.ni_ref);

      cPar.WriteMatrix(S, "S");
      cPar.WriteMatrix(Vq, "Vq");
      cPar.WriteVector(q, "q");
      cPar.WriteVector(s, "size");

      gsl_matrix_safe_free(S);
      gsl_matrix_safe_free(Vq);
      // gsl_matrix_safe_free (V);
      // gsl_matrix_safe_free (Vslope);
      gsl_vector_safe_free(q);
      gsl_vector_safe_free(s_study);
      gsl_vector_safe_free(s_ref);
      gsl_vector_safe_free(s);
    } else {
      gsl_matrix *Y = gsl_matrix_safe_alloc(cPar.ni_test, cPar.n_ph);
      gsl_matrix *W = gsl_matrix_safe_alloc(Y->size1, cPar.n_cvt);
      gsl_matrix *G = gsl_matrix_safe_alloc(Y->size1, Y->size1 * cPar.n_vc);

      // set covariates matrix W and phenotype matrix Y
      // an intercept should be included in W,
      cPar.CopyCvtPhen(W, Y, 0);

      // read kinship matrices
      if (!(cPar.file_mk).empty()) {
        ReadFile_mk(cPar.file_mk, cPar.indicator_idv, cPar.mapID2num,
                    cPar.k_mode, cPar.error, G);
        if (cPar.error == true) {
          cout << "error! fail to read kinship/relatedness file. " << endl;
          return;
        }

        // center matrix G, and obtain v_traceG
        double d = 0;
        (cPar.v_traceG).clear();
        for (size_t i = 0; i < cPar.n_vc; i++) {
          gsl_matrix_view G_sub =
              gsl_matrix_submatrix(G, 0, i * G->size1, G->size1, G->size1);
          CenterMatrix(&G_sub.matrix);
          d = 0;
          for (size_t j = 0; j < G->size1; j++) {
            d += gsl_matrix_get(&G_sub.matrix, j, j);
          }
          d /= (double)G->size1;
          (cPar.v_traceG).push_back(d);
        }
      } else if (!(cPar.file_kin).empty()) {
        ReadFile_kin(cPar.file_kin, cPar.indicator_idv, cPar.mapID2num,
                     cPar.k_mode, cPar.error, G);
        if (cPar.error == true) {
          cout << "error! fail to read kinship/relatedness file. " << endl;
          return;
        }

        // center matrix G
        CenterMatrix(G);
        validate_K(G);

        (cPar.v_traceG).clear();
        double d = 0;
        for (size_t j = 0; j < G->size1; j++) {
          d += gsl_matrix_get(G, j, j);
        }
        d /= (double)G->size1;
        (cPar.v_traceG).push_back(d);
      }
      // fit multiple variance components
      if (cPar.n_ph == 1) {
        //		  if (cPar.n_vc==1) {
        //		  } else {
        gsl_vector_view Y_col = gsl_matrix_column(Y, 0);
        VC cVc;
        cVc.CopyFromParam(cPar);
        if (cPar.a_mode == 61) {
          cVc.CalcVChe(G, W, &Y_col.vector);
        } else if (cPar.a_mode == 62) {
          cVc.CalcVCreml(cPar.noconstrain, G, W, &Y_col.vector);
        } else {
          cVc.CalcVCacl(G, W, &Y_col.vector);
        }
        cVc.CopyToParam(cPar);
        // obtain pve from sigma2
        // obtain se_pve from se_sigma2

        //}
      }
    }
  }

  // compute confidence intervals with additional summary statistics
  // we do not check the sign of z-scores here, but they have to be matched with
  // the genotypes
  if (cPar.a_mode == 66 || cPar.a_mode == 67) {
    // read reference file first
    gsl_matrix *S = gsl_matrix_safe_alloc(cPar.n_vc, cPar.n_vc);
    gsl_matrix *Svar = gsl_matrix_safe_alloc(cPar.n_vc, cPar.n_vc);
    gsl_vector *s_ref = gsl_vector_safe_alloc(cPar.n_vc);

    gsl_matrix_set_zero(S);
    gsl_matrix_set_zero(Svar);
    gsl_vector_set_zero(s_ref);

    if (!cPar.file_ref.empty()) {
      ReadFile_ref(cPar.file_ref, S, Svar, s_ref, cPar.ni_ref);
    } else {
      ReadFile_mref(cPar.file_mref, S, Svar, s_ref, cPar.ni_ref);
    }

    // need to obtain a common set of SNPs between beta file and the genotype
    // file; these are saved in mapRS2wA and mapRS2wK
    // normalize the weight in mapRS2wK to have an average of one; each element
    // of mapRS2wA is 1
    set<string> setSnps_beta;
    ReadFile_snps_header(cPar.file_beta, setSnps_beta);

    // obtain the weights for wA, which contains the SNP weights for SNPs used
    // in the model
    map<string, double> mapRS2wK;
    cPar.ObtainWeight(setSnps_beta, mapRS2wK);

    // set up matrices and vector
    gsl_matrix *Xz = gsl_matrix_safe_alloc(cPar.ni_test, cPar.n_vc);
    gsl_matrix *XWz = gsl_matrix_safe_alloc(cPar.ni_test, cPar.n_vc);
    gsl_matrix *XtXWz =
        gsl_matrix_safe_alloc(mapRS2wK.size(), cPar.n_vc * cPar.n_vc);
    gsl_vector *w = gsl_vector_safe_alloc(mapRS2wK.size());
    gsl_vector *w1 = gsl_vector_safe_alloc(mapRS2wK.size());
    gsl_vector *z = gsl_vector_safe_alloc(mapRS2wK.size());
    gsl_vector *s_vec = gsl_vector_safe_alloc(cPar.n_vc);

    vector<size_t> vec_cat, vec_size;
    vector<double> vec_z;

    map<string, double> mapRS2z, mapRS2wA;
    map<string, string> mapRS2A1;
    string file_str;

    // update s_vec, the number of snps in each category
    for (size_t i = 0; i < cPar.n_vc; i++) {
      vec_size.push_back(0);
    }

    for (map<string, double>::const_iterator it = mapRS2wK.begin();
         it != mapRS2wK.end(); ++it) {
      vec_size[cPar.mapRS2cat[it->first]]++;
    }

    for (size_t i = 0; i < cPar.n_vc; i++) {
      gsl_vector_set(s_vec, i, vec_size[i]);
    }

    // update mapRS2wA using v_pve and s_vec
    if (cPar.a_mode == 66) {
      for (map<string, double>::const_iterator it = mapRS2wK.begin();
           it != mapRS2wK.end(); ++it) {
        mapRS2wA[it->first] = 1;
      }
    } else {
      cPar.UpdateWeight(0, mapRS2wK, cPar.ni_test, s_vec, mapRS2wA);
    }

    // read in z-scores based on allele 0, and save that into a vector
    ReadFile_beta(cPar.file_beta, mapRS2wA, mapRS2A1, mapRS2z);

    // update snp indicator, save weights to w, save z-scores to vec_z, save
    // category label to vec_cat
    // sign of z is determined by matching alleles
    cPar.UpdateSNPnZ(mapRS2wA, mapRS2A1, mapRS2z, w, z, vec_cat);

    // compute an n by k matrix of X_iWz
    cout << "Calculating Xz ... " << endl;

    gsl_matrix_set_zero(Xz);
    gsl_vector_set_all(w1, 1);

    if (!cPar.file_bfile.empty()) {
      file_str = cPar.file_bfile + ".bed";
      PlinkXwz(file_str, cPar.d_pace, cPar.indicator_idv, cPar.indicator_snp,
               vec_cat, w1, z, 0, Xz);
    } else if (!cPar.file_geno.empty()) {
      BimbamXwz(cPar.file_geno, cPar.d_pace, cPar.indicator_idv,
                cPar.indicator_snp, vec_cat, w1, z, 0, Xz);
    } else if (!cPar.file_mbfile.empty()) {
      MFILEXwz(1, cPar.file_mbfile, cPar.d_pace, cPar.indicator_idv,
               cPar.mindicator_snp, vec_cat, w1, z, Xz);
    } else if (!cPar.file_mgeno.empty()) {
      MFILEXwz(0, cPar.file_mgeno, cPar.d_pace, cPar.indicator_idv,
               cPar.mindicator_snp, vec_cat, w1, z, Xz);
    }
    if (cPar.a_mode == 66) {
      gsl_matrix_safe_memcpy(XWz, Xz);
    } else if (cPar.a_mode == 67) {
      cout << "Calculating XWz ... " << endl;

      gsl_matrix_set_zero(XWz);

      if (!cPar.file_bfile.empty()) {
        file_str = cPar.file_bfile + ".bed";
        PlinkXwz(file_str, cPar.d_pace, cPar.indicator_idv, cPar.indicator_snp,
                 vec_cat, w, z, 0, XWz);
      } else if (!cPar.file_geno.empty()) {
        BimbamXwz(cPar.file_geno, cPar.d_pace, cPar.indicator_idv,
                  cPar.indicator_snp, vec_cat, w, z, 0, XWz);
      } else if (!cPar.file_mbfile.empty()) {
        MFILEXwz(1, cPar.file_mbfile, cPar.d_pace, cPar.indicator_idv,
                 cPar.mindicator_snp, vec_cat, w, z, XWz);
      } else if (!cPar.file_mgeno.empty()) {
        MFILEXwz(0, cPar.file_mgeno, cPar.d_pace, cPar.indicator_idv,
                 cPar.mindicator_snp, vec_cat, w, z, XWz);
      }
    }
    // compute an p by k matrix of X_j^TWX_iWz
    cout << "Calculating XtXWz ... " << endl;
    gsl_matrix_set_zero(XtXWz);

    if (!cPar.file_bfile.empty()) {
      file_str = cPar.file_bfile + ".bed";
      PlinkXtXwz(file_str, cPar.d_pace, cPar.indicator_idv, cPar.indicator_snp,
                 XWz, 0, XtXWz);
    } else if (!cPar.file_geno.empty()) {
      BimbamXtXwz(cPar.file_geno, cPar.d_pace, cPar.indicator_idv,
                  cPar.indicator_snp, XWz, 0, XtXWz);
    } else if (!cPar.file_mbfile.empty()) {
      MFILEXtXwz(1, cPar.file_mbfile, cPar.d_pace, cPar.indicator_idv,
                 cPar.mindicator_snp, XWz, XtXWz);
    } else if (!cPar.file_mgeno.empty()) {
      MFILEXtXwz(0, cPar.file_mgeno, cPar.d_pace, cPar.indicator_idv,
                 cPar.mindicator_snp, XWz, XtXWz);
    }
    // compute confidence intervals
    CalcCIss(Xz, XWz, XtXWz, S, Svar, w, z, s_vec, vec_cat, cPar.v_pve,
             cPar.v_se_pve, cPar.pve_total, cPar.se_pve_total, cPar.v_sigma2,
             cPar.v_se_sigma2, cPar.v_enrich, cPar.v_se_enrich);
    assert(!has_nan(cPar.v_se_pve));

    gsl_matrix_safe_free(S);
    gsl_matrix_safe_free(Svar);
    gsl_vector_safe_free(s_ref);

    gsl_matrix_safe_free(Xz);
    gsl_matrix_safe_free(XWz);
    gsl_matrix_safe_free(XtXWz);
    gsl_vector_safe_free(w);
    gsl_vector_safe_free(w1);
    gsl_vector_safe_free(z);
    gsl_vector_safe_free(s_vec);
  }

  // LMM or mvLMM or Eigen-Decomposition
  if (cPar.a_mode == M_LMM1 || cPar.a_mode == M_LMM2 || cPar.a_mode == M_LMM3 ||
      cPar.a_mode == M_LMM4 || cPar.a_mode == M_LMM5 ||
      cPar.a_mode == M_EIGEN) { // Fit LMM or mvLMM or eigen
    gsl_matrix *Y = gsl_matrix_safe_alloc(cPar.ni_test, cPar.n_ph);
    enforce_msg(Y, "allocate Y"); // just to be sure
    gsl_matrix *W = gsl_matrix_safe_alloc(Y->size1, cPar.n_cvt);
    gsl_matrix *B = gsl_matrix_safe_alloc(Y->size2, W->size2); // B is a d by c
                                                          // matrix
    gsl_matrix *se_B = gsl_matrix_safe_alloc(Y->size2, W->size2);
    gsl_matrix *G = gsl_matrix_safe_alloc(Y->size1, Y->size1);
    gsl_matrix *U = gsl_matrix_safe_alloc(Y->size1, Y->size1);
    gsl_matrix *UtW = gsl_matrix_calloc(Y->size1, W->size2);
    gsl_matrix *UtY = gsl_matrix_calloc(Y->size1, Y->size2);
    gsl_vector *eval = gsl_vector_calloc(Y->size1);
    gsl_vector *env = gsl_vector_safe_alloc(Y->size1);
    gsl_vector *weight = gsl_vector_safe_alloc(Y->size1);
    assert_issue(is_issue(26), UtY->data[0] == 0.0);

    // set covariates matrix W and phenotype matrix Y
    // an intercept should be included in W,
    cPar.CopyCvtPhen(W, Y, 0);
    if (!cPar.file_gxe.empty()) {
      cPar.CopyGxe(env);
    }

    // read relatedness matrix G
    if (!(cPar.file_kin).empty()) {
      ReadFile_kin(cPar.file_kin, cPar.indicator_idv, cPar.mapID2num,
                   cPar.k_mode, cPar.error, G);
      if (cPar.error == true) {
        cout << "error! fail to read kinship/relatedness file. " << endl;
        return;
      }

      // center matrix G
      CenterMatrix(G);
      validate_K(G);

      // is residual weights are provided, then
      if (!cPar.file_weight.empty()) {
        cPar.CopyWeight(weight);
        double d, wi, wj;
        for (size_t i = 0; i < G->size1; i++) {
          wi = gsl_vector_get(weight, i);
          for (size_t j = i; j < G->size2; j++) {
            wj = gsl_vector_get(weight, j);
            d = gsl_matrix_get(G, i, j);
            if (wi <= 0 || wj <= 0) {
              d = 0;
            } else {
              d /= safe_sqrt(wi * wj);
            }
            gsl_matrix_set(G, i, j, d);
            if (j != i) {
              gsl_matrix_set(G, j, i, d);
            }
          }
        }
      }

      // eigen-decomposition and calculate trace_G
      cout << "Start Eigen-Decomposition..." << endl;
      time_start = clock();

      if (cPar.a_mode == M_EIGEN) {
        cPar.trace_G = EigenDecomp_Zeroed(G, U, eval, 1);
      } else {
        cPar.trace_G = EigenDecomp_Zeroed(G, U, eval, 0);
      }
      // write(eval,"eval");

      if (!cPar.file_weight.empty()) {
        double wi;
        for (size_t i = 0; i < U->size1; i++) {
          wi = gsl_vector_get(weight, i);
          if (wi <= 0) {
            wi = 0;
          } else {
            wi = safe_sqrt(wi);
          }
          gsl_vector_view Urow = gsl_matrix_row(U, i);
          gsl_vector_scale(&Urow.vector, wi);
        }
      }

      cPar.time_eigen =
          (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);
    } else {
      ReadFile_eigenU(cPar.file_ku, cPar.error, U);
      if (cPar.error == true) {
        cout << "error! fail to read the U file. " << endl;
        return;
      }

      ReadFile_eigenD(cPar.file_kd, cPar.error, eval);
      if (cPar.error == true) {
        cout << "error! fail to read the D file. " << endl;
        return;
      }

      cPar.trace_G = 0.0;
      for (size_t i = 0; i < eval->size; i++) {
        if (gsl_vector_get(eval, i) < 1e-10) {
          gsl_vector_set(eval, i, 0);
        }
        cPar.trace_G += gsl_vector_get(eval, i);
      }
      cPar.trace_G /= (double)eval->size;
    }
    // write(eval,"eval2");

    if (cPar.a_mode == M_EIGEN) {
      cPar.WriteMatrix(U, "eigenU");
      cPar.WriteVector(eval, "eigenD");
    } else if (!cPar.file_gene.empty()) {
      // calculate UtW and Uty
      CalcUtX(U, W, UtW);
      CalcUtX(U, Y, UtY);

      assert_issue(is_issue(26), ROUND(UtY->data[0]) == -16.6143);

      LMM cLmm;
      cLmm.CopyFromParam(cPar);

      gsl_vector_view Y_col = gsl_matrix_column(Y, 0);
      gsl_vector_view UtY_col = gsl_matrix_column(UtY, 0);

      cLmm.AnalyzeGene(U, eval, UtW, &UtY_col.vector, W,
                       &Y_col.vector); // y is the predictor, not the phenotype

      cLmm.WriteFiles();
      cLmm.CopyToParam(cPar);
    } else {
      // calculate UtW and Uty
      CalcUtX(U, W, UtW);
      CalcUtX(U, Y, UtY);
      assert_issue(is_issue(26), ROUND(UtY->data[0]) == -16.6143);

      // calculate REMLE/MLE estimate and pve for univariate model
      if (cPar.n_ph == 1) { // one phenotype
        gsl_vector_view beta = gsl_matrix_row(B, 0);
        gsl_vector_view se_beta = gsl_matrix_row(se_B, 0);
        gsl_vector_view UtY_col = gsl_matrix_column(UtY, 0);

        assert_issue(is_issue(26), ROUND(UtY->data[0]) == -16.6143);

        CalcLambda('L', eval, UtW, &UtY_col.vector, cPar.l_min, cPar.l_max,
                   cPar.n_region, cPar.l_mle_null, cPar.logl_mle_H0);
        assert(!isnan(UtY->data[0]));

        CalcLmmVgVeBeta(eval, UtW, &UtY_col.vector, cPar.l_mle_null,
                        cPar.vg_mle_null, cPar.ve_mle_null, &beta.vector,
                        &se_beta.vector);

        assert(!isnan(UtY->data[0]));

        cPar.beta_mle_null.clear();
        cPar.se_beta_mle_null.clear();
        assert(!isnan(B->data[0]));
        assert(!isnan(se_B->data[0]));
        for (size_t i = 0; i < B->size2; i++) {
          cPar.beta_mle_null.push_back(gsl_matrix_get(B, 0, i));
          cPar.se_beta_mle_null.push_back(gsl_matrix_get(se_B, 0, i));
        }
        assert(!isnan(UtY->data[0]));
        assert(!isnan(cPar.beta_mle_null.front()));
        assert(!isnan(cPar.se_beta_mle_null.front()));

        // the following functions do not modify eval
        CalcLambda('R', eval, UtW, &UtY_col.vector, cPar.l_min, cPar.l_max,
                   cPar.n_region, cPar.l_remle_null, cPar.logl_remle_H0);
        CalcLmmVgVeBeta(eval, UtW, &UtY_col.vector, cPar.l_remle_null,
                        cPar.vg_remle_null, cPar.ve_remle_null, &beta.vector,
                        &se_beta.vector);

        cPar.beta_remle_null.clear();
        cPar.se_beta_remle_null.clear();
        assert(!isnan(B->data[0]));
        assert(!isnan(se_B->data[0]));

        for (size_t i = 0; i < B->size2; i++) {
          cPar.beta_remle_null.push_back(gsl_matrix_get(B, 0, i));
          cPar.se_beta_remle_null.push_back(gsl_matrix_get(se_B, 0, i));
        }

        CalcPve(eval, UtW, &UtY_col.vector, cPar.l_remle_null, cPar.trace_G,
                cPar.pve_null, cPar.pve_se_null);
        cPar.PrintSummary();

        // calculate and output residuals
        if (cPar.a_mode == M_LMM5) {
          gsl_vector *Utu_hat = gsl_vector_safe_alloc(Y->size1);
          gsl_vector *Ute_hat = gsl_vector_safe_alloc(Y->size1);
          gsl_vector *u_hat = gsl_vector_safe_alloc(Y->size1);
          gsl_vector *e_hat = gsl_vector_safe_alloc(Y->size1);
          gsl_vector *y_hat = gsl_vector_safe_alloc(Y->size1);

          // obtain Utu and Ute
          gsl_vector_safe_memcpy(y_hat, &UtY_col.vector);
          gsl_blas_dgemv(CblasNoTrans, -1.0, UtW, &beta.vector, 1.0, y_hat);

          double d, u, e;
          for (size_t i = 0; i < eval->size; i++) {
            d = gsl_vector_get(eval, i);
            u = cPar.l_remle_null * d / (cPar.l_remle_null * d + 1.0) *
                gsl_vector_get(y_hat, i);
            e = 1.0 / (cPar.l_remle_null * d + 1.0) * gsl_vector_get(y_hat, i);
            gsl_vector_set(Utu_hat, i, u);
            gsl_vector_set(Ute_hat, i, e);
          }

          // obtain u and e
          gsl_blas_dgemv(CblasNoTrans, 1.0, U, Utu_hat, 0.0, u_hat);
          gsl_blas_dgemv(CblasNoTrans, 1.0, U, Ute_hat, 0.0, e_hat);

          // output residuals
          cPar.WriteVector(u_hat, "residU");
          cPar.WriteVector(e_hat, "residE");

          gsl_vector_safe_free(u_hat);
          gsl_vector_safe_free(e_hat);
          gsl_vector_safe_free(y_hat);
        }
      }

      // Fit LMM or mvLMM (w. LOCO)
      if (cPar.a_mode == 1 || cPar.a_mode == 2 || cPar.a_mode == 3 ||
          cPar.a_mode == 4) {
        if (cPar.n_ph == 1) {
          LMM cLmm;
          cLmm.CopyFromParam(cPar);

          // if (is_check_mode()) disable_segfpe(); // disable fast NaN checking for now

          gsl_vector_view Y_col = gsl_matrix_column(Y, 0);
          gsl_vector_view UtY_col = gsl_matrix_column(UtY, 0);

          if (!cPar.file_bfile.empty()) {
            // PLINK analysis
            if (cPar.file_gxe.empty()) {
              cLmm.AnalyzePlink(U, eval, UtW, &UtY_col.vector, W,
                                &Y_col.vector, cPar.setGWASnps);
            }
          }
          
          cLmm.WriteFiles();
          cLmm.CopyToParam(cPar);
        } else {
          MVLMM cMvlmm;
          cMvlmm.CopyFromParam(cPar);

          // if (is_check_mode()) disable_segfpe(); // disable fast NaN checking

          // write(eval,"eval3");

          if (!cPar.file_bfile.empty()) {
            if (cPar.file_gxe.empty()) {
              cMvlmm.AnalyzePlink(U, eval, UtW, UtY);
            } 
          } 

          cMvlmm.WriteFiles();
          cMvlmm.CopyToParam(cPar);
        }
      }
    }

    // release all matrices and vectors
    gsl_matrix_safe_free(Y);
    gsl_matrix_safe_free(W);
    gsl_matrix_warn_free(B); // sometimes unused
    gsl_matrix_warn_free(se_B);
    gsl_matrix_warn_free(G);
    gsl_matrix_safe_free(U);
    gsl_matrix_safe_free(UtW);
    gsl_matrix_safe_free(UtY);
    gsl_vector_safe_free(eval);
    gsl_vector_free(env); // sometimes unused
  }

  // BSLMM
  if (cPar.a_mode == 11 || cPar.a_mode == 12 || cPar.a_mode == 13) {
    gsl_vector *y = gsl_vector_safe_alloc(cPar.ni_test);
    gsl_matrix *W = gsl_matrix_safe_alloc(y->size, cPar.n_cvt);
    gsl_matrix *G = gsl_matrix_safe_alloc(y->size, y->size);
    gsl_matrix *UtX = gsl_matrix_safe_alloc(y->size, cPar.ns_test);

    // set covariates matrix W and phenotype vector y
    // an intercept should be included in W,
    cPar.CopyCvtPhen(W, y, 0);

    // center y, even for case/control data
    cPar.pheno_mean = CenterVector(y);

    // run bvsr if rho==1
    if (cPar.rho_min == 1 && cPar.rho_max == 1) {
      // read genotypes X (not UtX)
      cPar.ReadGenotypes(UtX, G, false);

      // perform BSLMM analysis
      BSLMM cBslmm;
      cBslmm.CopyFromParam(cPar);
      time_start = clock();
      cBslmm.MCMC(UtX, y);
      cPar.time_opt = (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);
      cBslmm.CopyToParam(cPar);
      // else, if rho!=1
    } else {
      gsl_matrix *U = gsl_matrix_safe_alloc(y->size, y->size);
      gsl_vector *eval = gsl_vector_safe_alloc(y->size);
      gsl_matrix *UtW = gsl_matrix_safe_alloc(y->size, W->size2);
      gsl_vector *Uty = gsl_vector_safe_alloc(y->size);

      // read relatedness matrix G
      if (!(cPar.file_kin).empty()) {
        cPar.ReadGenotypes(UtX, G, false);

        // read relatedness matrix G
        ReadFile_kin(cPar.file_kin, cPar.indicator_idv, cPar.mapID2num,
                     cPar.k_mode, cPar.error, G);
        if (cPar.error == true) {
          cout << "error! fail to read kinship/relatedness file. " << endl;
          return;
        }

        // center matrix G
        CenterMatrix(G);
        validate_K(G);
      } else {
        cPar.ReadGenotypes(UtX, G, true);
      }

      // eigen-decomposition and calculate trace_G
      cout << "Start Eigen-Decomposition..." << endl;
      time_start = clock();

      cPar.trace_G = EigenDecomp_Zeroed(G, U, eval, 0);

      cPar.time_eigen =
          (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

      // calculate UtW and Uty
      CalcUtX(U, W, UtW);
      CalcUtX(U, y, Uty);

      // calculate REMLE/MLE estimate and pve
      CalcLambda('L', eval, UtW, Uty, cPar.l_min, cPar.l_max, cPar.n_region,
                 cPar.l_mle_null, cPar.logl_mle_H0);
      CalcLambda('R', eval, UtW, Uty, cPar.l_min, cPar.l_max, cPar.n_region,
                 cPar.l_remle_null, cPar.logl_remle_H0);
      CalcPve(eval, UtW, Uty, cPar.l_remle_null, cPar.trace_G, cPar.pve_null,
              cPar.pve_se_null);

      cPar.PrintSummary();

      // Creat and calcualte UtX, use a large memory
      cout << "Calculating UtX..." << endl;
      time_start = clock();
      CalcUtX(U, UtX);
      cPar.time_UtX = (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

      // perform BSLMM or BSLMMDAP analysis
      if (cPar.a_mode == 11 || cPar.a_mode == 12 || cPar.a_mode == 13) {
        BSLMM cBslmm;
        cBslmm.CopyFromParam(cPar);
        time_start = clock();
        if (cPar.a_mode == 12) { // ridge regression
          cBslmm.RidgeR(U, UtX, Uty, eval, cPar.l_remle_null);
        } else { // Run MCMC
          cBslmm.MCMC(U, UtX, Uty, eval, y);
        }
        cPar.time_opt =
            (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);
        cBslmm.CopyToParam(cPar);
      } else {
      }

      // release all matrices and vectors
      gsl_matrix_safe_free(G);
      gsl_matrix_safe_free(U);
      gsl_matrix_safe_free(UtW);
      gsl_vector_safe_free(eval);
      gsl_vector_safe_free(Uty);
    }
    gsl_matrix_safe_free(W);
    gsl_vector_safe_free(y);
    gsl_matrix_safe_free(UtX);
  }

  // BSLMM-DAP
  if (cPar.a_mode == 14 || cPar.a_mode == 15 || cPar.a_mode == 16) {
    if (cPar.a_mode == 14) {
      gsl_vector *y = gsl_vector_safe_alloc(cPar.ni_test);
      gsl_matrix *W = gsl_matrix_safe_alloc(y->size, cPar.n_cvt);
      gsl_matrix *G = gsl_matrix_safe_alloc(y->size, y->size);
      gsl_matrix *UtX = gsl_matrix_safe_alloc(y->size, cPar.ns_test);

      // set covariates matrix W and phenotype vector y
      // an intercept should be included in W,
      cPar.CopyCvtPhen(W, y, 0);

      // center y, even for case/control data
      cPar.pheno_mean = CenterVector(y);

      // run bvsr if rho==1
      if (cPar.rho_min == 1 && cPar.rho_max == 1) {
        // read genotypes X (not UtX)
        cPar.ReadGenotypes(UtX, G, false);

        // perform BSLMM analysis
        BSLMM cBslmm;
        cBslmm.CopyFromParam(cPar);
        time_start = clock();
        cBslmm.MCMC(UtX, y);
        cPar.time_opt =
            (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);
        cBslmm.CopyToParam(cPar);
        // else, if rho!=1
      } else {
        gsl_matrix *U = gsl_matrix_safe_alloc(y->size, y->size);
        gsl_vector *eval = gsl_vector_safe_alloc(y->size);
        gsl_matrix *UtW = gsl_matrix_safe_alloc(y->size, W->size2);
        gsl_vector *Uty = gsl_vector_safe_alloc(y->size);

        // read relatedness matrix G
        if (!(cPar.file_kin).empty()) {
          cPar.ReadGenotypes(UtX, G, false);

          // read relatedness matrix G
          ReadFile_kin(cPar.file_kin, cPar.indicator_idv, cPar.mapID2num,
                       cPar.k_mode, cPar.error, G);
          if (cPar.error == true) {
            cout << "error! fail to read kinship/relatedness file. " << endl;
            return;
          }

          // center matrix G
          CenterMatrix(G);
          validate_K(G);

        } else {
          cPar.ReadGenotypes(UtX, G, true);
        }

        // eigen-decomposition and calculate trace_G
        cout << "Start Eigen-Decomposition..." << endl;
        time_start = clock();
        cPar.trace_G = EigenDecomp_Zeroed(G, U, eval, 0);
        cPar.time_eigen =
            (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

        // calculate UtW and Uty
        CalcUtX(U, W, UtW);
        CalcUtX(U, y, Uty);

        // calculate REMLE/MLE estimate and pve
        CalcLambda('L', eval, UtW, Uty, cPar.l_min, cPar.l_max, cPar.n_region,
                   cPar.l_mle_null, cPar.logl_mle_H0);
        CalcLambda('R', eval, UtW, Uty, cPar.l_min, cPar.l_max, cPar.n_region,
                   cPar.l_remle_null, cPar.logl_remle_H0);
        CalcPve(eval, UtW, Uty, cPar.l_remle_null, cPar.trace_G, cPar.pve_null,
                cPar.pve_se_null);

        cPar.PrintSummary();

        // Creat and calcualte UtX, use a large memory
        cout << "Calculating UtX..." << endl;
        time_start = clock();
        CalcUtX(U, UtX);
        cPar.time_UtX =
            (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);

        // perform analysis; assume X and y are already centered
        BSLMMDAP cBslmmDap;
        cBslmmDap.CopyFromParam(cPar);
        time_start = clock();
        cBslmmDap.DAP_CalcBF(U, UtX, Uty, eval, y);
        cPar.time_opt =
            (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);
        cBslmmDap.CopyToParam(cPar);

        // release all matrices and vectors
        gsl_matrix_safe_free(G);
        gsl_matrix_safe_free(U);
        gsl_matrix_safe_free(UtW);
        gsl_vector_safe_free(eval);
        gsl_vector_safe_free(Uty);
      }

      gsl_matrix_safe_free(W);
      gsl_vector_safe_free(y);
      gsl_matrix_safe_free(UtX);
    } else if (cPar.a_mode == 15) {
      // perform EM algorithm and estimate parameters
      vector<string> vec_rs;
      vector<double> vec_sa2, vec_sb2, wab;
      vector<vector<vector<double>>> BF;

      // read hyp and bf files (functions defined in BSLMMDAP)
      ReadFile_hyb(cPar.file_hyp, vec_sa2, vec_sb2, wab);
      ReadFile_bf(cPar.file_bf, vec_rs, BF);

      cPar.ns_test = vec_rs.size();
      if (wab.size() != BF[0][0].size()) {
        cout << "error! hyp and bf files dimension do not match" << endl;
      }

      // load annotations
      gsl_matrix *Ac = NULL;
      gsl_matrix_int *Ad = NULL;
      gsl_vector_int *dlevel = NULL;
      size_t kc, kd;
      if (!cPar.file_cat.empty()) {
        ReadFile_cat(cPar.file_cat, vec_rs, Ac, Ad, dlevel, kc, kd);
      } else {
        kc = 0;
        kd = 0;
      }

      cout << "## number of blocks = " << BF.size() << endl;
      cout << "## number of analyzed SNPs/var = " << vec_rs.size() << endl;
      cout << "## grid size for hyperparameters = " << wab.size() << endl;
      cout << "## number of continuous annotations = " << kc << endl;
      cout << "## number of discrete annotations = " << kd << endl;

      // DAP_EstimateHyper (const size_t kc, const size_t kd, const
      // vector<string> &vec_rs, const vector<double> &vec_sa2, const
      // vector<double> &vec_sb2, const vector<double> &wab, const
      // vector<vector<vector<double> > > &BF, gsl_matrix *Ac, gsl_matrix_int
      // *Ad, gsl_vector_int *dlevel);

      // perform analysis
      BSLMMDAP cBslmmDap;
      cBslmmDap.CopyFromParam(cPar);
      time_start = clock();
      cBslmmDap.DAP_EstimateHyper(kc, kd, vec_rs, vec_sa2, vec_sb2, wab, BF, Ac,
                                  Ad, dlevel);
      cPar.time_opt = (clock() - time_start) / (double(CLOCKS_PER_SEC) * 60.0);
      cBslmmDap.CopyToParam(cPar);

      gsl_matrix_safe_free(Ac);
      gsl_matrix_int_free(Ad);
      gsl_vector_int_free(dlevel);
    } else {
      //
    }
  }

  cPar.time_total = (clock() - time_begin) / (double(CLOCKS_PER_SEC) * 60.0);

  return;
}

// #include "Eigen/Dense"

void GEMMA::WriteLog(int argc, char **argv, PARAM &cPar) {
  string file_str;
  file_str = cPar.path_out + "/" + cPar.file_out;
  file_str += ".log.txt";

  ofstream outfile(file_str.c_str(), ofstream::out);
  if (!outfile) {
    cout << "error writing log file: " << file_str.c_str() << endl;
    return;
  }

  outfile << "##" << endl;
  outfile << "## GEMMAX Version    = " << version << " (" << date << ")" << endl;
  outfile << "## Build profile    = " << GEMMA_PROFILE << endl ;
  #ifdef __GNUC__
  outfile << "## GCC version      = " << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__ << endl;
  #endif
  outfile << "## GSL Version      = " << GSL_VERSION << endl;
  // outfile << "## Eigen Version    = " << EIGEN_WORLD_VERSION << "." << EIGEN_MAJOR_VERSION << "." << EIGEN_MINOR_VERSION << endl;
#ifdef OPENBLAS

  #ifndef OPENBLAS_LEGACY
  outfile << "## OpenBlas         =" << OPENBLAS_VERSION << " - " << openblas_get_config() << endl;
  outfile << "##   arch           = " << openblas_get_corename() << endl;
  outfile << "##   threads        = " << openblas_get_num_threads() << endl;
  #else
  outfile << "## OpenBlas         = " << openblas_get_config() << endl;
  #endif
  string* pStr = new string[4] { "sequential", "threaded", "openmp" };
  outfile << "##   parallel type  = " << pStr[openblas_get_parallel()] << endl;
#endif

  outfile << "##" << endl;
  outfile << "## Command Line Input = ";
  for (int i = 0; i < argc; i++) {
    outfile << argv[i] << " ";
  }
  outfile << endl;

  outfile << "##" << endl;
  time_t rawtime;
  time(&rawtime);
  tm *ptm = localtime(&rawtime);

  outfile << "## Date = " << asctime(ptm);

  outfile << "##" << endl;
  outfile << "## Summary Statistics:" << endl;
  if (!cPar.file_cor.empty() || !cPar.file_study.empty() ||
      !cPar.file_mstudy.empty()) {
    outfile << "## number of total individuals in the sample = "
            << cPar.ni_study << endl;
    outfile << "## number of total individuals in the reference = "
            << cPar.ni_ref << endl;
    outfile << "## number of variance components = " << cPar.n_vc << endl;

    outfile << "## pve estimates = ";
    for (size_t i = 0; i < cPar.v_pve.size(); i++) {
      outfile << "  " << cPar.v_pve[i];
    }
    outfile << endl;

    outfile << "## se(pve) = ";
    for (size_t i = 0; i < cPar.v_se_pve.size(); i++) {
      outfile << "  " << cPar.v_se_pve[i];
    }
    outfile << endl;
    assert(!has_nan(cPar.v_se_pve));

    if (cPar.n_vc > 1) {
      outfile << "## total pve = " << cPar.pve_total << endl;
      outfile << "## se(total pve) = " << cPar.se_pve_total << endl;
    }

    outfile << "## sigma2 per snp = ";
    for (size_t i = 0; i < cPar.v_sigma2.size(); i++) {
      outfile << "  " << cPar.v_sigma2[i];
    }
    outfile << endl;

    outfile << "## se(sigma2 per snp) = ";
    for (size_t i = 0; i < cPar.v_se_sigma2.size(); i++) {
      outfile << "  " << cPar.v_se_sigma2[i];
    }
    outfile << endl;

    outfile << "## enrichment = ";
    for (size_t i = 0; i < cPar.v_enrich.size(); i++) {
      outfile << "  " << cPar.v_enrich[i];
    }
    outfile << endl;

    outfile << "## se(enrichment) = ";
    for (size_t i = 0; i < cPar.v_se_enrich.size(); i++) {
      outfile << "  " << cPar.v_se_enrich[i];
    }
    outfile << endl;
  } else if (!cPar.file_beta.empty() &&
             (cPar.a_mode == 61 || cPar.a_mode == 62)) {
    outfile << "## number of total individuals in the sample = "
            << cPar.ni_study << endl;
    outfile << "## number of total individuals in the reference = "
            << cPar.ni_total << endl;
    outfile << "## number of total SNPs/var in the sample = " << cPar.ns_study
            << endl;
    outfile << "## number of total SNPs/var in the reference panel = "
            << cPar.ns_total << endl;
    outfile << "## number of analyzed SNPs/var = " << cPar.ns_test << endl;
    outfile << "## number of variance components = " << cPar.n_vc << endl;
  } else if (!cPar.file_beta.empty() &&
             (cPar.a_mode == 66 || cPar.a_mode == 67)) {
    outfile << "## number of total individuals in the sample = "
            << cPar.ni_total << endl;
    outfile << "## number of total individuals in the reference = "
            << cPar.ni_ref << endl;
    outfile << "## number of total SNPs/var in the sample = " << cPar.ns_total
            << endl;
    outfile << "## number of analyzed SNPs/var = " << cPar.ns_test << endl;
    outfile << "## number of variance components = " << cPar.n_vc << endl;

    outfile << "## pve estimates = ";
    for (size_t i = 0; i < cPar.v_pve.size(); i++) {
      outfile << "  " << cPar.v_pve[i];
    }
    outfile << endl;

    outfile << "## se(pve) = ";
    for (size_t i = 0; i < cPar.v_se_pve.size(); i++) {
      outfile << "  " << cPar.v_se_pve[i];
    }
    outfile << endl;

    if (cPar.n_vc > 1) {
      outfile << "## total pve = " << cPar.pve_total << endl;
      outfile << "## se(total pve) = " << cPar.se_pve_total << endl;
    }

    outfile << "## sigma2 per snp = ";
    for (size_t i = 0; i < cPar.v_sigma2.size(); i++) {
      outfile << "  " << cPar.v_sigma2[i];
    }
    outfile << endl;

    outfile << "## se(sigma2 per snp) = ";
    for (size_t i = 0; i < cPar.v_se_sigma2.size(); i++) {
      outfile << "  " << cPar.v_se_sigma2[i];
    }
    outfile << endl;

    outfile << "## enrichment = ";
    for (size_t i = 0; i < cPar.v_enrich.size(); i++) {
      outfile << "  " << cPar.v_enrich[i];
    }
    outfile << endl;

    outfile << "## se(enrichment) = ";
    for (size_t i = 0; i < cPar.v_se_enrich.size(); i++) {
      outfile << "  " << cPar.v_se_enrich[i];
    }
    outfile << endl;
  } else {
    outfile << "## number of total individuals = " << cPar.ni_total << endl;

    if (cPar.a_mode == 43) {
      outfile << "## number of analyzed individuals = " << cPar.ni_cvt << endl;
      outfile << "## number of individuals with full phenotypes = "
              << cPar.ni_test << endl;
    } else if (cPar.a_mode != 27 && cPar.a_mode != 28) {
      outfile << "## number of analyzed individuals = " << cPar.ni_test << endl;
    }

    outfile << "## number of covariates = " << cPar.n_cvt << endl;
    outfile << "## number of phenotypes = " << cPar.n_ph << endl;
    if (cPar.a_mode == 43) {
      outfile << "## number of observed data = " << cPar.np_obs << endl;
      outfile << "## number of missing data = " << cPar.np_miss << endl;
    }
    if (cPar.a_mode == 25 || cPar.a_mode == 26 || cPar.a_mode == 27 ||
        cPar.a_mode == 28 || cPar.a_mode == 61 || cPar.a_mode == 62 ||
        cPar.a_mode == 63 || cPar.a_mode == 66 || cPar.a_mode == 67) {
      outfile << "## number of variance components = " << cPar.n_vc << endl;
    }

    if (!(cPar.file_gene).empty()) {
      outfile << "## number of total genes = " << cPar.ng_total << endl;
      outfile << "## number of analyzed genes = " << cPar.ng_test << endl;
    } else if (cPar.file_epm.empty()) {
      outfile << "## number of total SNPs/var = " << cPar.ns_total << endl;
      outfile << "## number of analyzed SNPs/var = " << cPar.ns_test << endl;
    } else {
      outfile << "## number of analyzed SNPs/var = " << cPar.ns_test << endl;
    }

    if (cPar.a_mode == 13) {
      outfile << "## number of cases = " << cPar.ni_case << endl;
      outfile << "## number of controls = " << cPar.ni_control << endl;
    }
  }

  if ((cPar.a_mode == 61 || cPar.a_mode == 62 || cPar.a_mode == 63) &&
      cPar.file_cor.empty() && cPar.file_study.empty() &&
      cPar.file_mstudy.empty()) {
    //	        outfile<<"## REMLE log-likelihood in the null model =
    //"<<cPar.logl_remle_H0<<endl;
    if (cPar.n_ph == 1) {
      outfile << "## pve estimates = ";
      for (size_t i = 0; i < cPar.v_pve.size(); i++) {
        outfile << "  " << cPar.v_pve[i];
      }
      outfile << endl;

      outfile << "## se(pve) = ";
      for (size_t i = 0; i < cPar.v_se_pve.size(); i++) {
        outfile << "  " << cPar.v_se_pve[i];
      }
      outfile << endl;

      if (cPar.n_vc > 1) {
        outfile << "## total pve = " << cPar.pve_total << endl;
        outfile << "## se(total pve) = " << cPar.se_pve_total << endl;
      }

      outfile << "## sigma2 estimates = ";
      for (size_t i = 0; i < cPar.v_sigma2.size(); i++) {
        outfile << "  " << cPar.v_sigma2[i];
      }
      outfile << endl;

      outfile << "## se(sigma2) = ";
      for (size_t i = 0; i < cPar.v_se_sigma2.size(); i++) {
        outfile << "  " << cPar.v_se_sigma2[i];
      }
      outfile << endl;

      if (!cPar.file_beta.empty()) {
        outfile << "## enrichment = ";
        for (size_t i = 0; i < cPar.v_enrich.size(); i++) {
          outfile << "  " << cPar.v_enrich[i];
        }
        outfile << endl;

        outfile << "## se(enrichment) = ";
        for (size_t i = 0; i < cPar.v_se_enrich.size(); i++) {
          outfile << "  " << cPar.v_se_enrich[i];
        }
        outfile << endl;
      }
    }
  }

  if (cPar.a_mode == 1 || cPar.a_mode == 2 || cPar.a_mode == 3 ||
      cPar.a_mode == 4 || cPar.a_mode == 5 || cPar.a_mode == 11 ||
      cPar.a_mode == 12 || cPar.a_mode == 13) {
    outfile << "## REMLE log-likelihood in the null model = "
            << cPar.logl_remle_H0 << endl;
    outfile << "## MLE log-likelihood in the null model = " << cPar.logl_mle_H0
            << endl;
    if (cPar.n_ph == 1) {
      // outfile<<"## lambda REMLE estimate in the null (linear mixed) model =
      // "<<cPar.l_remle_null<<endl;
      // outfile<<"## lambda MLE estimate in the null (linear mixed) model =
      // "<<cPar.l_mle_null<<endl;
      outfile << "## pve estimate in the null model = " << cPar.pve_null
              << endl;
      outfile << "## se(pve) in the null model = " << cPar.pve_se_null << endl;
      outfile << "## vg estimate in the null model = " << cPar.vg_remle_null
              << endl;
      outfile << "## ve estimate in the null model = " << cPar.ve_remle_null
              << endl;
      outfile << "## beta estimate in the null model = ";
      for (size_t i = 0; i < cPar.beta_remle_null.size(); i++) {
        outfile << "  " << cPar.beta_remle_null[i];
      }
      outfile << endl;
      outfile << "## se(beta) = ";
      for (size_t i = 0; i < cPar.se_beta_remle_null.size(); i++) {
        outfile << "  " << cPar.se_beta_remle_null[i];
      }
      outfile << endl;

    } else {
      size_t c;
      outfile << "## REMLE estimate for Vg in the null model: " << endl;
      for (size_t i = 0; i < cPar.n_ph; i++) {
        for (size_t j = 0; j <= i; j++) {
          c = (2 * cPar.n_ph - min(i, j) + 1) * min(i, j) / 2 + max(i, j) -
              min(i, j);
          outfile << tab(j) << cPar.Vg_remle_null[c];
        }
        outfile << endl;
      }
      outfile << "## se(Vg): " << endl;
      for (size_t i = 0; i < cPar.n_ph; i++) {
        for (size_t j = 0; j <= i; j++) {
          c = (2 * cPar.n_ph - min(i, j) + 1) * min(i, j) / 2 + max(i, j) -
              min(i, j);
          outfile << tab(j) << safe_sqrt(cPar.VVg_remle_null[c]);
        }
        outfile << endl;
      }
      outfile << "## REMLE estimate for Ve in the null model: " << endl;
      for (size_t i = 0; i < cPar.n_ph; i++) {
        for (size_t j = 0; j <= i; j++) {
          c = (2 * cPar.n_ph - min(i, j) + 1) * min(i, j) / 2 + max(i, j) -
              min(i, j);
          outfile << tab(j) << cPar.Ve_remle_null[c];
        }
        outfile << endl;
      }
      outfile << "## se(Ve): " << endl;
      for (size_t i = 0; i < cPar.n_ph; i++) {
        for (size_t j = 0; j <= i; j++) {
          c = (2 * cPar.n_ph - min(i, j) + 1) * min(i, j) / 2 + max(i, j) -
              min(i, j);
          outfile << tab(j) << safe_sqrt(cPar.VVe_remle_null[c]);
        }
        outfile << endl;
      }

      outfile << "## MLE estimate for Vg in the null model: " << endl;
      for (size_t i = 0; i < cPar.n_ph; i++) {
        for (size_t j = 0; j < cPar.n_ph; j++) {
          c = (2 * cPar.n_ph - min(i, j) + 1) * min(i, j) / 2 + max(i, j) -
              min(i, j);
          outfile << tab(j) << cPar.Vg_mle_null[c];
        }
        outfile << endl;
      }
      outfile << "## se(Vg): " << endl;
      for (size_t i = 0; i < cPar.n_ph; i++) {
        for (size_t j = 0; j <= i; j++) {
          c = (2 * cPar.n_ph - min(i, j) + 1) * min(i, j) / 2 + max(i, j) -
              min(i, j);
          outfile << tab(j) << safe_sqrt(cPar.VVg_mle_null[c]);
        }
        outfile << endl;
      }
      outfile << "## MLE estimate for Ve in the null model: " << endl;
      for (size_t i = 0; i < cPar.n_ph; i++) {
        for (size_t j = 0; j < cPar.n_ph; j++) {
          c = (2 * cPar.n_ph - min(i, j) + 1) * min(i, j) / 2 + max(i, j) -
              min(i, j);
          outfile << tab(j) << cPar.Ve_mle_null[c];
        }
        outfile << endl;
      }
      outfile << "## se(Ve): " << endl;
      for (size_t i = 0; i < cPar.n_ph; i++) {
        for (size_t j = 0; j <= i; j++) {
          c = (2 * cPar.n_ph - min(i, j) + 1) * min(i, j) / 2 + max(i, j) -
              min(i, j);
          outfile << tab(j) << safe_sqrt(cPar.VVe_mle_null[c]);
        }
        outfile << endl;
      }
      outfile << "## estimate for B (d by c) in the null model (columns "
                 "correspond to the covariates provided in the file): "
              << endl;
      for (size_t i = 0; i < cPar.n_ph; i++) {
        for (size_t j = 0; j < cPar.n_cvt; j++) {
          c = i * cPar.n_cvt + j;
          outfile << tab(j) << cPar.beta_remle_null[c];
        }
        outfile << endl;
      }
      outfile << "## se(B): " << endl;
      for (size_t i = 0; i < cPar.n_ph; i++) {
        for (size_t j = 0; j < cPar.n_cvt; j++) {
          c = i * cPar.n_cvt + j;
          outfile << tab(j) << cPar.se_beta_remle_null[c];
        }
        outfile << endl;
      }
    }
  }

  if (cPar.a_mode == 11 || cPar.a_mode == 12 || cPar.a_mode == 13 ||
      cPar.a_mode == 14 || cPar.a_mode == 16) {
    outfile << "## estimated mean = " << cPar.pheno_mean << endl;
  }

  if (cPar.a_mode == 11 || cPar.a_mode == 13) {
    outfile << "##" << endl;
    outfile << "## MCMC related:" << endl;
    outfile << "## initial value of h = " << cPar.cHyp_initial.h << endl;
    outfile << "## initial value of rho = " << cPar.cHyp_initial.rho << endl;
    outfile << "## initial value of pi = " << exp(cPar.cHyp_initial.logp)
            << endl;
    outfile << "## initial value of |gamma| = " << cPar.cHyp_initial.n_gamma
            << endl;
    outfile << "## random seed = " << cPar.randseed << endl;
    outfile << "## acceptance ratio = "
            << (double)cPar.n_accept /
                   (double)((cPar.w_step + cPar.s_step) * cPar.n_mh)
            << endl;
  }

  outfile << "##" << endl;
  outfile << "## Computation Time:" << endl;
  outfile << "## total computation time = " << cPar.time_total << " min "
          << endl;
  outfile << "## computation time break down: " << endl;
  if (cPar.a_mode == M_KIN || cPar.a_mode == M_KIN2 || cPar.a_mode == 11 ||
      cPar.a_mode == 13 || cPar.a_mode == 14 || cPar.a_mode == 16) {
    outfile << "##      time on calculating relatedness matrix = "
            << cPar.time_G << " min " << endl;
  }
  if (cPar.a_mode == M_EIGEN) {
    outfile << "##      time on eigen-decomposition = " << cPar.time_eigen
            << " min " << endl;
  }
  if (cPar.a_mode == M_LMM1 || cPar.a_mode == M_LMM2 || cPar.a_mode == M_LMM3 ||
      cPar.a_mode == M_LMM4 || cPar.a_mode == M_LMM5 || cPar.a_mode == 11 ||
      cPar.a_mode == 12 || cPar.a_mode == 13 || cPar.a_mode == 14 ||
      cPar.a_mode == 16) {
    outfile << "##      time on eigen-decomposition = " << cPar.time_eigen
            << " min " << endl;
    outfile << "##      time on calculating UtX = " << cPar.time_UtX << " min "
            << endl;
  }
  if ((cPar.a_mode >= 1 && cPar.a_mode <= 4) ||
      (cPar.a_mode >= 51 && cPar.a_mode <= 54)) {
    outfile << "##      time on optimization = " << cPar.time_opt << " min "
            << endl;
  }
  if (cPar.a_mode == 11 || cPar.a_mode == 13) {
    outfile << "##      time on proposal = " << cPar.time_Proposal << " min "
            << endl;
    outfile << "##      time on mcmc = " << cPar.time_opt << " min " << endl;
    outfile << "##      time on Omega = " << cPar.time_Omega << " min " << endl;
  }
  if (cPar.a_mode == 41 || cPar.a_mode == 42) {
    outfile << "##      time on eigen-decomposition = " << cPar.time_eigen
            << " min " << endl;
  }
  if (cPar.a_mode == 43) {
    outfile << "##      time on eigen-decomposition = " << cPar.time_eigen
            << " min " << endl;
    outfile << "##      time on predicting phenotypes = " << cPar.time_opt
            << " min " << endl;
  }
  outfile << "##" << endl;

  outfile.close();
  outfile.clear();

  info_msg("Done");
  return;
}
