library(MCView)
create_project("Dnmt3ab_EB")

import_dataset_metacell1(
  project = "Dnmt3ab_EB",
  dataset = "WT",
  scdb = "/net/mraid14/export/tgdata/users/atanay/proj/ebdnmt/scrna_db/",
  matrix = "eseb_07_wt",
  mc = "eseb_07_wt_bs500f_got_s2m",
  mc2d = "eseb_07_wt_bs500f_got_s2m",
  metacell_types_file = "metacell_types_wt.csv",
  cell_type_colors_file = "cell_type_colors.csv",  
  time_bin_field = "EB_day"
)

import_dataset_metacell1(
  project = "Dnmt3ab_EB",
  dataset = "3A",
  scdb = "/net/mraid14/export/tgdata/users/atanay/proj/ebdnmt/scrna_db/",
  matrix = "eseb_07_3a",
  mc = "eseb_07_3a_bs500f_got_s2m",
  mc2d = "eseb_07_3a_bs500f_got_s2m",
  metacell_types_file = "metacell_types_3a.csv",
  cell_type_colors_file = "cell_type_colors.csv",  
  time_bin_field = "EB_day"
)

import_dataset_metacell1(
  project = "Dnmt3ab_EB",
  dataset = "3B",
  scdb = "/net/mraid14/export/tgdata/users/atanay/proj/ebdnmt/scrna_db/",
  matrix = "eseb_07_3b",
  mc = "eseb_07_3b_bs500f_got_s2m",
  mc2d = "eseb_07_3b_bs500f_got_s2m",
  metacell_types_file = "metacell_types_3b.csv",
  cell_type_colors_file = "cell_type_colors.csv",  
  time_bin_field = "EB_day"
)

import_dataset_metacell1(
  project = "Dnmt3ab_EB",
  dataset = "DKO",
  scdb = "/net/mraid14/export/tgdata/users/atanay/proj/ebdnmt/scrna_db/",
  matrix = "eseb_07_3ab",
  mc = "eseb_07_3ab_bs500f_got_s2m",
  mc2d = "eseb_07_3ab_bs500f_got_s2m",
  metacell_types_file = "metacell_types_3ab.csv",
  cell_type_colors_file = "cell_type_colors.csv",  
  time_bin_field = "EB_day"
)

import_dataset_metacell1(
  project = "Dnmt3ab_EB",
  dataset = "TKO",
  scdb = "/net/mraid14/export/tgdata/users/atanay/proj/ebdnmt/scrna_db/",
  matrix = "eseb_07_TKO",
  mc = "eseb_07_TKO_bs500f_got_s2m",
  mc2d = "eseb_07_TKO_bs500f_got_s2m",
  metacell_types_file = "metacell_types_TKO.csv",
  cell_type_colors_file = "cell_type_colors.csv",  
  time_bin_field = "EB_day"
)