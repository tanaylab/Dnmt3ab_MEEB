library(MCView)
dir.create("Dnmt3ab_MEEB")

import_dataset_metacell1(
  project = "Dnmt3ab_MEEB/WT",
  dataset = "WT",
  title = "MEEB WT",
  selected_gene1 = "Dnmt3a",
  selected_gene2 = "Dnmt3b",
  scdb = "/net/mraid14/export/tgdata/users/atanay/proj/ebdnmt/scrna_db/",
  matrix = "eseb_07_wt",
  mc = "eseb_07_wt_bs500f_got_s2m",
  mc2d = "eseb_07_wt_bs500f_got_s2m",
  metacell_types_file = "metacell_types_wt.csv",
  cell_type_colors_file = "cell_type_colors.csv",  
  time_bin_field = "EB_day",
  metadata_fields = c("EB_day", "exp_state_r"),
  categorical = c("exp_state_r")
)

create_bundle("Dnmt3ab_MEEB/WT", name = "WT", path = "/home/aviezerl/tanaywiz/apps/zohar_mukamel/Dnmt3ab_MEEB", overwrite = TRUE)

import_dataset_metacell1(
  project = "Dnmt3ab_MEEB/3A",
  dataset = "3A",
  title = "MEEB Dnmt3a-/-",
  selected_gene1 = "Dnmt3a",
  selected_gene2 = "Dnmt3b",
  scdb = "/net/mraid14/export/tgdata/users/atanay/proj/ebdnmt/scrna_db/",
  matrix = "eseb_07_3a",
  mc = "eseb_07_3a_bs500f_got_s2m",
  mc2d = "eseb_07_3a_bs500f_got_s2m",
  metacell_types_file = "metacell_types_3a.csv",
  cell_type_colors_file = "cell_type_colors.csv",  
  time_bin_field = "EB_day",
  metadata_fields = c("EB_day", "exp_state_r"),
  categorical = c("exp_state_r")
)

create_bundle("Dnmt3ab_MEEB/3A", name = "3A", path = "/home/aviezerl/tanaywiz/apps/zohar_mukamel/Dnmt3ab_MEEB", overwrite = TRUE)

import_dataset_metacell1(
  project = "Dnmt3ab_MEEB/3B",
  dataset = "3B",
  title = "MEEB Dnmt3b-/-",
  selected_gene1 = "Dnmt3a",
  selected_gene2 = "Dnmt3b",
  scdb = "/net/mraid14/export/tgdata/users/atanay/proj/ebdnmt/scrna_db/",
  matrix = "eseb_07_3b",
  mc = "eseb_07_3b_bs500f_got_s2m",
  mc2d = "eseb_07_3b_bs500f_got_s2m",
  metacell_types_file = "metacell_types_3b.csv",
  cell_type_colors_file = "cell_type_colors.csv",  
  time_bin_field = "EB_day",
  metadata_fields = c("EB_day", "exp_state_r"),
  categorical = c("exp_state_r")
)

create_bundle("Dnmt3ab_MEEB/3B", name = "3B", path = "/home/aviezerl/tanaywiz/apps/zohar_mukamel/Dnmt3ab_MEEB", overwrite = TRUE)

import_dataset_metacell1(
  project = "Dnmt3ab_MEEB/DKO",
  dataset = "DKO",
  title = "MEEB DKO",
  selected_gene1 = "Dnmt3a",
  selected_gene2 = "Dnmt3b",
  scdb = "/net/mraid14/export/tgdata/users/atanay/proj/ebdnmt/scrna_db/",
  matrix = "eseb_07_3ab",
  mc = "eseb_07_3ab_bs500f_got_s2m",
  mc2d = "eseb_07_3ab_bs500f_got_s2m",
  metacell_types_file = "metacell_types_3ab.csv",
  cell_type_colors_file = "cell_type_colors.csv",  
  time_bin_field = "EB_day",
  metadata_fields = c("EB_day", "exp_state_r"),
  categorical = c("exp_state_r")
)

create_bundle("Dnmt3ab_MEEB/DKO", name = "DKO", path = "/home/aviezerl/tanaywiz/apps/zohar_mukamel/Dnmt3ab_MEEB", overwrite = TRUE)

import_dataset_metacell1(
  project = "Dnmt3ab_MEEB/TKO",
  dataset = "TKO",
  title = "MEEB TKO",
  selected_gene1 = "Dnmt3a",
  selected_gene2 = "Dnmt3b",
  scdb = "/net/mraid14/export/tgdata/users/atanay/proj/ebdnmt/scrna_db/",
  matrix = "eseb_07_TKO",
  mc = "eseb_07_TKO_bs500f_got_s2m",
  mc2d = "eseb_07_TKO_bs500f_got_s2m",
  metacell_types_file = "metacell_types_TKO.csv",
  cell_type_colors_file = "cell_type_colors.csv",  
  time_bin_field = "EB_day",
  metadata_fields = c("EB_day", "exp_state_r"),
  categorical = c("exp_state_r")
)

create_bundle("Dnmt3ab_MEEB/TKO", name = "TKO", path = "/home/aviezerl/tanaywiz/apps/zohar_mukamel/Dnmt3ab_MEEB", overwrite = TRUE)

file.copy("index.html", "/home/aviezerl/tanaywiz/apps/zohar_mukamel/Dnmt3ab_MEEB")

