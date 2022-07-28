# script for downloading and initializing the single-cell RNA-seq database for the metacell package
options(timeout = 1e9)

download_raw_umi_tables <- function() {
    download.file("https://dnmt3ab-eb.s3.eu-west-1.amazonaws.com/umi_tables.tar.gz", "umi_tables.tar.gz")

    system("tar -xzvf umi_tables.tar.gz")

    file.remove("umi_tables.tar.gz")
}


download_rna_data <- function() {
    download.file("https://dnmt3ab-eb.s3.eu-west-1.amazonaws.com/rna.tar.gz", "rna.tar.gz")
    
    system("tar -xzvf rna.tar.gz")

    file.remove("rna.tar.gz")
}

download_methylation_data <- function() {
    download.file("https://dnmt3ab-eb.s3.eu-west-1.amazonaws.com/methylation_data.tar.gz", "methylation_data.tar.gz")

    system("tar -xzvf methylation_data.tar.gz")

    file.remove("methylation_data.tar.gz")
}

download_all_data <- function() {
    download_raw_umi_tables()    
    download_rna_data()
    download_methylation_data()    
}


if (!dir.exists("rna/figs")) {
    dir.create("rna/figs")
}

if (!dir.exists("rna/paper_figs")) {
    dir.create("rna/paper_figs")
}

if (!dir.exists("methylation/output")) {
    dir.create("methylation/output")
}
