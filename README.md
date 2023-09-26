# VaccineResponse_GWAS
**Background**: The impact of genetics on the heterogeneous immune response and risks of breakthrough infections remains unclear for the different COVID-19 vaccines. Understanding the role of genetic variants on humoral response and subsequent breakthrough outcomes is key to reveal mechanisms underlying vaccine effectiveness. 

**Objective**: To perform four genome wide association studies among vaccinated participants for COVID-19 vaccine antibody response and breakthrough outcomes. 

    
For more information: [LINK TO THE ARTICLE]

# Running the analysis
Previously, you must ensure that you have a directory "${folder directory}" with the following files:
 - breakthrough.rds [RDS file]: File with breakthrough information. Columns required include
   
        -> eid: individual identification
        -> event_dt: vaccination date
        -> case_control_index: if it was infected within 84 days after vaccination (Cases) // not infected within 84 days after vaccination (Control)
        -> breakthrough_infection_type: Omicron / No information (NA) / Pre_Omicron
        -> origin: if it was hospitalised (1) / not hospitalised (0) / no information (NA)
   
 - immune_response.tab [TAB FILE]: File with immune response information. Columns required include

        -> eid
        -> Antibody_test_result: Field 27981.0.0 from UK Biobank
        -> Date_antibody_test_performed: Field 27982.0.0 from UK Biobank
        -> Received_first_COVID19_vaccination: Field 27983.0.0 from UK Biobank
        -> Date_first_COVID19_vaccination: Field 27984.0.0 from UK Biobank
        -> Received_second_COVID19_vaccination: Field 27985.0.0 from UK Biobank
        -> Date_second_COVID19_vaccination: Field 27986.0.0 from UK Biobank
        -> Covid_antibody_test_results: Field 27990.0.0 from UK Biobank
   
 - ukb669864_covariates.csv [CSV file]: File with covariates information. Columns required include

        -> Sex: Field 31 from UK Biobank
        -> BMI: Field 21001.0.0 from UK Biobank
        -> PC1,...,PC10: Field 22009 from UK Biobank
        -> Year_of_birth: Field 34 from UK Biobank
        -> Genetic batch: Field 22000 from UK Biobank
        -> Caucasian: Field 22006 from UK Biobank
        -> sex_chromosome_aneuploidy: Field 22019 from UK Biobank
        -> kinship_to_other_participants: Field 22021 from UK Biobank
        -> genetic_sex: Field 22001 from UK Biobank

Also, the following variables must be specified:

    -> dir_data <- "...": The path to the data directory.
    -> dir_results <- "...": The path to the results directory

Once you perform the GWAS with FUMA, you must download the following files to continue with the analysis:
 - All GWAS results in a directory located in ${dir_results/GWAS/} and with the following names: breakthrough_ImputedData_covidSeverity.txt, breakthrough_ImputedData_covidSusceptibility.txt, immuneResponse_ImputedData_one_dose_cohort.txt, immuneResponse_ImputedData_two_dose_cohort.txt
 - SNPs (annotations) files in a directory located in ${dir_results/Mapping/}
    
