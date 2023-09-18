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
   


