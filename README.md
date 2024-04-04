# Git Repository for Identifying Chemical Signatures of Household Products Using Suspect Screening Analysis

This repo provides the data, code (functions and a vignette), results, and figures that were a product of a research project conducted on household consumer product samples from 5 different product categories (cotton clothing, fabric upholstery, shampoo, baby soap, and silicone kitchen utensils).  The project was carried out by Dr. Zach Stanfield while working as a postdoctoral researcher at the U.S. EPA.  **_The views expressed in this repository are those of the author and do not necessarily represent the views or policies of the U.S. EPA.  Mention of trade names or commercial products do not constitute endorsement or recommendation for use._**

## Description
Consumer products are a major source of chemical exposure and therefore potential risk. Knowing what chemicals are typically present in different types of products is needed for risk evaluation. 
 Assessing similarity of new products with existing ones can identify any uncommon chemical ingredients.  As an ongoing effort in exposure science, researchers are attemping to fully characterize sources of chemical exposure, one major pathway of which is via use of consumer products in the home.  Different companies and manufacturers have different production processes and chemical ingredients in the products they produce, even if those products are highly similar (e.g. both are a cotton shirt).  Furthermore, these companies and manufacturers are not required to report the full list of chemicals and their concentrations within their products.  To fully characterize the chemical make-up of these products for the sake of risk assessment, exposure to chemical mixtures, etc., suspect screening analysis is being used.  Suspect Screening analysis or SSA is an analytical chemistry method to detect all chemical constituents in a sample and annotate those chemicals using a database of chemicals with known/previously confirmed mass spectra.  
 
For this work, a list of household consumer product categories and products of certain types were selected using multiple sources of information.  A range of products from five product categories were purchased from  various department stores, resulting in 118 unique products to be analyzed using two-dimensional gas chromatography time-of-flight mass spectrometry (GC x GC-TOFMS).  Some samples were duplicates of the same item or repeats (one sample from each of two identical product items), for a total of 170 samples, in order to asses variability of the analytical method and of the individaul products themselves.  Products were extracted with dichloromethane (DCM). After addition of an internal standard at 1 ppm, each extraction was analyzed via GC X GC-TOFMS to obtain its mass spectra. The spectra were matched to the 2017 NIST database and analytical standards were used to confirm a subset of the chemical identifications.  A total of 489 unique chemicals were confirmed or tentatively identified. Chemicals were annotated by reported or predicted functional uses and structural classification via ClassyFire. This study provides a baseline set of chemical ingredients (that is, representative mixtures) across common types of consumer products, which will help in evaluating new and existing products. Separating constituent chemicals into typical and atypical might inform exposure assessment, in vitro bioactivity screening, and ultimately the risk related to using such products.

## Repo Organization
* raw_data/
  + SSA data and files provided by the task order contractor   
* data/
  + Intermediate instances of the data during processing steps
  + All files needed to complete data processing or downstream analyses
* R/
  + Functions
* results/
  + Analysis outcomes
  + Final data tables for figures
* Rmd/
  + The single vignette file needed to complete all analyses, titled mainAnalysisPipeline.Rmd
* plots/
  + All figures grouped by analysis type
* other
  + Work in progress poster presented at the International Society of Exposure Science (ISES) 2023 conference
  + A file containing references and relevant papers for the poster presented at ISES 2023
 

