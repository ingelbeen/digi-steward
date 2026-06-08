# digi-steward
Sample size estimation for a 12-month stepped wedge study. There are two substudies:
1. The stepped wedge study to estimate the AI-generated videos effect on correct AMU/treatment adherence, evaluated at D7 after consultation (primary outcome)
   Assumptions:
- Per study site, at each provider monthly 20 patients surveys (~1 day)
- Of 20 surveyed patients, 10 (50%) were dispensed an antibiotic and of those, 4 (40%) for a correct indication. Based on:
      Sulis et al Plos Medicine 2020
      Valia et al JAC 2024
      Ingelbeen et al Lancet Infect Dis 2023
      Ardillon et al Plos Medicine 2023
- 2 provider(s) per month are added to the intervention group
- 50% correct antibiotic intake at baseline, with an expected increase in correct antibiotic use to 65% due to the intervention (ARR=15%)
- intraclass correlation coefficient (ICC) of 0.05 based on antibiotic use captured in exit surveys at different community-level providers in Ingelbeen, Valia et al Lancet Infect Dis 2026

2. A longitudinal cohort of the general population to evaluate the indirect intervention effect on healthcare utilisation and the per capita rate of antibiotic use (exploratory outcome; limited to HDSS in Nanoro, Burkina Faso, and Agincourt, South Africa)
   Assumptions
- at baseline 2.0 outpatient visits per capita per year (Burkina Faso), 3.58 outpatient visits per capita per year (South Africa). Based on 2016 estimates from appendix Moses et al Lancet Public Health 2018
- at baseline 73% of outpatients visits to primary care clinics, 23% to pharmacies or other medicine vendors. Based on Table S1 in appendix Ingelbeen, Valia et al Lancet Infect Dis 2026
- 50% AMU prevalence suring primary care clinic visits, based on Sulis et al Plos MEdicine 2020; 25% AMU during pharmacy visits, based on Valia et al JAC 2024
- simple random sample (random selection of community members from the HDSS database)
