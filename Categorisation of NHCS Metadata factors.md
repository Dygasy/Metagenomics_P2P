Extensive list of metadata factors- in this case we have anthropometric, clinical, biochemical, echocardiographic, lifestyle, pyschological, metabolic, dietary, and genetic data all combined. 

anthropometric
* BMI, Weight, Height, WaistCircumference, BSA, WHR
* Lean mass, Fat mass (SMM, BFM, PBF, ALM, Lean_LA, Lean_RA, Lean_LL, Lean_RL, Lean_T)
* BMR, Fitness_score

Blood Pressre and Cardiovascular hemodynamics
* SBP, DBP, Pulse, MAP, PP, RPP, PWVms, AugmentationIndex, PeripheralResistance

Echocardiographic Parameters:
* LV dimensions: IVSD, IVSS, LVIDD, LVIDS, LVPWD, LVPWS
*LV function: LVEF, LVFS, LVmass, LVmass_index
*LA volumes, MV E/A, TR Vmax, PASP, RAP
*Aortic root diameters: sinus, AO, LVOT
*Diastolic function: MV E/A ratio, IVRT, MV_DT
*RV function: RVTA
*Regurgitation & valvular measures: AR, MR, TR, PR

Respiratory/ Exercise
* VO2Max, O2 delivery indices

Metabolic/Lipids/ Amino acids
*Detailed lipid species: C2, C3, C4, ... C28, and complex derivatives (e.g. C203OHC183DC)

Amino acids: Gly, Ala, Ser, Pro, Val, Leu, Ile, Orn, Met, His, Phe, Arg, Cit, Tyr, Asp, Glu, Trp

Lifestyle & Diet
* PhysicalActivity: frequency, intensity, duration
* Diet: meat, pork, beef, chicken, fruit/vegetable, dairy, overall diet scores
* Sleep: SleepQ1-7, recoded scores
* Self-empowerment, motivation
  
Psychological & Quality of Life Scores
* SF-36 like scales: PF, RP, BP, GH, VT, SF, RE, MH
* MentalHealthScore, PhysicalHealthScore
* Pain interference, Social activities, Emotional domains
* NBS (norm-based) scores, SF6D utility

Additonal Biochemistry/Urine
* Creatinine, albumin urine, microalbumin ratio

Composite Cardiovascular Risk
* Risk, Risk_exact, VascularAge, GENPOP_PCS, GENPOP_MCS

We can use this groups for dimensional reduction (PCA/ clustering) or adjusting covariates
create composite scores (eg: echocardiographic score, metabolic score)
Use lifestyle/diet/psychological scores to stratify groups for outcome analysis





