# Author: Thibaud Derippe
# Set up file containing the main Lindner Model, some initial parameter and 
# some additional configuration information for using PaSM for creating VCs.

# Part 1: write the model equations ---------------------------------------

#### Write the equations of your model, followng RxODE format


model_RxODE <-RxODE({
  
  
  Bcl2_I(0) <- Bcl2_I0
  Bclxl_I(0) <- Bclxl_I0
  Mcl1_I(0) <- Mcl1_I0
  VDAC2(0) <- BAK0
  BAK(0) <- BAK0
  BAXc(0) <- BAXc0
  NOXA(0) <- NOXA0
  PUMA(0) <- PUMA0
  BIM(0) <- BIM0
  Mcl1(0) <- Mcl10
  Bclxl(0) <- Bclxl0
  Bcl2(0) <- Bcl20
  
  
  
  kdeg_Bcl2 <- 0.139
  kprod_Bcl2 <- Bcl20 * kdeg_Bcl2
  kdeg_Bclxl <- 0.139
  kprod_Bclxl <- Bclxl0 * kdeg_Bclxl
  kdeg_Mcl1 <- 0.00068 * 3600
  kprod_Mcl1 <- Mcl10 * kdeg_Mcl1
  kdeg_BIM <- 0.173
  kdeg_tBID <- 0.554
  kdeg_PUMA <- 0.204
  kdeg_NOXA <- 0.695
  kdeg_Bcl2_BIM <- 0.554
  kforward_Bcl2_BIM <- 0.108
  kbackward_Bcl2_BIM <- 0.504
  kdeg_Bclxl_BIM <- 0.554
  kforward_Bclxl_BIM <- 1.98
  kbackward_Bclxl_BIM <- 1.58
  kdeg_Mcl1_BIM <- 0.277
  kforward_Mcl1_BIM <- 4.68
  kbackward_Mcl1_BIM <- 0.936
  kdeg_Bcl2_tBID <- 0.554
  kforward_Bcl2_tBID <- 0.126
  kbackward_Bcl2_tBID <- 0.504
  kdeg_Bclxl_tBID <- 0.554
  kforward_Bclxl_tBID <- 0.0583
  kbackward_Bclxl_tBID <- 1.58
  kdeg_Mcl1_tBID <- 0.554
  kforward_Mcl1_tBID <- 0.0947
  kbackward_Mcl1_tBID <- 0.936
  kdeg_Bcl2_PUMA <- 0.554
  kforward_Bcl2_PUMA <- 0.028
  kbackward_Bcl2_PUMA <- 0.504
  kdeg_Bclxl_PUMA <- 0.554
  kforward_Bclxl_PUMA <- 0.311
  kbackward_Bclxl_PUMA <- 1.58
  kdeg_Mcl1_PUMA <- 0.277
  kforward_Mcl1_PUMA <- 0.493
  kbackward_Mcl1_PUMA <- 0.936
  kdeg_Bcl2_NOXA <- 0.554
  kforward_Bcl2_NOXA <- 0.00262
  kbackward_Bcl2_NOXA <- 0.504
  kdeg_Bclxl_NOXA <- 0.554
  kforward_Bclxl_NOXA <- 0.000158
  kbackward_Bclxl_NOXA <- 1.58
  kdeg_Mcl1_NOXA <- 0.925
  kforward_Mcl1_NOXA <- 0.0237
  kbackward_Mcl1_NOXA <- 0.936
  kforward_Bcl2_BAKa <- 5.04e-05
  kbackward_Bcl2_BAKa <- 0.504
  kforward_Bcl2_BAXma <- 0.0336
  kbackward_Bcl2_BAXma <- 0.504
  kforward_Bclxl_BAKa <- 0.0198
  kbackward_Bclxl_BAKa <- 1.58
  kforward_Bclxl_BAXma <- 0.00186
  kbackward_Bclxl_BAXma <- 1.58
  kforward_Mcl1_BAKa <- 0.117
  kbackward_Mcl1_BAKa <- 0.936
  kforward_Mcl1_BAXma <- 9.36e-05
  kbackward_Mcl1_BAXma <- 0.936
  kforward_BIM_BAXc <- 0.00925
  kbackward_BIM_BAXc <- 0.925
  kforward_tBID_BAXc <- 0.00925
  kbackward_tBID_BAXc <- 0.925
  kforward_PUMA_BAXc <- 0.00925
  kbackward_PUMA_BAXc <- 0.925
  kforward_BIM_BAK <- 0.00925
  kbackward_BIM_BAK <- 0.925
  kforward_tBID_BAK <- 0.00925
  kbackward_tBID_BAK <- 0.925
  kforward_PUMA_BAK <- 0.00925
  kbackward_PUMA_BAK <- 0.925
  k_BIM_BAK <- 418
  k_tBID_BAK <- 418
  k_PUMA_BAK <- 418
  k_BIM_BAXc <- 418
  k_tBID_BAXc <- 418
  k_PUMA_BAXc <- 418
  k_BAXca <- 418
  kforward_BAK_VDAC2 <- 0.00832
  kbackward_BAK_VDAC2 <- 8.32
  kforward_BAK2 <- 0.0461
  kbackward_BAK2 <- 0.695
  kforward_BAK4 <- 0.0461
  kbackward_BAK4 <- 0.695
  kforward_BAK6 <- 0.0461
  kbackward_BAK6 <- 0.695
  kforward_BAK8 <- 0.0461
  kbackward_BAK8 <- 0.695
  kforward_BAK10 <- 0.0461
  kbackward_BAK10 <- 0.695
  kforward_BAK12 <- 0.0461
  kbackward_BAK12 <- 0.695
  kforward_BAK8 <- 0.0461
  kbackward_BAK8 <- 0.695
  kforward_BAK10 <- 0.0461
  kbackward_BAK10 <- 0.695
  kforward_BAK12 <- 0.0461
  kbackward_BAK12 <- 0.695
  kforward_BAK12 <- 0.0461
  kbackward_BAK12 <- 0.695
  kforward_BAX2 <- 0.0461
  kbackward_BAX2 <- 0.695
  kforward_BAX4 <- 0.0461
  kbackward_BAX4 <- 0.695
  kforward_BAX6 <- 0.0461
  kbackward_BAX6 <- 0.695
  kforward_BAX8 <- 0.0461
  kbackward_BAX8 <- 0.695
  kforward_BAX10 <- 0.0461
  kbackward_BAX10 <- 0.695
  kforward_BAX12 <- 0.0461
  kbackward_BAX12 <- 0.695
  kforward_BAX8 <- 0.0461
  kbackward_BAX8 <- 0.695
  kforward_BAX10 <- 0.0461
  kbackward_BAX10 <- 0.695
  kforward_BAX12 <- 0.0461
  kbackward_BAX12 <- 0.695
  kforward_BAX12 <- 0.0461
  kbackward_BAX12 <- 0.695
  R1_deg_Bcl2 <- -kdeg_Bcl2 * Bcl2
  R2_prod_Bcl2 <- kprod_Bcl2
  R3_deg_Bclxl <- -kdeg_Bclxl * Bclxl
  R4_prod_Bclxl <- kprod_Bclxl
  R5_deg_Mcl1 <- -kdeg_Mcl1 * Mcl1
  R6_prod_Mcl1 <- kprod_Mcl1 *  (1 - A15 ** hillA15 / (A15  ** hillA15 + EC50A15 ** hillA15))
  R7_deg_BIM <- -kdeg_BIM * BIM
  R8_deg_tBID <- -kdeg_tBID * tBID
  R9_deg_PUMA <- -kdeg_PUMA * PUMA
  R10_deg_NOXA <- -kdeg_NOXA * NOXA
  R11_deg_Bcl2_BIM <- -kdeg_Bcl2_BIM * Bcl2_BIM
  R12_complex_Bcl2_BIM <- kforward_Bcl2_BIM * Bcl2 * BIM - kbackward_Bcl2_BIM * Bcl2_BIM
  R13_deg_Bclxl_BIM <- -kdeg_Bclxl_BIM * Bclxl_BIM
  R14_complex_Bclxl_BIM <- kforward_Bclxl_BIM * Bclxl * BIM - kbackward_Bclxl_BIM * Bclxl_BIM
  R15_deg_Mcl1_BIM <- -kdeg_Mcl1_BIM * Mcl1_BIM
  R16_complex_Mcl1_BIM <- kforward_Mcl1_BIM * Mcl1 * BIM - kbackward_Mcl1_BIM * Mcl1_BIM
  R17_deg_Bcl2_tBID <- -kdeg_Bcl2_tBID * Bcl2_tBID
  R18_complex_Bcl2_tBID <- kforward_Bcl2_tBID * Bcl2 * tBID - kbackward_Bcl2_tBID * Bcl2_tBID
  R19_deg_Bclxl_tBID <- -kdeg_Bclxl_tBID * Bclxl_tBID
  R20_complex_Bclxl_tBID <- kforward_Bclxl_tBID * Bclxl * tBID - kbackward_Bclxl_tBID * Bclxl_tBID
  R21_deg_Mcl1_tBID <- -kdeg_Mcl1_tBID * Mcl1_tBID
  R22_complex_Mcl1_tBID <- kforward_Mcl1_tBID * Mcl1 * tBID - kbackward_Mcl1_tBID * Mcl1_tBID
  R23_deg_Bcl2_PUMA <- -kdeg_Bcl2_PUMA * Bcl2_PUMA
  R24_complex_Bcl2_PUMA <- kforward_Bcl2_PUMA * Bcl2 * PUMA - kbackward_Bcl2_PUMA * Bcl2_PUMA
  R25_deg_Bclxl_PUMA <- -kdeg_Bclxl_PUMA * Bclxl_PUMA
  R26_complex_Bclxl_PUMA <- kforward_Bclxl_PUMA * Bclxl * PUMA - kbackward_Bclxl_PUMA * Bclxl_PUMA
  R27_deg_Mcl1_PUMA <- -kdeg_Mcl1_PUMA * Mcl1_PUMA
  R28_complex_Mcl1_PUMA <- kforward_Mcl1_PUMA * Mcl1 * PUMA - kbackward_Mcl1_PUMA * Mcl1_PUMA
  R29_deg_Bcl2_NOXA <- -kdeg_Bcl2_NOXA * Bcl2_NOXA
  R30_complex_Bcl2_NOXA <- kforward_Bcl2_NOXA * Bcl2 * NOXA - kbackward_Bcl2_NOXA * Bcl2_NOXA
  R31_deg_Bclxl_NOXA <- -kdeg_Bclxl_NOXA * Bclxl_NOXA
  R32_complex_Bclxl_NOXA <- kforward_Bclxl_NOXA * Bclxl * NOXA - kbackward_Bclxl_NOXA * Bclxl_NOXA
  R33_deg_Mcl1_NOXA <- -kdeg_Mcl1_NOXA * Mcl1_NOXA
  R34_complex_Mcl1_NOXA <- kforward_Mcl1_NOXA * Mcl1 * NOXA - kbackward_Mcl1_NOXA * Mcl1_NOXA
  R35_complex_Bcl2_BAKa <- kforward_Bcl2_BAKa * Bcl2 * BAKa - kbackward_Bcl2_BAKa * Bcl2_BAKa
  R36_complex_Bcl2_BAXma <- kforward_Bcl2_BAXma * Bcl2 * BAXma - kbackward_Bcl2_BAXma * Bcl2_BAXma
  R37_complex_Bclxl_BAKa <- kforward_Bclxl_BAKa * Bclxl * BAKa - kbackward_Bclxl_BAKa * Bclxl_BAKa
  R38_complex_Bclxl_BAXma <- kforward_Bclxl_BAXma * Bclxl * BAXma - kbackward_Bclxl_BAXma * Bclxl_BAXma
  R39_complex_Mcl1_BAKa <- kforward_Mcl1_BAKa * Mcl1 * BAKa - kbackward_Mcl1_BAKa * Mcl1_BAKa
  R40_complex_Mcl1_BAXma <- kforward_Mcl1_BAXma * Mcl1 * BAXma - kbackward_Mcl1_BAXma * Mcl1_BAXma
  R41_complex_BIM_BAXc <- kforward_BIM_BAXc * BIM * BAXc - kbackward_BIM_BAXc * BIM_BAXc
  R42_complex_tBID_BAXc <- kforward_tBID_BAXc * tBID * BAXc - kbackward_tBID_BAXc * tBID_BAXc
  R43_complex_PUMA_BAXc <- kforward_PUMA_BAXc * PUMA * BAXc - kbackward_PUMA_BAXc * PUMA_BAXc
  R44_complex_BIM_BAK <- kforward_BIM_BAK * BIM * BAK - kbackward_BIM_BAK * BIM_BAK
  R45_complex_tBID_BAK <- kforward_tBID_BAK * tBID * BAK - kbackward_tBID_BAK * tBID_BAK
  R46_complex_PUMA_BAK <- kforward_PUMA_BAK * PUMA * BAK - kbackward_PUMA_BAK * PUMA_BAK
  R47_disso_BIM_BAK <- -k_BIM_BAK * BIM_BAK
  R48_disso_tBID_BAK <- -k_tBID_BAK * tBID_BAK
  R49_disso_PUMA_BAK <- -k_PUMA_BAK * PUMA_BAK
  R50_disso_BIM_BAXc <- -k_BIM_BAXc * BIM_BAXc
  R51_disso_tBID_BAXc <- -k_tBID_BAXc * tBID_BAXc
  R52_disso_PUMA_BAXc <- -k_PUMA_BAXc * PUMA_BAXc
  R53_disso_BAXca <- -k_BAXca * BAXca
  R54_complex_BAK_VDAC2 <- kforward_BAK_VDAC2 * BAK * VDAC2 - kbackward_BAK_VDAC2 * BAK_VDAC2
  R55_complex_BAK2 <- kforward_BAK2 * BAKa * BAKa - kbackward_BAK2 * BAK2
  R56_complex_BAK4 <- kforward_BAK4 * BAK2 * BAK2 - kbackward_BAK4 * BAK4
  R57_complex_BAK6 <- kforward_BAK6 * BAK2 * BAK4 - kbackward_BAK6 * BAK6
  R58_complex_BAK8 <- kforward_BAK8 * BAK2 * BAK6 - kbackward_BAK8 * BAK8
  R59_complex_BAK10 <- kforward_BAK10 * BAK2 * BAK8 - kbackward_BAK10 * BAK10
  R60_complex_BAK12 <- kforward_BAK12 * BAK2 * BAK10 - kbackward_BAK12 * BAK12
  R61_complex_BAK8 <- kforward_BAK8 * BAK4 * BAK4 - kbackward_BAK8 * BAK8
  R62_complex_BAK10 <- kforward_BAK10 * BAK4 * BAK6 - kbackward_BAK10 * BAK10
  R63_complex_BAK12 <- kforward_BAK12 * BAK4 * BAK8 - kbackward_BAK12 * BAK12
  R64_complex_BAK12 <- kforward_BAK12 * BAK6 * BAK6 - kbackward_BAK12 * BAK12
  R65_complex_BAX2 <- kforward_BAX2 * BAXma * BAXma - kbackward_BAX2 * BAX2
  R66_complex_BAX4 <- kforward_BAX4 * BAX2 * BAX2 - kbackward_BAX4 * BAX4
  R67_complex_BAX6 <- kforward_BAX6 * BAX2 * BAX4 - kbackward_BAX6 * BAX6
  R68_complex_BAX8 <- kforward_BAX8 * BAX2 * BAX6 - kbackward_BAX8 * BAX8
  R69_complex_BAX10 <- kforward_BAX10 * BAX2 * BAX8 - kbackward_BAX10 * BAX10
  R70_complex_BAX12 <- kforward_BAX12 * BAX2 * BAX10 - kbackward_BAX12 * BAX12
  R71_complex_BAX8 <- kforward_BAX8 * BAX4 * BAX4 - kbackward_BAX8 * BAX8
  R72_complex_BAX10 <- kforward_BAX10 * BAX4 * BAX6 - kbackward_BAX10 * BAX10
  R73_complex_BAX12 <- kforward_BAX12 * BAX4 * BAX8 - kbackward_BAX12 * BAX12
  R74_complex_BAX12 <- kforward_BAX12 * BAX6 * BAX6 - kbackward_BAX12 * BAX12
  
  
  
  
  # Drugs starts only after 700 hours
  
  if( t >= 700){
    
    drug_act <- 1
    
  }else{
    
    drug_act <- 0
  }
  
  ka_Veneto <- 0.856
  Cl_Veneto <- 0.449
  Vd_Veneto <- 3.54
  Ke_Veneto <- Cl_Veneto / Vd_Veneto
  d/dt(Veneto_gut) <- - ka_Veneto * Veneto_gut * drug_act
  d/dt(Veneto_central) <- (ka_Veneto * Veneto_gut - Ke_Veneto * Veneto_central ) * drug_act
  Veneto_plasma <-  Veneto_central / Vd_Veneto # in ng/L ; MM 868,44 g/mol
  Veneto_plasma_microMolai <- Veneto_plasma * 1000 / 868.44
  Veneto_tumor <- Veneto_plasma * ratioTumor  
  
  d/dt(A15) <- - A15 * ke_Mcl1_I
 
  
  
  d/dt(Bcl2_I) <- -ke_BCl2_I * Bcl2_I  * drug_act
  d/dt(Bclxl_I) <- -ke_BClxl_I * Bclxl_I  * drug_act
  d/dt(Mcl1_I) <- -ke_Mcl1_I * Mcl1_I  * drug_act
  d/dt(Bcl2) <- - Veneto_tumor * k2_Bcl2_I * Bcl2  * drug_act+  R1_deg_Bcl2 + R2_prod_Bcl2 - R12_complex_Bcl2_BIM - R18_complex_Bcl2_tBID - R24_complex_Bcl2_PUMA - R30_complex_Bcl2_NOXA - R35_complex_Bcl2_BAKa - R36_complex_Bcl2_BAXma - Bcl2_I * k2_Bcl2_I * Bcl2  * drug_act
  d/dt(Bclxl) <- R3_deg_Bclxl + R4_prod_Bclxl - R14_complex_Bclxl_BIM - R20_complex_Bclxl_tBID - R26_complex_Bclxl_PUMA - R32_complex_Bclxl_NOXA - R37_complex_Bclxl_BAKa - R38_complex_Bclxl_BAXma - Bclxl_I * k2_Bclxl_I * Bclxl  * drug_act
  d/dt(Mcl1) <- R5_deg_Mcl1 + R6_prod_Mcl1 - R16_complex_Mcl1_BIM - R22_complex_Mcl1_tBID - R28_complex_Mcl1_PUMA - R34_complex_Mcl1_NOXA - R39_complex_Mcl1_BAKa - R40_complex_Mcl1_BAXma - Mcl1_I * k2_Mcl1_I * Mcl1  * drug_act
  d/dt(BIM) <- R7_deg_BIM - R12_complex_Bcl2_BIM - R14_complex_Bclxl_BIM - R16_complex_Mcl1_BIM - R41_complex_BIM_BAXc - R44_complex_BIM_BAK - R47_disso_BIM_BAK - R50_disso_BIM_BAXc + BIM0 * kdeg_BIM * switchE
  d/dt(tBID) <- R8_deg_tBID - R18_complex_Bcl2_tBID - R20_complex_Bclxl_tBID - R22_complex_Mcl1_tBID - R42_complex_tBID_BAXc - R45_complex_tBID_BAK - R48_disso_tBID_BAK - R51_disso_tBID_BAXc
  d/dt(PUMA) <- R9_deg_PUMA - R24_complex_Bcl2_PUMA - R26_complex_Bclxl_PUMA - R28_complex_Mcl1_PUMA - R43_complex_PUMA_BAXc - R46_complex_PUMA_BAK - R49_disso_PUMA_BAK - R52_disso_PUMA_BAXc + PUMA0 * kdeg_PUMA * switchE
  d/dt(NOXA) <- R10_deg_NOXA - R30_complex_Bcl2_NOXA - R32_complex_Bclxl_NOXA - R34_complex_Mcl1_NOXA + NOXA0 * kdeg_NOXA * switchE
  d/dt(Bcl2_BIM) <- R11_deg_Bcl2_BIM + R12_complex_Bcl2_BIM
  d/dt(Bclxl_BIM) <- R13_deg_Bclxl_BIM + R14_complex_Bclxl_BIM
  d/dt(Mcl1_BIM) <- R15_deg_Mcl1_BIM + R16_complex_Mcl1_BIM
  d/dt(Bcl2_tBID) <- R17_deg_Bcl2_tBID + R18_complex_Bcl2_tBID
  d/dt(Bclxl_tBID) <- R19_deg_Bclxl_tBID + R20_complex_Bclxl_tBID
  d/dt(Mcl1_tBID) <- R21_deg_Mcl1_tBID + R22_complex_Mcl1_tBID
  d/dt(Bcl2_PUMA) <- R23_deg_Bcl2_PUMA + R24_complex_Bcl2_PUMA
  d/dt(Bclxl_PUMA) <- R25_deg_Bclxl_PUMA + R26_complex_Bclxl_PUMA
  d/dt(Mcl1_PUMA) <- R27_deg_Mcl1_PUMA + R28_complex_Mcl1_PUMA
  d/dt(Bcl2_NOXA) <- R29_deg_Bcl2_NOXA + R30_complex_Bcl2_NOXA
  d/dt(Bclxl_NOXA) <- R31_deg_Bclxl_NOXA + R32_complex_Bclxl_NOXA
  d/dt(Mcl1_NOXA) <- R33_deg_Mcl1_NOXA + R34_complex_Mcl1_NOXA
  d/dt(Bcl2_BAKa) <- R35_complex_Bcl2_BAKa - Bcl2_BAKa * kelimBAXBAX * switchE
  d/dt(BAKa) <- -R35_complex_Bcl2_BAKa - R37_complex_Bclxl_BAKa - R39_complex_Mcl1_BAKa - R47_disso_BIM_BAK - R48_disso_tBID_BAK - R49_disso_PUMA_BAK - R55_complex_BAK2 - R55_complex_BAK2 - BAKa * kelimBAXBAX * switchE
  d/dt(Bcl2_BAXma) <- R36_complex_Bcl2_BAXma - Bcl2_BAXma * kelimBAXBAX * switchE
  d/dt(BAXma) <- -R36_complex_Bcl2_BAXma - R38_complex_Bclxl_BAXma - R40_complex_Mcl1_BAXma - R53_disso_BAXca - R65_complex_BAX2 - R65_complex_BAX2 - BAXma * kelimBAXBAX * switchE
  d/dt(Bclxl_BAKa) <- R37_complex_Bclxl_BAKa - R36_complex_Bcl2_BAXma - R38_complex_Bclxl_BAXma - R40_complex_Mcl1_BAXma - R53_disso_BAXca - R65_complex_BAX2 - R65_complex_BAX2 - Bclxl_BAKa * kelimBAXBAX * switchE
  d/dt(Bclxl_BAXma) <- R38_complex_Bclxl_BAXma - Bclxl_BAXma * kelimBAXBAX * switchE
  d/dt(Mcl1_BAKa) <- R39_complex_Mcl1_BAKa - Mcl1_BAKa * kelimBAXBAX * switchE
  d/dt(Mcl1_BAXma) <- R40_complex_Mcl1_BAXma - Mcl1_BAXma * kelimBAXBAX * switchE
  d/dt(BIM_BAXc) <- R41_complex_BIM_BAXc + R50_disso_BIM_BAXc - BIM_BAXc * kelimBAXBAX * switchE
  d/dt(BAXc) <- -R41_complex_BIM_BAXc - R42_complex_tBID_BAXc - R43_complex_PUMA_BAXc - BAXc * kelimBAXBAX * switchE + BAXc0 * kelimBAXBAX * switchE
  d/dt(tBID_BAXc) <- R42_complex_tBID_BAXc + R51_disso_tBID_BAXc - tBID_BAXc * kelimBAXBAX * switchE
  d/dt(PUMA_BAXc) <- R43_complex_PUMA_BAXc + R52_disso_PUMA_BAXc - PUMA_BAXc * kelimBAXBAX * switchE
  d/dt(BIM_BAK) <- R44_complex_BIM_BAK + R47_disso_BIM_BAK - BIM_BAK * kelimBAXBAX * switchE
  d/dt(BAK) <- -R44_complex_BIM_BAK - R45_complex_tBID_BAK - R46_complex_PUMA_BAK - R54_complex_BAK_VDAC2 - BAK * kelimBAXBAX * switchE + BAK0 * kelimBAXBAX * switchE
  d/dt(tBID_BAK) <- R45_complex_tBID_BAK + R48_disso_tBID_BAK - tBID_BAK * kelimBAXBAX * switchE
  d/dt(PUMA_BAK) <- R46_complex_PUMA_BAK + R49_disso_PUMA_BAK - PUMA_BAK * kelimBAXBAX * switchE
  d/dt(BAXca) <- -R50_disso_BIM_BAXc - R51_disso_tBID_BAXc - R52_disso_PUMA_BAXc + R53_disso_BAXca - BAXca * kelimBAXBAX * switchE
  d/dt(BAK_VDAC2) <- R54_complex_BAK_VDAC2 - BAK_VDAC2 * kelimBAXBAX * switchE
  d/dt(VDAC2) <- -R54_complex_BAK_VDAC2 - VDAC2 * kelimBAXBAX * switchE + BAK0 * kelimBAXBAX * switchE
  d/dt(BAK2) <- R55_complex_BAK2 - R56_complex_BAK4 - R56_complex_BAK4 - R57_complex_BAK6 - R58_complex_BAK8 - R59_complex_BAK10 - R60_complex_BAK12 - BAK2 * kelimBAXBAX * switchE
  d/dt(BAK4) <- R56_complex_BAK4 - R57_complex_BAK6 - R61_complex_BAK8 - R61_complex_BAK8 - R62_complex_BAK10 - R63_complex_BAK12 - BAK4 * kelimBAXBAX * switchE
  d/dt(BAK6) <- R57_complex_BAK6 - R58_complex_BAK8 - R62_complex_BAK10 - R64_complex_BAK12 - R64_complex_BAK12 - BAK6 * kelimBAXBAX * switchE
  d/dt(BAK8) <- R58_complex_BAK8 - R59_complex_BAK10 + R61_complex_BAK8 - R63_complex_BAK12 - BAK8 * kelimBAXBAX * switchE
  d/dt(BAK10) <- R59_complex_BAK10 - R60_complex_BAK12 + R62_complex_BAK10 - BAK10 * kelimBAXBAX * switchE
  d/dt(BAK12) <- R60_complex_BAK12 + R63_complex_BAK12 + R64_complex_BAK12 - BAK12 * kelimBAXBAX * switchE
  d/dt(BAX2) <- R65_complex_BAX2 - R66_complex_BAX4 - R66_complex_BAX4 - R67_complex_BAX6 - R68_complex_BAX8 - R69_complex_BAX10 - R70_complex_BAX12 - BAX2 * kelimBAXBAX * switchE
  d/dt(BAX4) <- R66_complex_BAX4 - R67_complex_BAX6 - R71_complex_BAX8 - R71_complex_BAX8 - R72_complex_BAX10 - R73_complex_BAX12 - BAX4 * kelimBAXBAX * switchE
  d/dt(BAX6) <- R67_complex_BAX6 - R68_complex_BAX8 - R72_complex_BAX10 - R74_complex_BAX12 - R74_complex_BAX12 - BAX6 * kelimBAXBAX * switchE
  d/dt(BAX8) <- R68_complex_BAX8 - R69_complex_BAX10 + R71_complex_BAX8 - R73_complex_BAX12 - BAX8 * kelimBAXBAX * switchE
  d/dt(BAX10) <- R69_complex_BAX10 - R70_complex_BAX12 + R72_complex_BAX10 - BAX10 * kelimBAXBAX * switchE
  d/dt(BAX12) <- R70_complex_BAX12 + R73_complex_BAX12 + R74_complex_BAX12 - BAX12 * kelimBAXBAX * switchE
  pctBAX <- 100 * (6 * BAX6 + 8 * BAX8 + 10 * BAX10 + 12 * BAX12)/BAXc0
  pctBAK <- 100 * (6 * BAK6 + 8 * BAK8 + 10 * BAK10 + 12 * BAK12)/BAK0
  Pore <- 100 * (6 * BAK6 + 8 * BAK8 + 10 * BAK10 + 12 * BAK12 + 6 * BAX6 + 8 * BAX8 + 10 * BAX10 + 12 * BAX12)/(BAK0 + BAXc0)

  
  ## Add death signal
  
  if(Pore > 10){
    
    TimeAboveRate <- 1
  }else{
    
    TimeAboveRate <- 0
    
  }
  
  d/dt(TimeAbove) <- TimeAboveRate
  
})



## Verification: perform verification number1


# Part 2: parameters, initial states and time measurement ---------------------------------------

## You need to fill the following section.
## Easiest way is to first eval "model_extract()" and copy paste the output
## To directly have pre-fille the right parameter defaults values and initial states
## You can of course do modifications if needed

parameters_default_values <- c(
  
  ratioTumor =  1,
  ke_BCl2_I = 0,
  ke_BClxl_I = 0,
  ke_Mcl1_I = 0,
  k2_Bcl2_I = 10,
  k2_Bclxl_I = 10,
  k2_Mcl1_I = 10,
  switchE = 1,
  kelimBAXBAX = 0.014,
  hillA15 = 5,
  EC50A15 = 1E-6
  
) # paremeters value to be used by default (NA if you want to force providing value)
initial_cmt_values <- c(
  
  Venetoclax_Bcl2 = 0
  
) # initial compartment values. At least one, and every missing cmt name would be set to 0

times <- c(0,100,200,400,600,680,690,seq(700,750,1)) # times you want to see your observations

## Verification: perform verification number2


# Part 3: determine the fate of the cell ----------------------------------

# res will be the results of the simulations (do "res <- simulations()" to help fill the file)
# criteria should be an expression working with res with final output being a TRUE of FALSE
# TRUE being cell death, FALSE being cell survival

criteria <- expr ({
  
  max(res$TimeAbove) > 0.1
  
})

## Verification: perform verification number3
protocols <- list( dose0 = tibble(cmt = "Venetoclax", time = 0, amt = 0))
# parametre / conc qui s'ils sont augmenté et mort alors morts, si diminué et survie alors survie
param_death <- c("A1210477_0", "A1155463_0", "Venetoclax_0", "BIM0", "PUMA0", "NOXA0" )
param_increase <- list(TimeAbove = c("Bcl2_I0", "Bclxl_I0", "Mcl1_I0", "BIM0", "PUMA0", "NOXA0" ))
# parametre / conc qui s'ils sont augmenté et mort alors morts, si diminué et survie alors survie
param_survive <- c("Bcl20", "Bclxl0", "Mcl10") # c("ke_Venetoclax")
param_reduce  <- list(TimeAbove = c("Bcl20", "Bclxl0", "Mcl10")) #


param_no_impact <- list(TimeAbove = character())
# Part 4: data and concentration to test -------------------------------------

# Part 4: data and concentration to test -------------------------------------

# data shoud have at least ID, Value and concX columns, X being replace by drug number (one col by drug concentration)
# Avoid any column starting with "conc" if it is not a drug concentration / dose column

data_VT <- read.table("D:/these/Second_project/QSP/VirtualTumor/datademo_rework.csv", sep = ";", header = T) %>%
  as_tibble


# Normally the following code will automatically detect numbers of drug and
# extract the different concentration levels

ndrug <- sum(grepl("^conc", names(data_VT)))

for(a in 1:ndrug){

drug <- paste0("conc", a)

expr(!!parse_expr(drug) <<- unique(data_VT[[!!drug]])) %>%
  eval

}

# You can always modify manually with such code
# ndrug <- 4
# conc1 <- c(0,0.08,0.16,0.32,0.64,1.3,2.60,5,10,20)
# conc2 <- c(0,0.08,0.16,0.32,0.64,1.3,2.60,5,10,20)
# conc3 <- c(0,0.08,0.16,0.32,0.64,1.3,2.60,5,10,20)
# conc4 <- c(0,5,10,15)





