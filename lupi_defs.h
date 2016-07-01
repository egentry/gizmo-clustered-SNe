#ifdef GRACKLE_OPTS
#define LUPI_SOLAR_MET grackle_data.SolarMetalFractionByMass
#define NUM_METAL_SPECIES 1
#endif

#ifdef GALSF_FB_LUPI
#define E_SN 			1.0e51
#define CHABRIER_MIN_MASS 	0.1 //Msun
#define CHABRIER_MAX_MASS 	100 //Msun
#define CHABRIER_NORM_PL 	0.238378 //(4.46556e-2/0.187331) HighMass Coeff/Global Normalisation

#define SNII_MIN_MASS		8.0 //Msun
#define SNII_MAX_MASS 		40.0 //Msun
#define SNI_MIN_AGE 		1.0e2 //Myr
#define SNI_MAX_AGE 		1.0e4 //Myr

#define SNI_STARMASS 		1.3e-3
#define MSNI_STARMASS 		5.97837e-3
#define SNI_NORM_MYR 		2.1714724e-1

#endif


#ifdef BH_LUPI
#define EDD_NORM 		(4*M_PI*C*GRAVITY*PROTONMASS/THOMPSON)
#define N_EDDINGTON		500 //Maximum eddington ratio allowed (for ETA_SLIMDISC model)
#define BH_FBK_COUPLE 		0.15
#define LUPI_BH_SPIN		0.9
#endif
