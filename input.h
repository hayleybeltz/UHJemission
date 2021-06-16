/*----------------------- input.h ------------------------------------

Author: Eliza Miller-Ricci (emillerricci@cfa.harvard.edu)
Last modified:

---------------------------------------------------------------------- */

#ifndef __INPUT_H__
#define __INPUT_H__

/* Input choices. */

#define NTAU 500             /* Number of optical depth points in grid */
#define NTEMP 46             /* Number of temperature points in grid   */
#define NPRESSURE 28         /* Number of pressure points in grid      */
#define NLAMBDA 13443        /* Number of wavelength points in grid    */

#define G 6.83               /* Planet Surface Gravity                 */
#define R_STAR 1.204e+09     /* Stellar radius                         */
#define ORB_SEP 4.936e+09    /* Orbital separation [m]                 */
#define STELLAR_TEMP 6250    /* Stellar blackbody temperature          */
#define INCIDENT_FRAC 0.0    /* Incident fraction of starlight for     
				2-stream                               */

#define FORMAT 2             /* FORMAT=1 -> small opacity table        */
                             /* FORMAT=2 -> large opacity table        */

/* File names */

#define T_P_FILE "DATA/wasp76tp.dat"
#define CHEM_FILE "DATA/eos_solar_gas.dat"
#define OUTPUT_FILE "wasp76.dat"

#define CH4_FILE "DATA/opacCH4.dat"
#define CO2_FILE "DATA/opacCO2.dat"
#define CO_FILE "DATA/opacCO.dat"
#define H2O_FILE "DATA/opacH2O.dat"
#define NH3_FILE "DATA/opacNH3.dat"
#define C2H2_FILE "DATA/opacC2H2.dat"
#define H2S_FILE "DATA/opacH2S.dat"
#define HCN_FILE "DATA/opacHCN.dat"
#define K_FILE "DATA/opacK.dat"
#define Na_FILE "DATA/opacNa.dat"
#define PH3_FILE "DATA/opacPH3.dat"
#define TiO_FILE "DATA/opacTiO.dat"
#define VO_FILE "DATA/opacVO.dat"
#define CIA_FILE "DATA/opacCIA.dat"

#endif /* !__INPUT_H__ */

/* ------- end ---------------------------- input.h  ----------------- */


