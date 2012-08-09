// stochastic model of CRISPR immune system based on ODE model of
// Joshua Weitz et al. Multiscale Model of CRISPR-induced coevolutionary dynamics: Diversification at the interface of Lamarch and Darwin

ctmc

//------- model parameters from Joshua's paper --------

//const int ALLELES       = 5;            // of one strain
const int VIRUS_STRAINS = 10;            // of one strain
const int CELL_STRAINS  = 10;            // of one strain

// Molecular
const double  p       = 1e-5;           // CRISPR failure probability
const double  q       = 1e-5;           // New spacer acquisition probability 

// Ecological constants
const double  r       = 1.0;            // Growth rate [1/h]
const double  K       = 20; //316227.7;       // Carrying capacity [1/mL]
const int     beta    = 50;             // Burst size
const double  phi     = 1e-7;           // Adsorbtion rate [mL/h]     // TODO: check whether it is one and the same for each cells and viruses!
const double  m       = 0.00001;//0.1;            // Virus decay rate [1/h]

// Evolutionary
const double  mu      = 5e-7;           // Virus mutation rate [1/h]
//const double  rho     = 0.1;            // Density cutoff [1/mL]

//------------------------------------------------------

formula N = c0+c1;                      // number of all cells

formula new_cell_strain = min 'i' such that c[i]=0
formula new_virus_strain = min 'i' such that v[i]=0

module rates
  // for every virus 'v' & cell 'c'
  [lysis_when_not_immune] c_immune[v]=0 & c>0 & v>0 -> (1.0-q)*(1-c_immune[v])*c*v*phi  : true;
  [lysis_when_immune]     c_immune[v]=1 & c>0 & v>0 -> p*c_immune[v]*c*v*phi            : true; 

  // for every cell 'c'
  [div]  (N<K) & (c>0) ->   r*c*(1.0-N/K) : true;
  [mort] (N>K) & (c>0) -> - r*c*(1.0-N/K) : true:

  // for every virus 'v'
  [infection]     c>0 & v>0 -> phi*c*v : true;
  [deactivation]  v>0       -> m*v     : true;
endmodules

// array of cell strains
module cell_strain[CELL_STRAINS]
  c : int;                             // number of cells of this strain
  c_immune[VIRUS_STRAINS] : [0..1];    // whether it has a spacer against virus strain 'i'

  // cell division and mortality
  [div]  true -> (c'=c+1);
  [mort] c>0  -> (c'=c-1);
 
  // for every virus
  [lysis_when_not_immune] c>0 -> (c'=c-1);
  [lysis_when_immune]     c>0 -> (c'=c-1);

  // new strain possible for every virus 'k'
  [lysis_when_not_immune] v_al[k]=0 & new_cell_strain<CELL_STRAINS -> q :
      cell_strain[new_cell_strain].c'=1 &
      cell_strain[new_cell_strain].c_immune'=(c_immune xor (1<<k)) &
      c'=c-1;
endmodule

// array of virus strains
module virus_strain[VIRUS_STRAINS]
  v : int;                                          // number of viruses of this strain
  v_almost_immune_cell: [0..CELL_STRAINS-1];     // whether cell strain 'i' can obtain immunity to this virus strain by only one mutation 

  // for every cell
  [lysis_when_not_immune] v>0 -> (v'=v+beta);
  [lysis_when_immune]     v>0 -> (v'=v+beta);

  [infection]             v>0 -> (v'=v-1);
  [deactivation]          v>0 -> (v'=v-1);

  // mutation for every cell 'k'
  [lysis_when_immune] true -> mu :
      virus_strain[new_virus_strain].v'=1 &
      virus_strain[new_virus_strain].v_almost_immune_cell'=k &
      v'=v-1);
endmodule

