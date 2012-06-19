#include <cstdio>
#include <map>
#include <utility>
#include <cstdlib>
#include <algorithm>
#include <vector>
using namespace std;

#define MP(a,b) make_pair(a,b)

// system parameters -- VERY IMPORTANT
const int CELLS               = 1000;           // cells in the pool
const int MOI                 = 10;             // multiplicity of infection
const int PHAGES              = CELLS * MOI;    // phages in the pool
const int PHAGES_BORN         = 0.01 * PHAGES;  // phages born from a single killed cell
const double NEW_SPACER_PROB  = 0.01;           // probability of obtaining a new spacer when a cell kills the attacking virus
const int STEPS               = 100000*CELLS;   // number of phage-cell attacks

// resistant phage mutation
double mutationProb[]         = { 0.001 };      // probability of mutation on each phage attack
double mutantKillness[]       = { 1.0 };        // probability of killing a cell which doesn't have the corresponding spacer

// the source of randomness
inline double randf() { return (double)rand() / RAND_MAX; }
inline int    randi(int max_val) { return rand() % max_val; }

// the state defines the system in each moment of time
class State {
 public:

  int nClades;            // number of different cell clades
  int nMutants;           // number of different phage mutants

  int clade[CELLS];       // 'clade[c]' is the clade of the cell 'c'
  int mutant[PHAGES];     // 'mutant[p]' is the mutant type of the phage 'p'

  map< pair<int,int>, double > KillProb;    // KillProb[ pair(clade,mutant) ]

  // methods
 public:
  void copy_mutant    (int mutant_curr, int mutant_new);
  void copy_clade     (int clade_curr, int clade_new);
  void create_mutant  (int phage, int cell, double killness);
  void create_clade   (int phage, int cell);
  void kill_cell      (int phage, int curr_cell);
  void kill_phage     (int phage);
} S;

// copies or makes a new phage mutant
void State::copy_mutant (int mutant_curr, int mutant_new)
{
  for(int cl=0; cl<nClades; cl++)
    KillProb[ MP(cl, mutant_new) ] = KillProb[ MP(cl, mutant_curr) ];
}

// copies or makes a new cell clade
void State::copy_clade (int clade_curr, int clade_new)
{
  for(int m=0; m<nMutants; m++)
    KillProb[ MP(clade_new, m) ] = KillProb[ MP(clade_curr, m) ];
}

// creates a new phage mutant and changes the type of mutant of the current phage to the new type
void State::create_mutant (int phage, int cell, double killness)
{
  int mutant_curr = mutant[phage];
  int mutant_new = nMutants;

  copy_mutant(mutant_curr, mutant_new);

  KillProb[ MP(S.clade[cell], mutant_new) ] = killness;
  
  mutant[phage] = mutant_new;
  nMutants++;
}

// creates a new cell clade and changes the clade of the current cell to the new one
void State::create_clade (int phage, int cell)
{
  int clade_curr = clade[cell];
  int clade_new = nClades;
  
  S.copy_clade(clade_curr, clade_new);

  KillProb[ MP(mutant[phage], clade_new) ] = 0;

  clade[cell] = clade_new;
  nClades++;
}

// replaces the current cell with a random other cell
// replaces random other phages with the current phage
void State::kill_cell (int phage, int cell)
{
  // a random cell replicates substituting the current cell
  int replicated_cell = (cell + randi(CELLS-1)) % CELLS;
  clade[cell] = clade[replicated_cell];
  //S.copy_clade(clade[replicated_cell], clade[cell]);

  // new PHAGES_BORN are born substituting some random PHAGES_BORN phages
  int new_phage[PHAGES_BORN];
  for(int p=0; p<PHAGES_BORN; p++) {
    new_phage[p] = randi(PHAGES);

    for(int q=0; q<p; q++)
      if(new_phage[q] == new_phage[p]) {
        p--; // generate a new random value
        break;
      }

    mutant[ new_phage[p] ] = mutant[phage];
  }
}

// replace the current phage with a random one
void State::kill_phage (int phage)
{
  int replicated_phage = (phage + randi(PHAGES-1)) % PHAGES;

  mutant[phage] = mutant[replicated_phage];
}

// initialization
void init ()
{
  S.nClades = 1;
  S.nMutants = 1;

  for(int c=0; c<CELLS; c++)
    S.clade[c] = 0;
  
  for(int p=0; p<CELLS; p++)
    S.mutant[p] = 0;

  S.KillProb[ MP(0,0) ] = 0.0;
}

// a phage attacks a cell
void attack (int phage, int cell)
{
  // mutations with only one level of killness
  if (randf() <= mutationProb[0]) {
    // 2)
    S.create_mutant(phage, cell, mutantKillness[0]);
    printf("create_mutant\n");
  }

  // whether the phage kills the cell
  if (randf() <= S.KillProb[ MP(S.mutant[phage],S.clade[cell]) ])
    // 3)
    S.kill_cell(phage, cell);
  else {
    // 4.1)
    if (randf() <= NEW_SPACER_PROB) {
      S.create_clade(phage, cell);
      printf("create_clade\n");
    }

    // 4.2)
    S.kill_phage(phage);
    //printf("phage killed\n");
  }
}

void stat (int step)
{
  vector<int> cladeHist(S.nClades);    // cladeHist[i] -- number of cells of clade 'i'  

  for(int c=0; c<CELLS; c++)
    ++cladeHist[ S.clade[c] ];

  //sort(cladeHist.begin(), cladeHist.end());

  int livingClades=0;
  for(int i=0; i<cladeHist.size(); i++)
    if (cladeHist[i])
      livingClades++;

  vector<int> mutantHist(S.nMutants);    // cladeHist[i] -- number of cells of clade 'i'  
  
  for(int m=0; m<PHAGES; m++)
    ++mutantHist[ S.mutant[m] ];
  
  int livingMutants=0;
  for(int i=0; i<mutantHist.size(); i++)
    if (mutantHist[i])
      livingMutants++;
  
//  printf("step #%d:\t%d clades(%d live)\t%d mutants(%d live)\n", step, S.nClades, livingClades, S.nMutants, livingMutants);
}

// good luck in live, cells and phages!
void live ()
{
  for(int step=1; step<=STEPS; step++) {
    //int c = randi(CELLS);
    int c = (step-1) % CELLS;
    int p = randi(PHAGES);

    attack(p, c);
    stat(step);
  }
}

int main ()
{
  srand(0);

  init();
  live();

  return 0;
}
