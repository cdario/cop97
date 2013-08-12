//============================================================================
// Name        : cop97.cpp
// Author      : cdario@github
// Version     : 1.0
// Copyright   : not set
// Description : My first Implementation of cop97 
// Location    : https://github.com/cdario/cop97
//============================================================================

#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>

using namespace std;

enum PIEnum {
  QUEUES, STOPS, DELAY
};
const int PI = DELAY;

const int red = 1;
const int mingreen = 2;
const int T = 10; //planning horizon or no. of discrete time steps
const int M = 8; //predetermined number of phases to compute -
static const char phaseSeq[] = {'A', 'B', 'C'};
vector<int> phases(phaseSeq, phaseSeq + sizeof (phaseSeq) /
		   sizeof (phaseSeq[0])); // A = 0, B = 1, C = 2

const int initialPhase = 2;
int idxCurrentPh = initialPhase; // set initial phase to C (2)
int permanentQueueLengths[T] = {};
int arrivalData[10][3] = {
  { 0, 0, 1},
  { 0, 0, 1},
  { 0, 0, 0},
  { 0, 1, 0},
  { 0, 1, 0},
  { 1, 1, 0},
  { 1, 0, 0},
  { 1, 1, 0},
  { 0, 1, 0},
  { 0, 0, 0}
};

std::vector< std::vector<int> > v; //v_j(s_j);
std::vector< std::vector<int> > x_star; // optimal solutions x*_j(s_j)
vector<vector<vector<int> > > Q; // permanent queue lengths Q_{phi, j}(s_j)
vector<vector<vector<int> > > L; // temporary queue lengths Q_{phi, j}(s_j, x_j)
vector<vector<vector<int> > > S; // temporary stopped  L_{sigma, j}(s_j, x_j)

vector<int> getFeasibleGreens(int sj, int j);
int getQ(int sj, int ph, int j);
int getArrivals(int a, int b, int phi);
int getB(int a, int b, int phi);
int getM(int phi, int xj); //simplified
int getT(int d, int phi); //simplified
int getSaturationFlow(int phi); //not implemented

void initMatrix_v(int init);
//void updateIndexPhase2Stage();
void printVector(vector<int> values);
void printMatrix(vector<vector<int> > values);
void printArray(int arry[], int sz);

void RunCOP();

int main(int argc, char** argv) {
  RunCOP();
  return 0;
}

void RunCOP() {

  std::vector<int> X[T] = {};

  v.resize(M);
  x_star.resize(M);

  for (int i = 0; i < M; ++i) {
    v[i].resize(T);
    x_star[i].resize(T);
  }

  initMatrix_v(-1);
  cout << "\n**cop97 started \n ... Matrix v initialised\n";
  cout << "\n**Computing DP stages";

  int j = 1;
  bool goOn = 1;

  do {

    // <editor-fold defaultstate="collapsed" desc="header stage">
    cout << endl << "\n\n\t\t\tStage " << j << " Calculations [" << phaseSeq[idxCurrentPh] << "]" << endl;
    cout << "--------------------------------------------------------------------" << endl;
    cout << "s" << j << "\tx*(s" << j << ")\tv(s" << j << ")\tQA\tQB\tQC\tXj(s" << j << ")\n";
    cout << "--------------------------------------------------------------------\n";
    // </editor-fold>

    for (int sj = red; sj <= T; sj++) {

      cout << " " << sj;

      X[j] = getFeasibleGreens(sj, j);
      int xSz = X[j].size();

      // <editor-fold defaultstate="collapsed" desc="header value function calculations">
      //                        cout << endl << "\n\n\t\t Value Function Calculations for Stage " << j << endl;
      //                        cout << "--------------------------------------------------------------------" << endl;
      //                        cout << "x" << j << "\tfj(sj, xj)\tsj\tvj(sj)\t\tfj + v(j-1)\n";
      //                        cout << "--------------------------------------------------------------------\n";
      // </editor-fold>

      //dimensions L and S : [ sj ] [ xj ] [ phi ]
      L.resize(T);
      S.resize(T);
      for (int i = 0; i < T; ++i) {
	L[i].resize(xSz);
	S[i].resize(xSz);

	for (int j = 0; j < xSz; ++j) {
	  L[i][j].resize(phases.size());
	  S[i][j].resize(phases.size());
	}
      }

      //dimensions Q : [ sj ] [ phi ] [ j ]
      Q.resize(T);
      for (int i = 0; i < T; ++i) {
	Q[i].resize(phases.size());

	for (int j = 0; j < phases.size(); ++j)
	  Q[i][j].resize(M);
      }

      int index_xj = 0;
      int currentValueFn = -1;
      int minValueFn = 99999; // arbitrary large value
      int optimal_x = -1;
      int optimal_index_x = -1;
      //int earliestArrival = 0;  // not used due to simplification

      for (vector<int>::iterator it = X[j].begin(); it != X[j].end(); ++it) {
	int xj = *it;

	int hj = (xj!=0) ? (xj+red) : 0; //transition value

	// at stage 0, no steps allocated
	int si = (j!=1) ? (sj-hj) : 0; // si equals s_{j-1}

	int tQueue = 0;
	int tStops = 0;
	int tDelay = 0;

	int index_sj = sj - red;
	int index_maxPh = -1;

	//performance index calculation Max Q Length
	// which phase has the longest temp queue?
	int pi_MaxQ = -1;
	int pi_NumStops = 0;
	int pi_Delay = 0;

	for (int index_p = 0; index_p < phases.size(); index_p++) {

	  if (index_p != idxCurrentPh) // phase w/o right-of-way
	    {

	      int arrival = arrivalData[sj][index_p];

	      /* //unused due to simplification
		 if (earliestArrival == 0 && arrival != 0)
		 earliestArrival = sj;
	      //*/

	      // temporary queues
	      tQueue = getQ(si, index_p, j - 1)
		+ getArrivals(si, sj, index_p);

	      // temporary stops
	      tStops = getArrivals(si, sj, index_p);

	      //delay
	      tDelay = getQ(si, index_p, j - 1)*(sj - si)
		+ getB(si, sj, index_p);

	    } else { //phase with right-of-way

	    // temporary queues
	    int queueTerm = getQ(si, idxCurrentPh, j - 1)
	      + getArrivals(si, si + xj, idxCurrentPh)
	      - getM(idxCurrentPh, xj);

	    tQueue = max(0, queueTerm)
	      + getArrivals(si + xj, sj, idxCurrentPh);

	    // temporary stops
	    int stopsTerm =
	      getArrivals(si, si + xj, idxCurrentPh)
	      - max(0, getM(idxCurrentPh, xj)
		    - getQ(si, idxCurrentPh, j - 1));

	    tStops = max(0, stopsTerm)
	      + getArrivals(si + xj, sj, idxCurrentPh);

	    /*
	     * t_{phi}:
	     * This function of s_j and x_j is only defined for phi = phi(j)
	     *
	     * Denotes the arrival time of the earliest vehicle that is required to stop
	     * while phase phi(j) has the right of way.
	     *
	     * When no arriving vehicles are stopped we put :
	     * t_{phi} = s_{j-1} + x_j, where phi = phi(j)
	     *
	     * comment: this function represents the time when a vehicle arrived and
	     * needs to stop, considering the maximum discharge for that phase (M/T),
	     * in the case of instantaneous queue clearance T(d) = 0 for all d,
	     * No arriving vehicles are stopped, as they are discharged instantaneously.
	     *
	     */
	    int tp = si + xj; // equals s_{j} - red

	    /*      //unused to simplification
		    if if (earliestArrival != 0)
		    tp = earliestArrival;
	    //*/

	    // delay
	    int delayTerm = min(getQ(si, idxCurrentPh, j - 1),
				getM(idxCurrentPh, xj));

	    tDelay = getT(delayTerm, idxCurrentPh)
	      + max(0, getQ(si, idxCurrentPh, j - 1) -
		    getM(idxCurrentPh, xj))*(sj - si)
	      + getB(tp, sj, idxCurrentPh);

	  }

	  // record temporary queue lengths and stops
	  L[index_sj][index_xj][index_p] = tQueue;
	  S[index_sj][index_xj][index_p] = tStops;

	  // PI Max Queue : use operator max
	  if (tQueue > pi_MaxQ) {
	    pi_MaxQ = tQueue;
	    index_maxPh = index_p;
	  }

	  // PI Stops & Delay : use operator +
	  pi_NumStops += tStops;
	  pi_Delay += tDelay;

	} //end phaseSequence cycle

	// index fix
	if (j != 1)
	  si -= red;

	switch (PI) {
	case QUEUES:
	  currentValueFn = max(pi_MaxQ, v[j - 1][si]);
	  break;
	case STOPS:
	  currentValueFn = pi_NumStops + v[j - 1][si];
	  break;
	case DELAY:
	  currentValueFn = pi_Delay + v[j - 1][si];
	  break;
	}

	//minimisation v_j : keep minimum value
	if (minValueFn > currentValueFn) {
	  minValueFn = currentValueFn;
	  optimal_x = xj;
	  optimal_index_x = index_xj;
	}

	//<editor-fold defaultstate="collapsed" desc="body value fn calculations">

	//                cout << "" << xj;
	//                cout << "\t" << maxQ;
	//                cout << "\t\t" << si;
	//                cout << "\t" << v[j - 1][si];
	//                cout << "\t\t" << vj << endl;

	//</editor-fold>

	index_xj++;
      } // end X[j] cycle

      // sj - red :  adjust value to column index
      v[j][sj - red] = minValueFn;
      x_star[j][sj - red] = optimal_x;

      // -red and -1 deal, reconcile indices

      int optIndeX = 0; // stage 1 simplification
      if (j != 1)
	optIndeX = optimal_index_x;

      // temporary to permanent queue lengths
      for (int pp = 0; pp < phases.size(); pp++)
	Q[sj - red][pp][j - 1] = L[sj - red][optIndeX][pp]; // -1 :index

      //<editor-fold defaultstate="collapsed" desc="body stage">
      cout << "\t" << optimal_x;
      cout << "\t" << v[j][sj - red];

      for (int pp = 0; pp < phases.size(); pp++)
	cout << "\t" << Q[sj - red][pp][j - 1];

      cout << setfill(' ') << setw(30 - 2 * L.size());
      printVector(X[j]);
      cout << endl;
      // </editor-fold>

      if (sj % 2 == 0)
	cout << endl;
    } //end sj cycle


    //stopping criterion

    if (j >= phases.size()) {
      for (int k = 1; k <= phases.size() - 1; k++) {
	goOn = goOn && (v[j - k][T- red] == v[j][T - red]);
      }

      goOn = !goOn;

      if (goOn)
	{
	  //Updates index of current phase in cycles
	  idxCurrentPh = idxCurrentPh==2 ? 0:idxCurrentPh + 1;
	  j++;
	}
    }
    else
      {
	idxCurrentPh = idxCurrentPh==2 ? 0:idxCurrentPh + 1;
	j++;
      }

  } while (goOn);

  cout << "\n**Stopping Criterion satisfied\n **Printing states and optimal control tables";
  cout << "\n\nTable states( v )\n\n";
  printMatrix(v);
  cout << "\nTable optimal control associated (x^star)\n\n";
  printMatrix(x_star);



  /*  Retrieval of Optimal Policy     */

  int jsize = j - (phases.size() - 1);
  //std::vector<int> s;
  int s_star= T;
  cout <<  "\n**Recovering optimal solution";
  cout <<  "\n\nMin control delay sequence ";


  int idxSeq = initialPhase;

  for(int jj= jsize; jj>=1; jj--)
    {
      int xx = x_star[jj][s_star-red];
      cout << phaseSeq[idxSeq] << ":" <<xx <<" ";

      if (jj > 1)
	{
	  int hj_star = (xx!=0) ? (xx+red) : 0; //transition value
	  s_star = s_star - hj_star;
	}

      idxSeq = idxSeq==2 ? 0:idxSeq + 1;

    }

  cout  << "\n **cop97 ended \n";

};


/**
 *  @brief  Return the set of feasible greens for state %sj
 *  @param sj number of steps currently assigned
 *  @param j  phase of the dynamic programme
 *  @return  %vector of %int values containing the feasible set
 *
 *   This function computes the set of greens given the state, minimum green
 *  and red times
 */

std::vector<int> getFeasibleGreens(int sj, int j) {

  std::vector< int > set;

  if (j == 1) //simplification removes the min green restriction for stage 1
    set.push_back(sj - red);
  else {
    set.push_back(0);
    int c = mingreen;

    if (!(sj - red < mingreen)) {
      do {
	set.push_back(c);
	c++;
      } while (c < (sj - red));
    }
  }

  return set;
};

/**
 *  @brief  Retrieves arrival data within a time interval
 *  @param  a time interval lower bound
 *  @param  b time interval upper bound
 *  @param  phi current phase
 *  @return number of arrivals
 *
 *  Sums arrivals within the interval for the given phase
 *
 */

int getArrivals(int a, int b, int phi) {

  if (a == b)
    return 0;

  int vehicles = 0;
  int cc = a;

  do {
    vehicles += arrivalData[cc][phi];
    cc++;
  } while (cc < b); // for [a,b), with a!=b

  return vehicles;
}

/**
 *  @brief  Retrieves min delay within a time interval
 *  @param  a time interval lower bound
 *  @param  b time interval upper bound
 *  @param  phi phase requested
 *  @return min delay within the interval
 *
 *  Sums arrivals within the interval [a, b)
 *
 */

int getB(int a, int b, int phi) {

  int dB = 0;
  int k = a; 

  do {
    if (arrivalData[k][phi] != 0) //if vehicle k requests phase phi
      dB += b - k;
    k++;
  } while (k < b); //a <= ak < b

  return dB;
}

/**
 *  @brief  maximum vehicle discharge give a green time
 *  @param  phi phase index to be discharged
 *  @param  xj green time
 *  @return maximum vehicle discharge given green

 *  Finds the maximum number of cars that can be discharged for a green
 *  time based on the phase saturation flow
 */

int getM(int phi, int xj) {

  /*  TODO: simplicity assumption: M_phi(x) = INF for all x > 0;
   *   i.e. instantaneous queue clearance
   */

  if (xj == 0)
    return 0;

  return 100000;
}

/**
 *  @brief  gets value from structure of permanent queue lengths
 *  @param  sj steps currently assigned
 *  @param  ph phase of this queue
 *  @param  j programme stage
 *  @return permanent queue length
 *
 *  Retrieves permanent queue for stage, time step and phase, assuming initial
 *  queues are zero
 */

int getQ(int sj, int ph, int j) {

  if (sj == 0)
    return 0; //assuming initial queues are zero

  //index fix
  return Q[sj - red][ph][j - 1];
}

/**
 *  @brief  function of the number of vehicles discharged
 *  @param  d number of vehicles to be discharged
 *  @param  ph phase for which vehicles are discharged
 *  @return total delay in discharging %d vehicles

 *  Computes the delay of discharging vehicles based on the phase and the
 *  saturation flow for the given phase. Returns zero, with the
 *  instantaneous queue clearance = 0 assumption
 */
int getT(int d, int phi) {

  // example assumption, instantaneous queue clearance = 0, for all d
  return 0;

  // TODO: implement return d/getSaturationFlow(phi);
}

/**
 *  @brief  obtain the saturation flow for a given phase
 *  @param  ph phase of the saturation flow
 *  @return saturation flow value (vehicles per-hour per-lane)

 *  Obtains saturation flow to be used in function getM and getT, if no
 *  instantaneous clearance assumption is implemented
 */
int getSaturationFlow(int phi) {
  return 0;
}

/*
************** other functions
*/

/**
 *  @brief  initialise matrix v to the given value
 *  @param  init %int value to use to initialise the matrix
 *
 *  appropriate matrix v initialisation
 */

void initMatrix_v(int init) {
  for (int i = 0; i < v.size(); ++i) {
    for (int j = 0; j < v[i].size(); ++j) {
      if (i == 0) {
	v[i][j] = 0;
      } else {
	v[i][j] = init;
      }
    }
  }
}

void printArray(int arry[], int sz) {

  cout << "[ ";
  for (int i = sz - 1; i >= 0; i--) {
    cout << arry[i] << " ";
  }
  cout << "]";
}

void printVector(vector<int> values) {

  cout << "[ ";
  for (vector<int>::iterator i = values.begin(); i != values.end(); ++i) {
    cout << *i << " ";
  }
  cout << "]";
}

void printMatrix(vector<vector<int> > values) {

  for (int i = 0; i < values.size(); ++i) {
    printVector(values[i]);
    cout << endl;
  }
}

