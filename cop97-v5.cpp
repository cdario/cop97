#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

enum PIEnum {
  QUEUES, STOPS, DELAY
};
const int PI = DELAY;

const int red = 1;
const int mingreen = 2;
const int T = 10; //planning horizon
const int M = 8; //predetermined number of phases to compute -
static const char phaseSeq[] = {'A', 'B', 'C'};
vector<int> phases(phaseSeq, phaseSeq + sizeof (phaseSeq) /
		   sizeof (phaseSeq[0])); // A = 0, B = 1, C = 2

int initialPhase; // = 2;
int idxCurrentPh; //= initialPhase; // set initial phase to C (2)
int permanentQueueLengths[T] = {};
int arrivalData[10][3];

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
int getArrivalEarliest(int si, int sj, int xj, int phi); //NEW


void initMatrices(int init);
void printVector(vector<int> values);
void printMatrix(vector<vector<int> > values);
void printArray(int arry[], int sz);
void printArrivals();

void RunCOP();
void loadFromFile(char* filename);

int main(int argc, char** argv) {

  loadFromFile(argv[1]);
  initialPhase  = *argv[2]- 48;
  idxCurrentPh = initialPhase;

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

  initMatrices(-1);
  unsigned int j = 1;
  bool criterion_flag = 1;

  do {

          // <editor-fold defaultstate="collapsed" desc="header stage">
      cout << endl << "\n\t\t\tStage " << j << " Calculations [" << phaseSeq[idxCurrentPh] << "]" << endl;
      cout << "--------------------------------------------------------------------" << endl;
      cout << "s" << j << "\tx*(s" << j << ")\tv(s" << j << ")\tQA\tQB\tQC\tXj(s" << j << ")\n";
      cout << "--------------------------------------------------------------------\n";
      // </editor-fold>

    for (int sj = red; sj <= T; sj++) {

         cout << " " << sj;


      X[j] = getFeasibleGreens(sj, j);
      int xSz = X[j].size();

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

     Q.resize(T);
     for (unsigned int i = 0; i < T; ++i) {
      Q[i].resize(phases.size());

      for (unsigned int j = 0; j < phases.size(); ++j)
        Q[i][j].resize(M);
    }

    int index_xj = 0;
    int currentValueFn = -1;
    int minValueFn = 99999;
    int optimal_x = -1;
    int optimal_index_x = -1;

    for (vector<int>::iterator it = X[j].begin(); it != X[j].end(); ++it) {
     int xj = *it;

	int hj = (xj!=0) ? (xj+red) : 0; //transition value

	// at stage 0, no steps allocated
	int si = (j!=1) ? (sj-hj) : 0; // si equals s_{j-1}
//int si = sj -hj;
//  cout << "\n hj, si: " <<hj << ", " << si<< "\n";

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

	for (unsigned int index_p = 0; index_p < phases.size(); index_p++) {

	  if (index_p != idxCurrentPh) // phase w/o right-of-way
   {

     int arrival = arrivalData[sj][index_p];


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


       //NEW
       /*
        Calculate tp, function of sj and xj
       */
       int tp = getArrivalEarliest(si, sj, xj, idxCurrentPh);

	//    int tp = si + xj; // equals s_{j} - red
  //    cout <<"(tp: "<< tp<<")";

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

	// index fix TODO: implications
	if (j != 1 && si >= red){
   si -= red;
   //cout << "   si = " << si << " \n"; 
 }

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
     for (unsigned int pp = 0; pp < phases.size(); pp++)
      	Q[sj - red][pp][j - 1] = L[sj - red][optIndeX][pp]; // -1 :index

      /**print*************************/

cout << "\t" << optimal_x;
        cout << "\t" << v[j][sj - red];

        for (unsigned int pp = 0; pp < phases.size(); pp++)
    cout << "\t" << Q[sj - red][pp][j - 1];

        cout << setfill(' ') << setw(30 - 2 * L.size());
        printVector(X[j]);
        cout << endl;
        // </editor-fold>

        if (sj % 2 == 0)
    cout << endl;


          /**print*************************/



    } //end sj cycle

    //************ STOPPING CRITERION ***********
//if(criterion_flag)
 //{ 
    if (j >= phases.size()) {
      for (unsigned int k = 1; k <= phases.size() - 1; k++) {
      	criterion_flag = criterion_flag && (v[j - k][T- red] == v[j][T - red]);
      }

      criterion_flag = !criterion_flag;
idxCurrentPh = idxCurrentPh==2 ? 0:idxCurrentPh + 1;
      if (criterion_flag)
      {
	  //Updates index of current phase in cycles
       
       j++;
     }
   }
   else
   {
     idxCurrentPh = idxCurrentPh==2 ? 0:idxCurrentPh + 1;
     j++;
   }
//}
 } while (criterion_flag);

     cout << "\n**Stopping Criterion satisfied\n **Printing states and optimal control tables";
    cout << "\n\nTable states( v )\n\n";
    printMatrix(v);
    cout << "\nTable optimal control associated (x^star)\n\n";
    printMatrix(x_star);

  /*  Retrieval of Optimal Policy     */
// cout << endl << "j :"<< j  <<endl;

 int jsize = j - (phases.size() - 1);
 int s_star= T;

 //cout << endl << "jsize :"<< jsize  <<endl;
 //cout << endl << "s* :"<< s_star  <<endl;

//new
 //int idxSeq = initialPhase;
int idxSeq = idxCurrentPh;
//string controlSeq = "[ ";
 //cout << "[ ";
    //int optimalControlSeq [jsize];
  int* optimalControlSeq = new int[jsize];


  for(int jj= jsize; jj>=1; jj--)
  {
    int xx = x_star[jj][s_star-red];

    //cout <<"\njj, s_star-red = "<< jj <<", "<<s_star-red << " : " << xx;
    optimalControlSeq[jj-1] = xx;

    if (jj > 1) {
     int hj_star = (xx!=0) ? (xx+red) : 0; 
     s_star = s_star - hj_star;
   }

//TODO: FISHY,  BACKTRACING? SHOULDN'T RUN BACKWARDS
   //idxSeq = idxSeq==2 ? 0:idxSeq + 1;
   idxSeq = idxSeq==0 ? 2:idxSeq - 1;

 }
 //cout << "]\n";

cout << endl; 
printArray(optimalControlSeq, jsize);
delete optimalControlSeq;
cout << endl;
printArrivals();

}; /**************** END MAIN*************/




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

//arrival time of the earliest vehicle required to stop when phi(j);
int getArrivalEarliest(int si, int sj, int xj, int phi)
{
  int timeArrival = 999;
  /*    //NEW
  for (int p=0; p< phases.size(); p++)
  {
    if (p != phi)
    {
      for (int i = sj; i >= 0; i--)
      { 
         if(arrivalData[i][p]!=0)
          {
            //if (timeArrival >= i)
              timeArrival = i;
          }
          else
          {
            break;
          }
      }
    } 
  }

  if (timeArrival == 999)  // */
  // no stops, similar for A  calculations
  timeArrival = si + xj;

  return timeArrival;
}

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

int getM(int phi, int xj) {

  /*  TODO: simplicity assumption: M_phi(x) = INF for all x > 0;
   *   i.e. instantaneous queue clearance
   */

   if (xj == 0)
    return 0;

  return 100000;
}

int getQ(int sj, int ph, int j) {

  if (sj == 0)
    return 0; //assuming initial queues are zero

  //index fix
  return Q[sj - red][ph][j - 1];
}

int getT(int d, int phi) {

  // example assumption, instantaneous queue clearance = 0, for all d
  return 0;

  // TODO: implement return d/getSaturationFlow(phi);
}

int getSaturationFlow(int phi) {
  return 0;
}

/*
************** OTHER FUNCTIONS
*/

void initMatrices(int init) {
  for (unsigned int i = 0; i < v.size(); ++i) {
    for (unsigned int j = 0; j < v[i].size(); ++j) {
      if (i == 0) {
       x_star[i][j] = v[i][j] = 0;

     } else {
       x_star[i][j] = v[i][j] = init;
     }
   }
 }
}

void printArray(int arry[], int sz) {

cout << "[ ";
  for (int i = 0; i < sz; i++) {
    cout << phaseSeq[(i+initialPhase)%phases.size()]<<":"<< arry[i] << " ";
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

for (unsigned int i = 0; i < values.size(); ++i) {
  printVector(values[i]);
  cout << endl;
}
}

void printArrivals() {
cout << endl; //TODO: dynamic sizes
for (unsigned int i = 0; i < 10; ++i) {
  for (unsigned int j = 0; j < 3; ++j) {
      cout << " " << arrivalData[i][j];
  }
  cout << endl;
}
cout << endl;
}


void loadFromFile(char* filename) {
int x, y;
ifstream in(filename);

if (!in) {
  cout << "Cannot open file.\n";
  return;
}

for (y = 0; y < 10; y++) {
  for (x = 0; x < 3; x++) {

   in >> arrivalData[y][x];
 }
}
in.close();
}