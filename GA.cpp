#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <time.h>
using namespace std;

const int cDebug = 0;

enum Xover { eRandom, eUniform, eOnePoint, eTwoPoint };

const Xover  CrossoverType = eTwoPoint;
const double cCrossoverRate = 0.95;
const double cMutationRate = 0.9;
const int    cNumGens = 1500;
const int    cPopSize = 100; // must be an even number
const int    cTournamentSize = 5;
const int    Seed = 1234;	//I replace Seed with time
const int    cTargetFitness = 10000; //desired distance of tour
int	cIndividualLength = 80;
int parSel = 0;	//parents selection : 0 - Tournament ; 1 - Roulette
int *longitude = NULL;
int *latitude = NULL;
int *cityType = NULL;
int *co = NULL; //uniq in crossover - the index of co is the city number
vector<int> co0;
vector<int> co2;
double weightTable[3][3] = {10, 7.5, 5, 
							7.5, 5, 2.5,
							5, 2.5, 1};
int distanceTable[1000][1000];		//distance table which is To avoid the large expense of 
//calculating the same distances during each fitness evaluation

void InitPop(int ***CrntPop, int ***NextPop, int **Fitness, int **BestMenber, double **rFitness);
void FreeMem(int **CrntPop, int **NextPop, int *Fitness, int *BestMember);
int Tournament(int *Fitness, int TournamentSize);
int Roulette(double *Fitness);
int EvaluateFitness(int *Member);
void Crossover(int *P1, int *P2, int *C1, int *C2);
void Copy(int *P1, int *P2, int *C1, int *C2);
void Mutate(int *Member);
double Rand01();    // 0..1
int RandInt(int n); // 0..n-1

int main(int argc, char *argv[]) {

	int **CrntPop, **NextPop; // the crnt & next population lives here
	// The possible longest distance between two city [0,999] - [999,0]. The distance is 1412.80.
	int *Fitness, BestFitness = 15000000, *BestMember; // fitness vars
	double *rFitness;
	int i, TargetReached = false;
	char fileName[100];

	//opening input file
	ifstream inFile;
	if(argc == 1){
		cout << "Please input file : ";
		cin >> fileName;
		inFile.open(fileName);
	}else{
		inFile.open(argv[1]);
	}
	if (!inFile) {
		cerr << "Unable to open file datafile.txt"; exit(1);
	}
	//Initiating arrays to save the positions and category of city
	inFile >> cIndividualLength;
	longitude = new int[cIndividualLength];
	latitude = new int[cIndividualLength];
	cityType = new int[cIndividualLength];
	co = new int[cIndividualLength];

	for(int i = 0; i < cIndividualLength; i++){
		inFile >> longitude[i] >> latitude[i] >> cityType[i];
	}
	inFile.close();

	//initiate the distanceTable
	for(int i = 0; i < cIndividualLength; i++){
		for(int j = 0; j < cIndividualLength; j++){
			int x = longitude[i] - longitude[j];
			int y = latitude[i] - latitude[j];
			double w = weightTable[cityType[i] - 1][cityType[j] - 1];
			distanceTable[i][j] = sqrt(x*x + y*y) * w;
		}
	}

	InitPop(&CrntPop, &NextPop, &Fitness, &BestMember, &rFitness);
	for (int Gen = 0; Gen<cNumGens; Gen++) {
		for (i = 0; i<cPopSize; i++) {
			// Evaluate the fitness of pop members
			Fitness[i] = EvaluateFitness(CrntPop[i]);
			if (BestFitness > Fitness[i]) { // save best member
				BestFitness = Fitness[i];
				for (int j = 0; j<cIndividualLength; j++)
					BestMember[j] = CrntPop[i][j];
				if (Fitness[i] <= cTargetFitness) {
					TargetReached = true;
					break;
				}
			}
		}
		if (TargetReached)break;
		
		//Calculate the Roulette wheel for each tour
		double WorstFitness = -1;
		long sumRFitness = 0;
		for (i = 0; i<cPopSize; i++) {
			if (WorstFitness < Fitness[i]) 
				WorstFitness = Fitness[i];
			rFitness[i] = Fitness[i];
		}
		for (i = 0; i<cPopSize; i++) {
			rFitness[i] -= WorstFitness;
			sumRFitness += rFitness[i];
		}
		for (i = 0; i<cPopSize; i++) {
			rFitness[i] /= sumRFitness;
		}

		// Produce the next population
		for (i = 0; i<cPopSize; i += 2) {
			int Parent1 = 0;
			int Parent2 = 0;

			if(parSel == 0){
				//Using tournament for parent selection
				Parent1 = Tournament(Fitness, cTournamentSize);
				Parent2 = Tournament(Fitness, cTournamentSize);
			}else{
				//Using roulette for parent selection
				Parent1 = Roulette(rFitness);
				Parent2 = Roulette(rFitness);
			}
			
			//Checking crossover
			if (cCrossoverRate>Rand01())
				Crossover(CrntPop[Parent1], CrntPop[Parent2], NextPop[i], NextPop[i + 1]);
			else
				Copy(CrntPop[Parent1], CrntPop[Parent2], NextPop[i], NextPop[i + 1]);
			if (cMutationRate<Rand01())Mutate(NextPop[i]);
			if (cMutationRate<Rand01())Mutate(NextPop[i + 1]);
		}
		int **Tmp = CrntPop; CrntPop = NextPop; NextPop = Tmp;
		//Print the fitness of every 20th generation 
		if(Gen % 20 == 0)
			cout << setw(3) << Gen << ':' << setw(5) << BestFitness << endl;
	}
	//Check if the target mentioned in line number 22 has been reached or not
	if (TargetReached) 
		cout << "Target fitness reached: " << BestFitness << " !\n";
	else 
		cout << "Target fitness not reached: " << BestFitness << " !\n";
	//Printing best path
	cout << "Best Individual: ";
	for (i = 0; i<cIndividualLength; i++)
		cout << BestMember[i]<<" "; 
	cout << endl;
	//Free memory
	FreeMem(CrntPop, NextPop, Fitness, BestMember);
	char s[20]; cin.getline(s, 20);
	return 0;
}

void InitPop(int ***CrntPop, int ***NextPop, int **Fitness, int **BestMember, double **rFitness) {
	int i, j, t, tmp;
	// load the initial pop members with random cities. 
	srand((int)time(NULL));
	*CrntPop = new int*[cPopSize];
	*NextPop = new int*[cPopSize];
	//Cities are represented with the number between 0 to popsize 
	for (i = 0; i<cPopSize; i++) {
		(*CrntPop)[i] = new int[cIndividualLength];
		(*NextPop)[i] = new int[cIndividualLength];
	}
	*Fitness = new int[cPopSize];
	*rFitness = new double[cPopSize];
	*BestMember = new int[cIndividualLength];
	//If fitness function can't be calculated or best member is still unknown
	if (Fitness == NULL || BestMember == NULL)exit(1);
	for (i = 0; i<cPopSize; i++) {
		for (j = 0; j<cIndividualLength; j++)
			(*CrntPop)[i][j] = j;
		for (j = 0; j<cIndividualLength; j++){
			tmp = RandInt(cIndividualLength);	//generate 0..cIndividualLength-1
			t = (*CrntPop)[i][j];
			(*CrntPop)[i][j] = (*CrntPop)[i][tmp];
			(*CrntPop)[i][tmp] = t;
		}
	}
}

void FreeMem(int **CrntPop, int **NextPop, int *Fitness, int *BestMenber) {
	for (int i = 0; i<cPopSize; i++) {
		delete[]CrntPop[i];
		delete[]NextPop[i];
	}
	delete CrntPop;
	delete NextPop;
	delete Fitness;
	delete BestMenber;
}

int EvaluateFitness(int *Member) {
	//Evaluate the fitness function from distance table
	int p1, p2;
	int TheFitness = 0;
	for(int i = 1; i < cIndividualLength; i++) {
		p1 = Member[i];
		p2 = Member[i-1];
		TheFitness += distanceTable[p1][p2];
	}
	return(TheFitness);
}

int Tournament(int *Fitness, int TournamentSize) {
	int WinFit = 15000000, Winner;
	//Perform tournament based parent selection
	for (int i = 0; i < TournamentSize; i++) {
		int j = RandInt(cPopSize);
		if (Fitness[j] < WinFit) {
			WinFit = Fitness[j];
			Winner = j;
		}
	}
	return Winner;
}

int Roulette(double *rFitness) {
	double RandomNumber = Rand01();
	double TempSum = 0;
	//Perform roulette based parent selection
	for(int i = 0; i < cPopSize; i++){
		TempSum += rFitness[i];
		if(TempSum > RandomNumber) return i;
	}
	return RandInt(cPopSize);
}

void Crossover(int *P1, int *P2, int *C1, int *C2) {
	int i, Left, Right;
	switch (CrossoverType) {
	case eRandom: // swap random genes
		for (i = 0; i<cIndividualLength; i++) {
			if (RandInt(2)) {
				C1[i] = P1[i]; C2[i] = P2[i];
			}
			else {
				C1[i] = P2[i]; C2[i] = P1[i];
			}
		}
		break;
	case eUniform: // swap odd/even genes
		for (i = 0; i<cIndividualLength; i++) {
			if (i % 2) {
				C1[i] = P1[i]; C2[i] = P2[i];
			}
			else {
				C1[i] = P2[i]; C2[i] = P1[i];
			}
		}
		break;
	case eOnePoint:  // perform 1 point x-over
		Left = RandInt(cIndividualLength);
		if (cDebug) {
			printf("Cut points: 0 <= %d <= %d\n", Left, cIndividualLength - 1);
		}
		for (i = 0; i <= Left; i++) {
			C1[i] = P1[i]; C2[i] = P2[i];
		}
		for (i = Left + 1; i<cIndividualLength; i++) {
			C1[i] = P2[i]; C2[i] = P1[i];
		}
		break;
	case eTwoPoint:  // perform 2 point x-over
		Left = RandInt(cIndividualLength - 1);
		Right = Left + 1 + RandInt(cIndividualLength - Left - 1);
		if (cDebug) {
			printf("Cut points: 0 <= %d < %d <= %d\n", Left, Right, cIndividualLength - 1);
		}
		for (i = 0; i <= Left; i++) {
			C1[i] = P1[i]; C2[i] = P2[i];
		}
		for (i = Left + 1; i <= Right; i++) {
			C1[i] = P2[i]; C2[i] = P1[i];
		}
		for (i = Right + 1; i<cIndividualLength; i++) {
			C1[i] = P1[i]; C2[i] = P2[i];
		}
		break;
	default:
		printf("Invalid crossover?\n");
		exit(1);
	}
	//uniq child C1
	//initiate the co arrays
	for (i = 0; i<cIndividualLength; i++)
		co[i] = 0;
	co0.clear();
	co2.clear();
	for (i = 0; i<cIndividualLength; i++) {
		co[C1[i]] += 1;
		if(co[C1[i]] == 2)
			co2.push_back(i);
	}
	for (i = 0; i<cIndividualLength; i++) {
		if (co[i] == 0)
			co0.push_back(i);
	}
	int s0 = co0.size();
	for (i = 0; i<s0; i++) {
		C1[co2[0]] = co0[0];
		co2.erase(co2.begin());
		co0.erase(co0.begin());
	}
	for (int i = 0; i < cIndividualLength; i++) {
		for (int j = i + 1; j < cIndividualLength; j++) {
			if (C1[i] == C1[j]) {
				cout << "C1" << endl;
				system("pause");
			}
		}
	}
	//uniq child C2
	//initiate the co arrays
	for (i = 0; i<cIndividualLength; i++)
		co[i] = 0;
	co0.clear();
	co2.clear();
	for (i = 0; i<cIndividualLength; i++) {
		co[C2[i]]++;
		if (co[C2[i]] == 2)
			co2.push_back(i);
	}
	for (i = 0; i<cIndividualLength; i++) {
		if (co[i] == 0)
			co0.push_back(i);
	}
	int s2 = co2.size();
	for (i = 0; i<s2; i++) {
		C2[co2[0]] = co0[0];
		co2.erase(co2.begin());
		co0.erase(co0.begin());
	}
	//Removing duplicate cities 
	for (int i = 0; i < cIndividualLength; i++) {
		for (int j = i + 1; j < cIndividualLength; j++) {
			if (C2[i] == C2[j]) {
				cout << "C2" << endl;
				system("pause");
			}
		}
	}
}

void Mutate(int *Member) {
	// We swap two  randomly selected cities in the member
	int num = (int)(cIndividualLength / 50);
	for (int i = 0; i < num; i++) {
		int Pick = RandInt(cIndividualLength);
		int Pick1 = RandInt(cIndividualLength);
		int t = Member[Pick];
		Member[Pick] = Member[Pick1];
		Member[Pick1] = t;
	}
}

void Copy(int *P1, int *P2, int *C1, int *C2) {
	for (int i = 0; i<cIndividualLength; i++) {
		C1[i] = P1[i]; C2[i] = P2[i];
	}
}

double Rand01() { // 0..1
	return(rand() / (double)(RAND_MAX));
}

int RandInt(int n) { // 0..n-1
	return int(rand() / (double(RAND_MAX) + 1) * n);
}
