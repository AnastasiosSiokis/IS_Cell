#include <iostream>
#include <vector>
#include <fstream>
#include <string>


//#include <random>
//#include <iomanip>
#include "math.h"
//#include "stdio.h"
//#include <sstream>
//#include <cstdlib>
//#include <ctime>
//#include <chrono>
//#include <array>
//#include <algorithm>

#define USE_DENS_LATTICES
/*
#ifdef USE_DENS_LATTICES
DensT
DensM
DensTM
CellKD
*/
// #ifdef USE_DENS_LATTICES
#define UseGraphInterface
//#define USE_Theta
#define USE_ThetaProb
//#define USE_Rounding
//#define USE_NeumannNeighbors
//#define Big_Lattice
//#define SMALL_lattice
#define dim 150  // 150 @ 70nm, 300 @ 35nm


#define USE_SHUFFLE


enum grid_states{nothing, border, TCR, LFA1, pMHC, ICAM1, TM, LI, Nucleation, NucleationPolymerized, Polymerization, CD45, emptyActin,
                 CD2, CD48, CD2_CD48, CD28, CD80, CD28_CD80, PD1, PDL1, PD1_PDL1, /*out_TCR, out_LFA1, out_pMHC, out_ICAM1, out_CD45, out_CD2, out_CD48,*/ NumberOfTypes}; // state of nodes, attribution to molecule types

enum agent_interactions{diffusion, complex_forces, binding_kinetics, internalization, externalization};

void TesteAll();

void TimerStart();
double TimerStop();

int newID();

class Vector2{
public:
    Vector2(void);
    Vector2(float X, float Y);
    Vector2(float sel_X, float sel_Y, float find_X, float find_Y);
    ~Vector2(void);

    float Length();

    Vector2 getNormalize();

    void Add(Vector2 vector);

    Vector2 getMultiplied(double weight);

    float DistanceTo(Vector2* vector);

    float X, Y;
}; // Vector2 ends here


class gridpoint{
public:
    void clear();
    gridpoint();
    //gridpoint(long int i, grid_states);
    gridpoint(int i, int j, grid_states, int IDinList);
    ~gridpoint();

    int coord_x; // coordinate x
    int coord_y; // coordinate x

    int ID;
    int ListIndex;

    inline void setID(int newID);
    inline int getID();

    inline void setListIndex(int newIndex);
    inline int getListIndex();

    //long int index; // index(deiktis) on the grid

    // coordinates of this point in the lattice
    //long x[2]; // because we have 2D lattice. for 3D-> x[3] etc

    //long int near_n[4]; // nearest neighbors (left, right, up, down)
    //long int diag_n[4]; // next-nearest neightbors (diagonal)

//------------- For the AGENTS ------------------------------------------------------------------

    //properties for the agents (Name and Number)
    grid_states agent; // to create the agents

    grid_states getAgent(); // get the name of the agent
    inline int getCoord_x();    // get coordinate x
    inline int getCoord_y();    // get coordinate y

    inline void setAgent(grid_states AgentType); // set the name of the agent
    inline void setCoord_x(int x); // set coordinate x
    inline void setCoord_y(int y); // set coordinate y

    bool isTracked;

    void copyAgent(gridpoint* originalAgent);

//-----------------------------------------------------------------------------------------------
    //grid_states status;

    Vector2 direction;

    std::string print();
}; // GRIDPOINT ends here.

class doubleLattice{
public:
    doubleLattice();
    std::vector< std::vector<double>* > latt; // fill the vector
    double access(int x , int y);
    void set(int x, int y, double label);
    void erase();
};


//template<int dim> // template to specify dimensions of vector. I defined earlier so that i can call "class ChooseAgent"...maybe not needed
class lattice{
public:
    lattice();
    std::vector< std::vector<gridpoint>* > latt; // fill the vector

    gridpoint* access(int x , int y);

    gridpoint* accessID(int _ID);

    std::map<gridpoint*, gridpoint*> Unmoved;

    //std::vector<gridpoint*> listOfAgents;

    bool StateChange(grid_states label, int x, int y);

    bool CreateAgent(grid_states label, int RandposX, int RandposY);

    grid_states Pick(int RandposX, int RandposY);

    bool Move(int oldX, int oldY, int newX,  int newY, bool refuseAddingList = false);

    bool Erase(grid_states label, int posX, int posY); // old one

    bool Erase(grid_states label, int pos);

    int CheckNeighbors(int x, int y);

    int Count(grid_states AgentID = NumberOfTypes);

    double RDF(grid_states fromWho, grid_states toWho, int posX, int posY);

    //int PossiblePositions(grid_states AgentID = NumberOfTypes, double d);

    double DistanceFromCenter(grid_states label, int x, int y);

    static void getDelta_Epsilon(int &delta, int &epsilon);

    static void getDelta_Epsilon_from_Angle(int &delta, int &epsilon, Vector2 &force);

    std::vector<gridpoint*>* GetListOfAgents(grid_states i); // pointer to the vector, faster than the previous

    double Force(double distance, double weight, double currentDistance);

    double Lennard_Jones(double Rmin, double PotentialWell, double currentDistance);

    void trackIt(int x, int y);

    void trackIt(int ID);

    void unTrackIt(int x, int y);

    void unTrackIt(int ID);

    lattice* TrackingLattice;

    void TrackRandom(grid_states label, int nb_to_track = 1);

    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
private:
    std::vector< std::vector <gridpoint* >* > listOfAgents; // vector of vector* of gridpoints -> this way we don't care if the size is different

    //std::vector<gridpoints*> GetListOfAgents(int i); // this gives a copy of the vector

    gridpoint* accessFromList(grid_states label, int posInList);

    void RemoveAgentFromList(grid_states label, int pos); // pos = position inside the list

    int CreateAgentInList(grid_states label, gridpoint* gp);

    void ChangeMemoryPosition(grid_states label, int pos, gridpoint* gp); // for the pointer

    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
public:
    void CleanLattice();
    void Honeycomb();
    void Square();
    void SquareandCircle();
    void CircleandSquare();
    void CircleandSquareandCircle();
    void CircleBig();
    void CircleMedium();
    void CircleSmall();
    void GridSizesOver25();
    void GridSizesOver15();
    void GridSizesOver10();
    void GridSizesOver5();
    void GridSizesOver3();
    void GridSizesOver2();
    void Horizontal();
    void Vertical();

    int nbAgents[NumberOfTypes];

    void checkListStructure();

}; // class lattice ends here

class simulation{
public:
    double chronos;
    double T; // = 0.012; // s WE SHOULD CHOOSE T=0.01225 IF WE WANT P_near_m=1
    double Dm; //= 0.10;  // mikro-meters^2 s^(-1)
    double Dc; // = 0.06;  // mikro-meters^2 s^(-1)
    //const double P_i = 0.3;   // Probability NOT π

    double R;  // mikro-meters
    double a; //= 0.07; // mikro-meters
    double a0;
    //double pi; //= 3.141592654; // π

    //const double Kaff_LI = pow(10.,7.);     // M^(-1)
    //const double Kaff_TM = 2 * pow(10.,5.); // M^(-1)

    double Wrep; //= -1.0;       // repulsion weight
    double Lrep; //= 0.4 / a;    // repulsion length
    double Watt; // = 0.02;      // attraction weight
    double Latt; //= 4.9 / a;    // attraction length

    double WrepPD1_PDL1; //= -1.0;       // repulsion weight
    double LrepPD1_PDL1; //= 0.4 / a;    // repulsion length

    double Att_CD2CD48;
    double Watt_CD2CD48;
    double Att_CD2TM;
    double Watt_CD2TM;

    double Att_TCRTM;
    double Watt_TCRTM;

    double Att_TMTM;
    double Watt_TMTM;

    double Rep_byTM_CD2CD48;
    double Wrep_byTM_CD2CD48;
    double Att_byTM_CD2CD48;
    double Watt_byTM_CD2CD48;

    double Att_CD2_CD28;
    double Watt_CD2_CD28;

    double Att_CD28_PD1;
    double Watt_CD28_PD1;
    double Att_CD2_PD1;
    double Watt_CD2_PD1;

    double Att_PD1_PD1;
    double Watt_PD1_PD1;


    double Watt_TM_PD1;
    double Att_TM_PD1;


    double Att_CD28CD80;
    double Watt_CD28CD80;
    double Rep_byTM_CD28CD80;
    double Wrep_byTM_CD28CD80;

    double CentrVecTM;
    double CentrVecLI;
    double CentrVecCD2CD58;
    double OutVecCD2CD58;
    double CentrVecCD28CD80;
    double CentrVecPD1_PDL1;

    double Watt_by_TM; // attraction weight of CD45 by TCR-pMHC
    double Att_by_TM;  // attraction of CD45 by TCR-pMHC
    double Watt_CD2_by_TM;
    double Att_CD2_by_TM;
    double Wrep_by_TM; // attraction weight of CD45 by TCR-pMHC
    double Rep_by_TM;
    double Wrep_by_LI; // repulsion weight of CD45 by LFA-1-ICAM-1
    double Rep_by_LI;  // repulsion of CD45 by LFA-1-ICAM-1
    double Wrep_by_CC;
    double Rep_by_CC;
    double AttCD45;
    double WattCD45;
    double RepCD45;
    double WrepCD45;

    int NbInternalizedAgents;
    int NbExternalizedAgents;

    int NbInternalizedTCRs;
    int NbInternalizedLFA1s;
    int NbInternalizedTMs;
    int NbInternalizedLIs;

    int nbInCD45;
    int nbInTCR;
    int nbInLFA;
    int nbInMHC;
    int nbInICAM;
    int nbInCD2;
    int nbInCD48;
    int nbOutCD45;
    int nbOutTCR;
    int nbOutLFA;
    int nbOutMHC;
    int nbOutICAM;
    int nbOutCD2;
    int nbOutCD48;

    std::ofstream cellKDValue;

    double nbTMs;
    double nbCCs;

    double N_A; //= 6.02214085774 * pow(10.,23.); // M^(-1) Avogadro number

    double L_TM; // meters-> Bond length TCR-pMHC
    double L_LI; // meters-> Bond length LFA1-ICAM1
    double Length_of_CD45;

//------------- For CD2-CD48 Interactions ------------------------------------------//
    double L_CD2_CD48;

    double Koff_CD2_CD48;
    double Kon_CD2_CD48;

    double Toff_CD2_CD48;
    double Ton_CD2_CD48;

    double Poff_CD2_CD48;
    double Pon_CD2_CD48;
//--------------------------------------------------------------------------------//

//------------- For CD28-CD80 Interactions ------------------------------------------//
    double L_CD28_CD80;

    double Koff_CD28_CD80;
    double Kon_CD28_CD80;

    double Toff_CD28_CD80;
    double Ton_CD28_CD80;

    double Poff_CD28_CD80;
    double Pon_CD28_CD80;
//--------------------------------------------------------------------------------//

    //------------- For PD1-PDL1 Interactions ------------------------------------------//
        double L_PD1_PDL1;

        double Koff_PD1_PDL1;
        double Kon_PD1_PDL1;

        double Toff_PD1_PDL1;
        double Ton_PD1_PDL1;

        double Poff_PD1_PDL1;
        double Pon_PD1_PDL1;
    //--------------------------------------------------------------------------------//

    // On- and Off-Probabilities calculation
    double Koff_TM; //= 0.1;   // s^(-1)
    double Koff_LI; //= 0.03;  // s^(-1)

    double Kon_TM; //= 2. * pow(10.,4.) * pow(10.,15.); // M^(-1) s^(-1)
    double Kon_LI; //= 3. * pow(10.,5.) * pow(10.,15.); // M^(-1) s^(-1)

    double Toff_TM; //= 1 / Koff_TM;
    double Toff_LI; //= 1 / Koff_LI;
    double Ton_TM; //= (L_TM * pow(a,2) * N_A) / Kon_TM;
    double Ton_LI; //= (L_LI * pow(a,2) * N_A) / Kon_LI;

    double Poff_TM; //= T/Toff_TM; // probability to unbind TM
    double Poff_LI; //= T/Toff_LI; // probability to unbind LI
    double Pon_TM; //= T/Ton_TM;  // probability to bind TM
    double Pon_LI; //= T/Ton_LI;  // probability to bind LI

    double Tinternalize; // Internalization
    double Pinternalize; // Probability of Internalization

    double TextTCR;  // EXTERNALIZATION of TCRs
    double TextLFA1; // EXTERNALIZATION of LFA1s
    double TextCD45;
    double TextCD2;
    double TextCD28;
    double TextMHC;
    double TextICAM;

    double PextTCR;  // PROBABILITY OF EXTERNALIZATION OF TCRs
    double PextLFA1; // PROBABILITY OF EXTERNALIZATION OF LFA1s
    double PextCD45;
    double PextCD2;
    double PextCD28;
    double PextMHC;
    double PextICAM;

    double P_outFluxCD45;
    double P_outFluxTCR;
    double P_outFluxLFA;
    double P_outFluxMHC;
    double P_outFluxICAM;
    double P_outFluxCD2;
    double P_outFluxCD48;

    double P_inFluxCD45;
    double P_inFluxTCR;
    double P_inFluxLFA;
    double P_inFluxMHC;
    double P_inFluxICAM;
    double P_inFluxCD2;
    double P_inFluxCD48;

    double RangeInternalization;
    double GridDistanceInternalization;


    double P_near_m; //= (4*T*Dm*sqrt(2)) / ((pow(a,2))*(sqrt(2)+1)); // molecule to nearest neighbor
    double P_near_c; //= (4*T*Dc*sqrt(2)) / ((pow(a,2))*(sqrt(2)+1)); // complex to nearest neighbor
    double P_diag_m; //= P_near_m / (sqrt(2));      // molecule to nearest neighbor
    double P_diag_c; //= P_near_c / (sqrt(2));      // complex to nearest neaighbor, basically it is LI complex

    int init_TM;
    int init_LI;

    int init_LFA1;
    int init_TCR;
    int init_pMHC;
    int init_ICAM;
    int init_CD45;
    int init_CD2;
    int init_CD48;
    int init_CD28;
    int init_CD80;
    int init_PD1;
    int init_PDL1;

    int p_LFA1;
    int p_TCR;
    int p_pMHC;
    int p_ICAM;
    int p_CD45;
    int p_CD2;
    int p_CD48;

    double Pnucleate;
    double Pdenucleate;
    double Ppolymerize;
    double Pdepolymerize;
    double Pnucl_polym; // nucleation point that gets polymerized
    double Pdenucl_polym; // nucleation point that gets de-polymerized
    double RadiusTM;
    double radiusLI;

    double Area;
    double AreaCell;
    double AreaWhole;
    double SynapseRadius;

    double NbInD_5, NbInD_10, NbInD_15, NbInD_20, NbInD_25, NbInD_30, NbInD_35, NbInD_40, NbInD_45, NbInD_50, NbInD_55, NbInD_60, NbInD_65, NbInD_70, NbInD_75;

    double PossiblePositionsTotal;

    lattice APC;   // create the lattice for the APC
    lattice Tcell; // create the lattice for the T cell
//    lattice Actin;
    lattice Colocalization;
    lattice Tracker;

    #ifdef USE_DENS_LATTICES
    void calculateDensities(double densityDistance);
    void calculateClusters();
//    doubleLattice DensT;
//    doubleLattice DensM;
//    doubleLattice DensTM;
//    doubleLattice CellKD;

    doubleLattice clusterComponents;
    #endif

    void Initialize(bool);

    void Simulate(double);
    simulation();
};

class observer{
    // Number of different MHCs
    // Average distance from initial pool
    // Average sequence distance within population
public:
    observer(std::string _filename): filename(_filename) {}
    ~observer(){}
    std::string filename;
    std::vector<double> variable;
    std::vector<double> time;
    void writeToFile(){
        std::ofstream outFile(filename.c_str(),std::ios::out);
        int row = variable.size();
        if (outFile){
            for (int i = 0; i<row; ++i)
                outFile << time[i] << "\t" << variable[i] << std::endl;
            outFile.close();
        }
        else
            std::cout << "Observer: the file cannot be created." << filename << std::endl;
    }
    void pushData(double value, double _time){
        variable.push_back(value);
        time.push_back(_time);
    }
    void clear(){
        variable.clear();
        time.clear();
    }
};

class observerTable{
public:
    observerTable(std::string _filename, int _nbCols): filename(_filename), nbCols(_nbCols) {if(nbCols > 10) std::cerr << "observerTable is not programmed for more than 10 columns/n";}
    ~observerTable(){}
    std::string filename;
    int nbCols;


    std::vector<std::string> colNames;
    std::vector< std::vector<double> > variable;

    void writeToFile(){
        int row = variable.size();
        if(row == 0) return;
        std::ofstream outFile(filename.c_str(),std::ios::out);
        if (outFile){
            for(int j = 0; j < nbCols; ++j){
                outFile << (j == 0 ? "" : "\t") << colNames[j] << std::endl;
            }
            outFile << "/n";
            for (int i = 0; i<row; ++i){
                if((int) variable[i].size() != nbCols) std::cerr << "AAARGH WTF the size of the observerTable is not square /n";
                for(int j = 0; j < nbCols; ++j){
                    outFile << (j == 0 ? "" : "\t") << variable[i][j] << std::endl;
                }
                outFile << "/n";
            }
            outFile.close();
        }
        else
            std::cout << "Observer: the file cannot be created." << filename << std::endl;
    }
    void pushData(double v1 = NAN, double v2 = NAN, double v3 = NAN, double v4 = NAN, double v5 = NAN, double v6 = NAN, double v7 = NAN, double v8 = NAN, double v9 = NAN, double v10 = NAN){
        std::vector<double> newLine = std::vector<double>();
        newLine.resize(nbCols);
        if(nbCols <= 1) newLine[0] = v1;
        if(nbCols <= 2) newLine[1] = v2;
        if(nbCols <= 3) newLine[2] = v3;
        if(nbCols <= 4) newLine[3] = v4;
        if(nbCols <= 5) newLine[4] = v5;
        if(nbCols <= 6) newLine[5] = v6;
        if(nbCols <= 7) newLine[6] = v7;
        if(nbCols <= 8) newLine[7] = v8;
        if(nbCols <= 9) newLine[8] = v9;
        if(nbCols <= 10) newLine[9] = v10;
        variable.push_back(newLine);
    }
    void giveColNames(std::string v1 = std::string(""), std::string v2 = std::string(""), std::string v3 = std::string(""), std::string v4 = std::string(""), std::string v5 = std::string(""), std::string v6 = std::string(""), std::string v7 = std::string(""), std::string v8 = std::string(""), std::string v9 = std::string(""), std::string v10 = std::string("")){
        if(nbCols <= 1) colNames[0] = v1;
        if(nbCols <= 2) colNames[1] = v2;
        if(nbCols <= 3) colNames[2] = v3;
        if(nbCols <= 4) colNames[3] = v4;
        if(nbCols <= 5) colNames[4] = v5;
        if(nbCols <= 6) colNames[5] = v6;
        if(nbCols <= 7) colNames[6] = v7;
        if(nbCols <= 8) colNames[7] = v8;
        if(nbCols <= 9) colNames[8] = v9;
        if(nbCols <= 10) colNames[9] = v10;
    }

    void clear(){
        variable.clear();
    }
};

class masterObserver{
public:
    masterObserver(std::string folder):
        TCRnumber(folder + std::string("TCR_number.txt")),
        LFA1number(folder + std::string("LFA1_number.txt")),
        CD45number(folder + std::string("CD45_number.txt")),
        pMHCnumber(folder + std::string("pMHC_number.txt")),
        ICAM1number(folder + std::string("ICAM1_number.txt")),
        TCR_pMHCnumber(folder + std::string("TCR_pMHC_number.txt")),
        LFA1_ICAM1number(folder + std::string("LFA1_ICAM1_number.txt")){}
        //agentMovement(folder + std::string("agent_movement.txt")),
        //moleculeConcentration(folder + std::string("agent_concentration.txt")){ }

    observer TCRnumber;
    observer LFA1number;
    observer CD45number;
    observer pMHCnumber;
    observer ICAM1number;
    observer TCR_pMHCnumber;
    observer LFA1_ICAM1number;
    //observer agentMovement;
    //observer moleculeConcentration;

    void clear(){
        TCRnumber.clear();
        LFA1number.clear();
        CD45number.clear();
        pMHCnumber.clear();
        ICAM1number.clear();
        TCR_pMHCnumber.clear();
        LFA1_ICAM1number.clear();
        //agentMovement.clear();
        //moleculeConcentration.clear();
    }

    void writeFiles(){
        TCRnumber.writeToFile();
        LFA1number.writeToFile();
        CD45number.writeToFile();
        pMHCnumber.writeToFile();
        ICAM1number.writeToFile();
        TCR_pMHCnumber.writeToFile();
        LFA1_ICAM1number.writeToFile();
        //agentMovement.writeToFile();
        //moleculeConcentration.writeToFile();
    }
};


int main_2(lattice &Tcell, lattice &APC, lattice &Actin);

