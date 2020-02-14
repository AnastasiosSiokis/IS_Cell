#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <iomanip>
#include "math.h"
#include "stdio.h"
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <array>
#include <algorithm>
#include <list>
#include <utility>
#include <complex>
#include <tuple>
#include <map>
#include <set>

#include "IS_cell.h"

#define pi 3.141592654
#define dr 0.21 // 2 gridpoints
std::chrono::time_point<std::chrono::system_clock> start, stop;
void TimerStart(){
    start = std::chrono::system_clock::now();
}
double TimerStop(){
    stop = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = stop-start;
    //std::cerr << "Elapsed time " << elapsed_seconds.count() << std::endl;
    return elapsed_seconds.count();
}

std::random_device *rd; // allocate a pointer
std::mt19937 *gen;
std::uniform_int_distribution<> *SpreadDistrib;          // integer distribution
std::poisson_distribution<> *SpreadDistrib_p;            // integer distribution
std::uniform_int_distribution<> *DirectionDistrib;       // integer distribution
std::uniform_int_distribution<> *DirDistrib;             // integer distribution
std::uniform_int_distribution<> *InteractionDistrib;     // integer distribution
std::uniform_int_distribution<> *DirectionDistrib2;      // integer distribution
std::uniform_real_distribution<> *RealDistrib;           // real distribution
std::normal_distribution<> *NormalDistribTCRpMHC;               // normal distribution
std::normal_distribution<> *NormalDistribLFAICAM;               // normal distribution
std::normal_distribution<> *NormalDistribCD2CD58;               // normal distribution
//std::uniform_real_distribution<> *RealDirection;           // real distribution
std::uniform_real_distribution<> *InternalizationDistrib;// real distribution

int RandGen(){return (*SpreadDistrib)(*gen);} // function for easier usage of the random engine
int RandGen_p(){return (*SpreadDistrib_p)(*gen);}
int RandDirection(){return (*DirectionDistrib)(*gen);}
int RandDir(){return (*DirDistrib)(*gen);}
int RandInteraction(){return (*InteractionDistrib)(*gen);}
double RandReal(){return (*RealDistrib)(*gen);}
double NormalDistTCRpMHC(){return (*NormalDistribTCRpMHC)(*gen);}
double NormalDistLFAICAM(){return (*NormalDistribLFAICAM)(*gen);}
double NormalDistCD2CD58(){return (*NormalDistribCD2CD58)(*gen);}
//double RandRealDirection(){return (*RealDirection)(*gen);}
double RandInterArea(){return (*InternalizationDistrib)(*gen);}


    gridpoint::gridpoint(){
       //index=0; // nothing --> apo to enum grid_states{nothing,....};
       //int i;
       //for (i=0; i<2; i++) x[i]=0;
       //for (i=0; i<4; i++) near_n[i]=0;
       //for (i=0; i<4; i++) diag_n[i]=0;
        coord_x = 0;
        coord_y = 0;
        agent = nothing;
        ListIndex = -1;
        ID = -1;
        isTracked = false;
    }
    /*gridpoint::gridpoint(long int i, grid_states s){
       index=i; // 0: nothing, 1: TCR, 2: LFA-1, 3: pMHC, 4: ICAM-1, 5: TCR-pMHC, 6: LFA-1-ICAM-1
       //int j;
       //for (j=0; j<2; j++) x[j]=0;
       //for (j=0; j<4; j++) near_n[j]=0;
       //for (j=0; j<4; j++) diag_n[j]=0;
       agent=s;
    }*/
    gridpoint::gridpoint(int i, int j, grid_states s, int IDinList){
        coord_x = i;
        coord_y = j;
        agent = s;
        ListIndex = IDinList;
        ID = -1;
        isTracked = false;
    }

    void gridpoint::copyAgent(gridpoint* originalAgent){

        coord_x = originalAgent->coord_x;
        coord_y = originalAgent->coord_y;
        agent = originalAgent->getAgent();
        ID = originalAgent->ID;
        isTracked = false;
    }

    gridpoint::~gridpoint(){}

    inline void gridpoint::setID(int newID){ID = newID;}
    inline int gridpoint::getID(){return ID;}
    inline void gridpoint::setListIndex(int newIndex){ListIndex = newIndex;}
    inline int gridpoint::getListIndex(){return ListIndex;}

    grid_states gridpoint::getAgent(){return agent;} // get the name of the agent
    inline void gridpoint::setAgent(grid_states AgentType){agent = AgentType;} // set the name of the agent

    inline int gridpoint::getCoord_x(){return coord_x;} // get coordinate x
    inline int gridpoint::getCoord_y(){return coord_y;} // get coordinate y
    inline void gridpoint::setCoord_x(int x){coord_x = x;} // set coordinate x
    inline void gridpoint::setCoord_y(int y){coord_y = y;} // set coordinate y
    std::string gridpoint::print(){
        std::stringstream res;
        res << "type " << agent << " pos " << coord_x << "," << coord_y << " ID " << ID << " listIndex " << ListIndex;
        return res.str(); // return the string stream this way
    }

    /*char operator==(const gridpoint& a, const gridpoint& b) {
      if ( (a.index!=b.index) || (a.status!=b.status) ) return 0;
      return 1;
    }*/

    /*char operator!=(const gridpoint& a, const gridpoint& b) {
      return !((a==b)== true);
    }*/


    doubleLattice::doubleLattice(){
        latt.resize(dim,NULL);
        for(int i=0;i<dim;i++){
            std::vector<double> *pointer =  new std::vector<double> (dim, 0);
            latt[i] = pointer;
        }
    }

        void doubleLattice::erase(){
            for(int i=0;i<dim;i++){
                latt[i]->clear();
                latt[i]->resize(dim,0);
            }
        }

    double doubleLattice::access(int x , int y){
        if((x < 0) || (x >= dim)) {std::cerr << "Drunk Hafiz brain overflow" << std::endl; }
        if((y < 0) || (y >= dim)) {std::cerr << "You just got tired of wrong indexes" << std::endl; }
        return (*latt[x])[y];
    }
    void doubleLattice::set(int x, int y, double label){
        if((x < 0) || (x >= dim)) {std::cerr << "Hafiz brain overflow" << std::endl; return;}
        if((y < 0) || (y >= dim)) {
            std::cerr << "You just got sick of wrong indexes" << std::endl; return;
        }
        (*latt[x])[y] = label;
    }


    lattice::lattice(){
        latt.resize(dim,NULL); //fixes the length of latt
        for(int i=0;i<dim;i++){
            std::vector<gridpoint> *pointer =  new std::vector<gridpoint> (dim);
            latt[i] = pointer;
        }
        for(int i=0;i<dim;i++){
            for(int j=0;j<dim;j++){
                access(i,j)->setCoord_x(i);
                access(i,j)->setCoord_y(j);
            }
        }

        for (int i=0; i<NumberOfTypes; ++i){
            nbAgents[i] = 0;
        }
        listOfAgents.resize(NumberOfTypes, NULL); // initialize with zeros // Structure
        for (int i=0; i<NumberOfTypes; ++i){
            listOfAgents[i] = new std::vector< gridpoint*> (); // fill all the lists
        }
        TrackingLattice = NULL;
    }

    gridpoint* lattice::access(int x , int y){ // access a certain position of a gridpoint
        return &(*latt[x])[y];
    }

    gridpoint* lattice::accessID(int _ID){ //
        for (int i=0; i<dim;++i){
            for (int j=0; j<dim;++j){
                if(_ID == this->access(i,j)->getID()){
                    return this->access(i,j);
                }
            }
        }
        return NULL;
    }

    gridpoint* lattice::accessFromList(grid_states label, int pos){
        int SizeOfList = (*listOfAgents[label]).size();
        if ((pos<0) || (pos > SizeOfList-1)){
            std::cerr << "lattice::accessFromList: NOOO you are outside the list state = " << label << " index = " << pos << " sizeList = " <<  (*listOfAgents[label]).size() << std::endl;
        }
        else{
            return (*listOfAgents[label])[pos];
        }
        return NULL;
    }

    bool lattice::CreateAgent(grid_states label, int x, int y){ // create agents with a label, random x and y, specific number
        gridpoint* position = this->access(x,y);
        //if((x == 0) && (label != border)) std::cerr << "Corrupting Borders " << std::endl;
        if (position->getAgent() == nothing){
            position->setAgent(label);
            position->setID(newID());
            ++nbAgents[position->getAgent()];
            int Index = CreateAgentInList(label, position);
            position->setListIndex(Index);
            return true;
        }
        //listOfAgents.push_back(position);
        return false;
    }

    bool lattice::Erase(grid_states label, int x, int y){ // delete the identity of the gridpoint and change it to nothing, OLD FUNCTION
        /*gridpoint* position = this->access(x,y);
        if (position->getAgent() != label){
            std::cerr<< "Erased wrong agent" <<std::endl;
            return false;
        }
        --nbAgents[position->getAgent()];
        position->setAgent(nothing);
        RemoveAgentFromList(label, position->getListIndex());
        position->setID(-1);
        position->setListIndex(-1);
        std::cerr<< "Don't use this old Erase function"<<std::endl;
        return true;*/
        gridpoint* position = this->access(x,y);
        int Index = position->getListIndex();
        return Erase(label, Index);
    }

    bool lattice::Erase(grid_states label, int pos){ // change the grid_state of the gridpoint to nothing
        gridpoint* position = this->accessFromList(label, pos);
        if (position->getAgent() != label){
            std::cerr<< "Erased wrong agent" <<std::endl;
            return false;
        }
        --nbAgents[position->getAgent()];
        position->setAgent(nothing);
        RemoveAgentFromList(label, pos);
        position->setID(-1);
        position->setListIndex(-1);
        return true;
    }

    /*double lattice::RDF(grid_states fromWho, grid_states toWho, int posX, int posY){
        this->access(posX, posY)->getAgent();

        //double distance = sqrt(pow(position->getCoord_x() -dim/2 , 2.0)+pow(position->getCoord_y()- dim/2., 2.0));


    }*/

    grid_states lattice::Pick(int RandPosX, int RandPosY){ // randomly pick an agent //grid_states Randlabel,
        return this->access(RandPosX, RandPosY)->getAgent(); // this->
    }

    bool lattice::Move(int oldX, int oldY, int newX,  int newY, bool refuseAddingList){  // N, W , S, E, NE, NW, SW, SE
        gridpoint* position = this->access(newX,newY);
        gridpoint* old_position = this->access(oldX,oldY);
        if (position->getAgent() == nothing){
            position->setAgent(old_position->getAgent());
            position->setID(old_position->getID());
            position->setListIndex(old_position->getListIndex());
            ++nbAgents[position->getAgent()];
            if ((old_position->isTracked) && (TrackingLattice)){
                TrackingLattice->access(newX,newY)->copyAgent(old_position);
            }
            position->isTracked = old_position->isTracked;
        }
        //------------- Outflux ----------------------------
        /*else if (position->getAgent() == border){
            ++nbOutMol[old_position->getAgent()];
            --nbInMol[old_position->getAgent()];
            old_position->setAgent(nothing);
        }*/
        //--------------------------------------------------
        else{
            if ((position->getAgent() != border) && (position->getAgent() != nothing)){
                //UnMoved.push_back(std::pair<gridpoint*, gridpoint*>(old_position, position));
                if(refuseAddingList) std::cerr << "Big Shit happens" << std::endl;

                if(!refuseAddingList){
                    if(Unmoved.find(old_position) != Unmoved.end()) std::cerr << "Adding two times an element" << std::endl;
                    Unmoved.insert(std::pair<gridpoint*, gridpoint*>(old_position, position));
                }
            }
            return false; // i can save the unmoved agents here....
        }

        if (access(oldX, oldY)->getAgent() == nothing){
            std::cerr << "error---can not delete emty gridpoint"<<std::endl;
            return false;
        }
        --nbAgents[old_position->getAgent()];
        old_position->setAgent(nothing);
        old_position->setID(-1);
        old_position->setListIndex(-1);
        old_position->isTracked = false;

        ChangeMemoryPosition(position->getAgent(), position->getListIndex(), position);
        return true;
    }

    /*bool lattice::Move(int oldX, int oldY, int newX,  int newY){  // N, W , S, E, NE, NW, SW, SE
        if (access(newX, newY)->getAgent() == nothing){
            CreateAgent(access(oldX, oldY)->getAgent(), newX, newY);
            if (access(oldX, oldY)->getAgent()==nothing)
                std::cerr<<"<error---can not delete emty gridpoint"<<std::endl;
            Erase(access(oldX, oldY)->getAgent(), oldX, oldY);
            return true;
        }
        return false;
    }*/

    bool lattice::StateChange(grid_states label, int x, int y){
        gridpoint* position = this->access(x,y);
        if ((position->getAgent() == nothing) || (position->getAgent() == border) ){
            std::cerr<< "Err: lattice::StateChange, Problem..."<<std::endl;
            return false;
        }
        --nbAgents[position->getAgent()]; // reduce old label
        RemoveAgentFromList(position->getAgent(),position->getListIndex());
        position->setAgent(label);
        int Index = CreateAgentInList(label, position);
        position->setListIndex(Index);
        ++nbAgents[position->getAgent()]; // increase new label
        if ((position->isTracked) && (TrackingLattice)){
            TrackingLattice->access(x, y)->copyAgent(position);
        }
        return true;
    }

    void lattice::TrackRandom(grid_states label, int nb_to_track){
        if(nb_to_track <= 0) return;

        //nbAgents[label];
        std::uniform_int_distribution<> dis(0,nbAgents[label]-1);
        int chosen = dis(*gen);
        gridpoint* spy = this->accessFromList(label, chosen);
        if(spy->isTracked) nb_to_track++;
        else spy->isTracked = true;

        TrackRandom(label, --nb_to_track);

    }

    void lattice::trackIt(int x, int y){
        this->access(x,y)->isTracked=true;
    }

    void lattice::trackIt(int ID){
        gridpoint* position = this->accessID(ID);
        if(position){
            position->isTracked=true;
        }
    }

    void lattice::unTrackIt(int x, int y){
        this->access(x,y)->isTracked=false;
    }

    void lattice::unTrackIt(int ID){
        gridpoint* position = this->accessID(ID);
        if(position){
            this->accessID(ID)->isTracked=false;
        }
    }

    void lattice::getDelta_Epsilon(int &delta, int &epsilon){
        //int newX, newY;
        double dir = RandReal();
        if (dir <= sqrt(2)/(sqrt(2)+1)){
            double v_or_h = RandReal(); // vertical or horizontal
            if (v_or_h<0.5){
                delta = 2*RandDir()-1;
                epsilon = 0;
            }
            else {
                delta = 0;
                epsilon = 2*RandDir()-1;
            }
        }
        else {
            delta = 2 * RandDir() - 1;
            epsilon = 2 * RandDir() - 1;
        }
        /*if ((delta != 0) && (epsilon != 0)){
            std::cerr<<"Diag"<<std::endl;
        }
        else
            std::cerr<<"Near"<<std::endl;*/
    }

    double lengthFromTheta(double theta){
        return sqrt(1.+tan(theta)*tan(theta));
    }

    void lattice::getDelta_Epsilon_from_Angle(int &delta, int &epsilon, Vector2 &force){
        float theta;
        float PI = 3.141592654;
        delta = 0;
        epsilon = 0;
        if((std::fabs(force.X) < 1e-12) && (std::fabs(force.Y) < 1e-12)) {
            delta = 0;
            epsilon = 0;
            return;
        }
        //int newX, newY;
        theta = (atan2(force.Y, force.X));
        //theta += (RandReal() -0.5) * 50.0;
        //std::cerr<<theta<<std::endl;

#ifdef USE_Theta
        theta = (atan2(force.Y, force.X)*180./PI);
        //theta += (RandReal() -0.5) * 50.0;
        if ((theta > -31.8198) && (theta < 31.8198)){
            delta = 1;
            epsilon = 0;
        }
        else if ((theta >= 31.8198) && (theta <= 58.1802)){
            delta = 1;
            epsilon = 1;
        }
        else if ((theta > 58.1802) && (theta < 121.8198)){
            delta = 0;
            epsilon = 1;
        }
        else if ((theta >= 121.8198) && (theta <= 148.1802)){
            delta = -1;
            epsilon = 1;
        }
        else if ((theta > 148.1802) /*&& (theta <= 180.)*/){
            delta = -1;
            epsilon = 0;
        }
        else if (/*(theta >= -180.) && */(theta < 211.8198-360.0)){
            delta = -1;
            epsilon = 0;
        }
        else if ((theta >= 211.8198-360.0) && (theta <= 238.1802-360.0)){
            delta = -1;
            epsilon = -1;
        }
        else if ((theta > 238.1802-360.0) && (theta < 301.8198-360.0)){
            delta = 0;
            epsilon = -1;
        }
        else if ((theta >= 301.8198-360.0) && (theta <= 328.1802-360.0)){
            delta = 1;
            epsilon = -1;
        }
        else{
        std::cerr<<"No value taken"<<std::endl;
        std::cerr<<"Theta = " << theta <<" \t d = "<< delta <<"\t e ="<< epsilon<< "\t SumX"<< force.X << "\t SumY"<< force.Y <<std::endl;
        }
        //std::cerr<<theta<<std::endl;
#endif

#ifdef USE_ThetaProb
        double length;
        if ( (theta >= -PI) && (theta <= -3.*PI/4.) ){
            length = lengthFromTheta(theta+PI);
            if (RandReal() <= 1./length){
                delta = -1.;
            }
            if (RandReal() < tan(theta)/length){
                epsilon = -1.;
            }
        }
        else if ( (theta > -3.*PI/4.) && (theta <= -PI/2.) ){
            length = lengthFromTheta(-PI/2-theta);
            if (RandReal() < -tan(theta+PI/2.)/length){
                delta = -1.;
            }
            if (RandReal() <= 1./length){
                epsilon = -1.;
            }
        }
        else if ( (theta > -PI/2.) && (theta <= -PI/4.) ){
            length = lengthFromTheta(PI/2+theta);
            if (RandReal() < tan(theta+PI/2.)/length){
                delta = 1.;
            }
            if (RandReal() <= 1./length){
                epsilon = -1.;
            }
        }
        else if ( (theta > -PI/4.) && (theta <= 0.0) ){
            length = lengthFromTheta(-theta);
            if (RandReal() <= 1./length){
                delta = 1.;
            }
            if (RandReal() < -tan(theta)/length){
                epsilon = -1;
            }
        }
        else if ( (theta > 0.0) && (theta <= PI/4.) ){
            length = lengthFromTheta(theta);
            if (RandReal() <= 1./length){
                delta = 1.;
            }
            if (RandReal() < tan(theta)/length){
                epsilon = 1.;
            }
        }
        else if ( (theta > PI/4.) && (theta <= PI/2.) ){
            length = lengthFromTheta(PI/2-theta);
            if (RandReal() < -tan(theta+PI/2.)/length){
                delta = 1.;
            }
            if (RandReal() <= 1./length){
                epsilon = 1.;
            }
        }
        else if ( (theta > PI/2.) && (theta <= 3*PI/4.) ){
            length = lengthFromTheta(-PI/2+theta);
            if (RandReal() < tan(theta+PI/2.)/length){
                delta = -1.;
            }
            if (RandReal() <= 1./length){
                epsilon = 1.;
            }
        }
        else if ( (theta > 3*PI/4) && (theta <= PI) ){
            length = lengthFromTheta(PI-theta);
            if (RandReal() <= 1./length){
                delta = -1.;
            }
            if (RandReal() < -tan(theta)/length){
                epsilon = 1.;
            }
        }
        //std::cerr << "length = " << length << std::endl;
#endif

#ifdef USE_Rounding
        delta = std::round(force.X);
        epsilon = std::round(force.Y);
#endif

#ifdef USE_NeumannNeighbors
        theta = (atan2(force.Y, force.X)*180./PI);
        if ((theta > -22.5) && (theta < 22.5)){
            delta = 1;
            epsilon = 0;
        }
        else if ((theta >= 22.5) && (theta <= 67.5)){
            delta = 1;
            epsilon = 1;
        }
        else if ((theta > 67.5) && (theta < 112.5)){
            delta = 0;
            epsilon = 1;
        }
        else if ((theta >= 112.5) && (theta <= 157.5)){
            delta = -1;
            epsilon = 1;
        }
        else if ((theta > 157.5) /*&& (theta <= 180.)*/){
            delta = -1;
            epsilon = 0;
        }
        else if (/*(theta >= -180.) && */(theta < -157.5)){
            delta = -1;
            epsilon = 0;
        }
        else if ((theta >= -157.5) && (theta <= -112.5)){
            delta = -1;
            epsilon = -1;
        }
        else if ((theta > -112.5) && (theta < -67.5)){
            delta = 0;
            epsilon = -1;
        }
        else if ((theta > -67.5) && (theta < -22.5)){
            delta = 1;
            epsilon = -1;
        }
        else{
        std::cerr<<"No value taken"<<std::endl;
        std::cerr<<"Theta = " << theta <<" \t d = "<< delta <<"\t e ="<< epsilon<< "\t SumX"<< force.X << "\t SumY"<< force.Y <<std::endl;
        }

        if (delta*epsilon != 0){
            if ( RandReal() > (1/sqrt(2)) ){
                delta = 0;
                epsilon =0;
            }
        }
#endif
    }

    int lattice::CheckNeighbors(int x, int y){ // check the 8 neighbors around the picked gridpoint
        int NumOccupiedNeighbors = 0;
        if ((x>0) && (x<dim-1) && (y>0) && (y<dim-1)){
            if(this->access(x,y)->getAgent()!=nothing){
                for (int i=-1; i<=1;++i){
                    for (int j=-1; j<=1;++j){
                        if( !((i==0) && (j==0)) && (this->access(x+i,y+j)->getAgent() != nothing)){
                            NumOccupiedNeighbors++;
                        }
                    }
                }
            }
        }
        else if ((x<=0) && (x>=dim-1) && (y<=0) && (y>=dim-1)){
            std::cerr<< "Danger Danger!!"<< std::endl;
        }
        return NumOccupiedNeighbors;
    }

    int lattice::Count(grid_states AgentID){
        if (AgentID != NumberOfTypes){
            return nbAgents[AgentID];
        }
        int nb = 0;
        for (int i=0; i< NumberOfTypes; ++i){
            nb += nbAgents[i];
        }
        return nb-nbAgents[nothing]-nbAgents[border];
        /*//int AgentNumber = 0;
        for (int i=0; i<dim; ++i){
            for (int j=0; j<dim; ++j){
                if ((NumberOfTypes != AgentID) && (this->access(i, j)->getAgent() == AgentID)){ // get the number of specific agent
                    ++AgentNumber;
                }
                if ((NumberOfTypes == AgentID) && (this->access(i, j)->getAgent() != nothing) && (this->access(i, j)->getAgent() != border)){
                    ++AgentNumber; // get the number of ALL the agents
                }
            }
        }
        //return AgentNumber;*/
    }

    /*int lattice::PossiblePositions(grid_states AgentID, int x, int y, double d){
        double xc = dim/2;
        double yc = dim/2;
        gridpoint* position = this->access(x,y);
        d = sqrt(pow(xc- (float)i, 2.0)+pow(yc- (float)j, 2.0));
        if (AgentID != NumberOfTypes){
            return nbAgents[AgentID];
        }
        int nb = 0;
        for (int i=0; i< NumberOfTypes; ++i){
            nb += nbAgents[i];
        }
        return nb-nbAgents[border];
    }*/

    /*double lattice::DistanceFromCenter(grid_states label, int x, int y){
        for(int i=0; i<NumberOfTypes; i++){
            gridpoint* position = this->access(x, y);
            if (position->getAgent() == label){
            //position->getID();
            double dist = sqrt(pow(position->getCoord_x() - dim/2., 2.0)+pow(position->getCoord_y()- dim/2., 2.0));
            return dist;
            }
            else{
                std::cerr<< "Problem with distance from the centre of " << this->access(x,y)->getAgent() << std::endl;
            }
        }
        //return 0;
    }*/

    void lattice::RemoveAgentFromList(grid_states label, int pos){
        int SizeOfList = (*listOfAgents[label]).size();
        if ((pos<0) || (pos > SizeOfList-1)){
            std::cerr << "lattice::RemoveAgentFromList: NOOO you are outside the list state = " << label << " index = " << pos << " sizeList = " <<  (*listOfAgents[label]).size() << std::endl;
        }
        else{
            (*listOfAgents[label])[pos] = (*listOfAgents[label])[SizeOfList-1]; // take the last position and store it to the position that i want to remove
            (*listOfAgents[label])[pos]->setListIndex(pos);
            (*listOfAgents[label]).resize(SizeOfList-1);
        }
    }

    int lattice::CreateAgentInList(grid_states label, gridpoint* gp){
        (*listOfAgents[label]).push_back(gp);
        return (*listOfAgents[label]).size()-1;
    }

    void lattice::ChangeMemoryPosition(grid_states label, int pos, gridpoint* gp){
        int SizeOfList = (*listOfAgents[label]).size();
        if ((pos<0) || (pos > SizeOfList-1)){
            std::cerr << "lattice::ChangeMemoryPosition: NOOO you are outside the list state = " << label << " index = " << pos << " sizeList = " <<  (*listOfAgents[label]).size() << std::endl;
        }
        else{
            (*listOfAgents[label])[pos] = gp; // new address(pointer) of the label
        }
    }

   void lattice::checkListStructure(){
        // checks the size
       for(int i = 1; i < NumberOfTypes; ++i){
           if((int) (*listOfAgents[i]).size() != nbAgents[i])
               std::cerr << "size of List " << (*listOfAgents[i]).size() << "of agents and nbAgents " << nbAgents[i] << " are not matching " << std::endl;
       }
   }

   double lattice::Force(double distance, double weight, double currentDistance){
       /*if (currentDistance <= 0.0){
           std::cerr<< "You are checking your self you stupid idiot"<<std::endl;
       }*/
       if (currentDistance>0){
           //return distance*weight/(currentDistance);
           //return std::min(distance*weight/(currentDistance),weight);
           //return std::min(weight*(std::sqrt(distance/currentDistance*currentDistance)),weight);
           //return distance*distance*weight/(currentDistance*currentDistance);
           return std::min(distance*distance*weight/(currentDistance*currentDistance),weight);
           //return weight*(std::sqrt(distance/currentDistance*currentDistance));
           //return 1-tanh(3.14 - distance);
       }
        return 0;
   }


   double lattice::Lennard_Jones(double Rmin, double PotentialWell, double currentDistance){
       /*if (currentDistance <= 0.0){
           std::cerr<< "You are checking your self you stupid idiot"<<std::endl;
       }*/
       if (currentDistance>0){
           //return 4*PotentialWell*(pow(Rmin/currentDistance,12)-pow(Rmin/currentDistance,6));
           //return std::min(4*PotentialWell*(pow(Rmin/currentDistance,12)-pow(Rmin/currentDistance,6)), PotentialWell);
           //return 24*PotentialWell*(pow(Rmin,6)/pow(currentDistance,7)-2*pow(Rmin,12)/pow(currentDistance,13));
           //return std::min(24*PotentialWell*pow(Rmin,6)/pow(currentDistance,7)-48*PotentialWell*pow(Rmin,12)/pow(currentDistance,13),PotentialWell);
           //return 48.0*PotentialWell/currentDistance*(pow(Rmin/currentDistance,12.0)-(1.0/2.0)*pow(Rmin/currentDistance,6.0));
           return std::min(48.0*PotentialWell/currentDistance*(pow(Rmin/currentDistance,12.0)-(1.0/2.0)*pow(Rmin/currentDistance,6.0)),PotentialWell);
       }
        return 0;
   }

   std::vector<gridpoint*>* lattice::GetListOfAgents(grid_states i){ // pointer to the vector, faster than the previous
        return listOfAgents[i];
   }

    void lattice::CleanLattice(){

        for (int i=0; i<dim; ++i){
            for (int j=0; j<dim; ++j){
                gridpoint* position = this->access(i, j);
                if (position->getAgent() != nothing){
                    position->setAgent(nothing);
                }
                //else
                //    std::cerr << "Already empty" << std::endl;
            }
        }
    }
/*#ifdef BIG_Lattice
    void lattice::Honeycomb(){
        for (int i=0; i<dim-1-100; ++i){ // honeycomb up-right part
            CreateAgent(border, dim/2+i, dim/2+i); // ok
            CreateAgent(border, dim/2+i, dim/2+1+i); // ok
            //CreateAgent(border, dim/2+i, 115+i); // ok
            //CreateAgent(border, dim/2+i, 116+i); // ok
            CreateAgent(border, dim/2+i, 130+i); // ok
            CreateAgent(border, dim/2+i, 131+i); // ok
            //CreateAgent(border, dim/2+i, 145+i); // ok
            //CreateAgent(border, dim/2+i, 146+i); // ok
            CreateAgent(border, dim/2+i, 160+i); // ok
            CreateAgent(border, dim/2+i, 161+i); // ok
        }
        for (int i=0; i<100; ++i){ // honeycomb up-left part
            CreateAgent(border, i, dim-i); // ok
            CreateAgent(border, i, dim+1-i); // ok
            //CreateAgent(border, i, dim+15-i); // ok
            //CreateAgent(border, i, dim+16-i); // ok
            CreateAgent(border, i, dim+30-i); // ok
            CreateAgent(border, i, dim+31-i); // ok
            //CreateAgent(border, i, dim+45-i); // ok
            //CreateAgent(border, i, dim+46-i); // ok
            CreateAgent(border, i, dim+60-i); // ok
            CreateAgent(border, i, dim+61-i); // ok
        }
        for (int i=0; i<dim-1-145; ++i){ // honeycomb down-right part
            CreateAgent(border, 115+i, 85+i); // ok
            CreateAgent(border, 115+i, 86+i); // ok
            CreateAgent(border, 130+i, 70+i); // ok
            CreateAgent(border, 130+i, 71+i); // ok
            CreateAgent(border, 145+i, 55+i); // ok
            CreateAgent(border, 145+i, 56+i); // ok
        }
        for (int i=0; i<85; ++i){ // honeycomb down-left part
            CreateAgent(border, i, dim-30-i); // ok
            CreateAgent(border, i, dim-31-i); // ok
        }
        for (int i=0; i<70; ++i){ // honeycomb down-left part
            CreateAgent(border, i, dim-60-i); // ok
            CreateAgent(border, i, dim-61-i); // ok
        }
        for (int i=0; i<55; ++i){ // honeycomb down-left part
            CreateAgent(border, i, dim-90-i); // ok
            CreateAgent(border, i, dim-91-i); // ok
        }
        for (int i=0; i<=dim/2; ++i){ // honeycomb middle vertical
            CreateAgent(border, dim/2, i); // ok
        }
        for (int i=0; i<=85; ++i){  // honeycomb vertical
            CreateAgent(border, 115, i); // ok
            CreateAgent(border, 85, i); // ok
        }
        for (int i=0; i<=70; ++i){ // honeycomb vertical
            CreateAgent(border, 130, i); // ok
            CreateAgent(border, 70, i); // ok
        }
        for (int i=0; i<=55; ++i){ // honeycomb vertical
            CreateAgent(border, 145, i); // ok
            CreateAgent(border, 55, i); // ok
        }
    }

    void lattice::Square(){
        for (int i=0; i<dim-1;++i){ // square borders
            CreateAgent(border, 2*dim/7, i);
            CreateAgent(border, 5*dim/7, i);
            CreateAgent(border, i, 2*dim/7);
            CreateAgent(border, i, 5*dim/7);
        }
    }

    void lattice::Horizontal(){
        for (int i=1; i<=dim-2; ++i){
            CreateAgent(border, i, 50);
            CreateAgent(border, i, 70);
            CreateAgent(border, i, 90);
            CreateAgent(border, i, 110);
            CreateAgent(border, i, 130);
            CreateAgent(border, i, 150);
        }
    }

    void lattice::Vertical(){
        for (int i=1; i<=dim-2; ++i){
            CreateAgent(border, 50, i);
            CreateAgent(border, 70, i);
            CreateAgent(border, 90, i);
            CreateAgent(border, 110, i);
            CreateAgent(border, 130, i);
            CreateAgent(border, 150, i);
        }
    }
#endif */

//#ifdef SMALL_Lattice
    void lattice::Honeycomb(){
        for (int i=0; i<dim-1-dim/2; ++i){ // honeycomb up-right part
            CreateAgent(border, dim/2+i, dim/2+i); // ok
            CreateAgent(border, dim/2+i, dim/2+1+i); // ok
            CreateAgent(border, dim/2+i, 97+i); // ok
            CreateAgent(border, dim/2+i, 98+i); // ok
            CreateAgent(border, dim/2+i, 120+i); // ok
            CreateAgent(border, dim/2+i, 121+i); // ok
        }
        for (int i=0; i<dim/2; ++i){ // honeycomb up-left part
            CreateAgent(border, i, dim-i); // ok
            CreateAgent(border, i, dim+1-i); // ok
            CreateAgent(border, i, dim+22-i); // ok
            CreateAgent(border, i, dim+23-i); // ok
            CreateAgent(border, i, dim+45-i); // ok
            CreateAgent(border, i, dim+46-i); // ok
        }
        for (int i=0; i<dim-1-100; ++i){ // honeycomb down-right part
            CreateAgent(border, 5*dim/8+i, 70+i); // ok
            CreateAgent(border, 5*dim/8+i, 71+i); // ok
        }
        for (int i=0; i<37; ++i){ // honeycomb down-right part
            CreateAgent(border, 6*dim/8+i, 65+i); // ok
            CreateAgent(border, 6*dim/8+i, 66+i); // ok
        }
        for (int i=0; i<19; ++i){ // honeycomb down-right part
            CreateAgent(border, 7*dim/8+i, 60+i); // ok
            CreateAgent(border, 7*dim/8+i, 61+i); // ok
        }
        for (int i=0; i<3*dim/8+1; ++i){ // honeycomb down-left part
            CreateAgent(border, i, dim-22-i); // ok
            CreateAgent(border, i, dim-23-i); // ok
        }
        for (int i=0; i<2*dim/8+1; ++i){ // honeycomb down-left part
            CreateAgent(border, i, dim-47-i); // ok
            CreateAgent(border, i, dim-48-i); // ok
        }
        for (int i=0; i<dim/8; ++i){ // honeycomb down-left part
            CreateAgent(border, i, dim-72-i); // ok
            CreateAgent(border, i, dim-73-i); // ok
        }
        for (int i=0; i<=dim/2; ++i){ // honeycomb middle vertical
            CreateAgent(border, dim/2, i); // ok
        }
        for (int i=0; i<=70; ++i){  // honeycomb vertical
            CreateAgent(border, 5*dim/8, i); // ok
            CreateAgent(border, 3*dim/8, i); // ok
        }
        for (int i=0; i<=65; ++i){ // honeycomb vertical
            CreateAgent(border, 6*dim/8, i); // ok
            CreateAgent(border, 2*dim/8, i); // ok
        }
        for (int i=0; i<=60; ++i){ // honeycomb vertical
            CreateAgent(border, 7*dim/8, i); // ok
            CreateAgent(border, dim/8, i); // ok
        }
    }

    void lattice::Square(){
        for (int i=0; i<dim;++i){ // square borders
            CreateAgent(border, 2*dim/7, i);
            CreateAgent(border, 5*dim/7, i);
            CreateAgent(border, i, 2*dim/7);
            CreateAgent(border, i, 5*dim/7);
        }
    }

    void lattice::SquareandCircle(){
        double xc = dim/2;
        double yc = dim/2;
        for (int i=0; i<dim; ++i){
            for (int j=0; j<dim; ++j){
                double distance = sqrt(pow(xc- (float)i, 2.0)+pow(yc- (float)j, 2.0));
                if ((distance >= 20) && (distance <= 21)){
                    CreateAgent(border, i, j);
                }
            }
        }
        for (int i=0; i<dim;++i){ // square borders
            CreateAgent(border, 2*dim/7, i);
            CreateAgent(border, 5*dim/7, i);
            CreateAgent(border, i, 2*dim/7);
            CreateAgent(border, i, 5*dim/7);
        }
    }

    void lattice::CircleandSquare(){
        double xc = dim/2;
        double yc = dim/2;
        for (int i=0; i<dim; ++i){
            for (int j=0; j<dim; ++j){
                double distance = sqrt(pow(xc- (float)i, 2.0)+pow(yc- (float)j, 2.0));
                if ((distance >= 50) && (distance <= 51)){
                    CreateAgent(border, i, j);
                }
            }
        }
        for (int i=0; i<dim;++i){ // square borders
            CreateAgent(border, 3*dim/7, i);
            CreateAgent(border, 4*dim/7, i);
            CreateAgent(border, i, 3*dim/7);
            CreateAgent(border, i, 4*dim/7);
        }
    }

    void lattice::CircleandSquareandCircle(){
        double xc = dim/2;
        double yc = dim/2;
        for (int i=0; i<dim; ++i){
            for (int j=0; j<dim; ++j){
                double distance = sqrt(pow(xc- (float)i, 2.0)+pow(yc- (float)j, 2.0));
                if ((distance >= 50) && (distance <= 51)){
                    CreateAgent(border, i, j);
                }
                if ((distance >= 25) && (distance <= 26)){
                    CreateAgent(border, i, j);
                }
            }
        }
        for (int i=0; i<dim;++i){ // square borders
            CreateAgent(border, 3*dim/7, i);
            CreateAgent(border, 4*dim/7, i);
            CreateAgent(border, i, 3*dim/7);
            CreateAgent(border, i, 4*dim/7);
        }
    }

    void lattice::CircleBig(){
        double xc = dim/2;
        double yc = dim/2;
        for (int i=0; i<dim; ++i){
            for (int j=0; j<dim; ++j){
                double distance = sqrt(pow(xc- (float)i, 2.0)+pow(yc- (float)j, 2.0));
                if ((distance >= 50) && (distance <= 52)){
                    CreateAgent(border, i, j);
                }
            }
        }
    }

    void lattice::CircleMedium(){
        double xc = dim/2;
        double yc = dim/2;
        for (int i=0; i<dim; ++i){
            for (int j=0; j<dim; ++j){
                double distance = sqrt(pow(xc- (float)i, 2.0)+pow(yc- (float)j, 2.0));
                if ((distance >= 30) && (distance <= 32)){
                    CreateAgent(border, i, j);
                }
            }
        }
    }

    void lattice::CircleSmall(){
        double xc = dim/2;
        double yc = dim/2;
        for (int i=0; i<dim; ++i){
            for (int j=0; j<dim; ++j){
                double distance = sqrt(pow(xc- (float)i, 2.0)+pow(yc- (float)j, 2.0));
                if ((distance >= 10) && (distance <= 11)){
                    CreateAgent(border, i, j);
                }
            }
        }
    }

    void lattice::GridSizesOver25(){ // a tiny mesh
        for (int i=0; i<dim; i+=dim/25){
            for (int j=0; j<dim-1; ++j){
                CreateAgent(border, i, j);
                CreateAgent(border, j, i);
            }
        }
    }

    void lattice::GridSizesOver15(){ // a tiny mesh
        for (int i=0; i<dim; i+=dim/15){
            for (int j=0; j<dim-1; ++j){
                CreateAgent(border, i, j);
                CreateAgent(border, j, i);
            }
        }
    }
    void lattice::GridSizesOver10(){ // a small mesh
        for (int i=0; i<dim; i+=dim/10){
            for (int j=0; j<dim-1; ++j){
                CreateAgent(border, i, j);
                CreateAgent(border, j, i);
            }
        }
    }

    void lattice::GridSizesOver5(){ // a medium mesh
        for (int i=0; i<dim; i+=dim/5){
            for (int j=0; j<dim-1; ++j){
                CreateAgent(border, i, j);
                CreateAgent(border, j, i);
            }
        }
    }

    void lattice::GridSizesOver3(){ // big square at the center
        for (int i=0; i<dim; i+=dim/3){
            for (int j=0; j<dim-1; ++j){
                CreateAgent(border, i, j);
                CreateAgent(border, j, i);
            }
        }
    }

    void lattice::GridSizesOver2(){ // cross in the middle
        for (int i=0; i<dim; i+=dim/2){
            for (int j=0; j<dim-1; ++j){
                CreateAgent(border, i, j);
                CreateAgent(border, j, i);
            }
        }
    }

    void lattice::Horizontal(){
        for (int i=1; i<=dim-1; ++i){
            CreateAgent(border, i, 1*dim/7);
            CreateAgent(border, i, 2*dim/7);
            CreateAgent(border, i, 3*dim/7);
            CreateAgent(border, i, 4*dim/7);
            CreateAgent(border, i, 5*dim/7);
            CreateAgent(border, i, 6*dim/7);
        }
    }

    void lattice::Vertical(){
        for (int i=1; i<=dim-2; ++i){
            CreateAgent(border, 1*dim/7, i);
            CreateAgent(border, 2*dim/7, i);
            CreateAgent(border, 3*dim/7, i);
            CreateAgent(border, 4*dim/7, i);
            CreateAgent(border, 5*dim/7, i);
            CreateAgent(border, 6*dim/7, i);
        }
    }
//#endif

//----------------------------------------------------------------------------------------------------------------
Vector2::Vector2(void){
    this->X = 0;
    this->Y = 0;
}

Vector2::Vector2(float X, float Y){
    this->X = X;
    this->Y = Y;
}

Vector2::Vector2(float sel_X, float sel_Y, float find_X, float find_Y){
    this->X = find_X - sel_X;
    this->Y = find_Y - sel_Y;
}

float Vector2::Length(){
    return sqrt(X * X + Y * Y);
}

Vector2 Vector2::getNormalize(){
    Vector2 vector= Vector2(0.,0.);
    float length = this->Length();
        if(length > 1e-12){
        vector.X = X/length;
        vector.Y = Y/length;
    }
    return vector;
}
void Vector2::Add(Vector2 u){
    X = u.X + X;
    Y = u.Y + Y;
}

Vector2 Vector2::getMultiplied(double weight){
    Vector2 newVector =  Vector2();
    newVector.X = weight * X;
    newVector.Y = weight * Y;
    return newVector;
}

float Vector2::DistanceTo(Vector2* v){
    float distance = sqrt(pow(v->X-X,2)+pow(v->Y-Y,2));
    return distance;
}

Vector2::~Vector2(void){}

//*************************************************************************************************************************************************

#ifndef UseGraphInterface
int main(){
    simulation* gamoto = new simulation();
    gamoto->Initialize();
    for(int i = 0; i < 30*60; ++i){
        gamoto-> Simulate(1);
    }
}
#endif

simulation::simulation() : APC(), Tcell(), Colocalization(), Tracker() {
    //std::cerr << "state constr 5,0 " << APC.access(0,5)->getAgent() << std::endl;
    Tcell.TrackingLattice = &Tracker;
    APC.TrackingLattice = &Tracker;
}


int newID(){ // create new ID every time that you create an agent
    static int ID = 0; // always remembers the previous ID
    return ++ID;
}

void simulation::Initialize(bool DefaultValues){

    //cellKDValue.open("cellKDValue.txt", std::ofstream::out | std::ofstream::app); // it is defined in .h

    //TimerStart();
//std::cerr << "state init" << APC.access(0,5)->getAgent() << std::endl;

    rd   = new std::random_device();
    gen  = new std::mt19937 ((*rd)());
    SpreadDistrib  = new std::uniform_int_distribution<> (1,dim-2);      // integer  distribution between 1 and dim-2
    SpreadDistrib_p  = new std::poisson_distribution<> (dim/2);
    DirectionDistrib = new std::uniform_int_distribution<> (-1,1);       // integer  distribution between -1 and 1
    DirDistrib = new std::uniform_int_distribution<> (0,1);              // integer  distribution between 0 and 1
    InteractionDistrib = new std::uniform_int_distribution<> (0,4);      // integer  distribution between 0 and 2
    RealDistrib = new std::uniform_real_distribution<> (0,1);            // real distribution between 0 and 1
    NormalDistribTCRpMHC = new std::normal_distribution<> (923,923*0.25);         //
    NormalDistribLFAICAM = new std::normal_distribution<> (2309,2309*0.25);         //
    NormalDistribCD2CD58 = new std::normal_distribution<> (2770,923*0.25);         //
    //InteractionDistrib = new std::uniform_real_distribution<> (0.6,1.5); // real distribution between 0.6 and 1.5

    a = 0.07;    // mikro-meters
    a0 = 0.07;   // mikro-meters
    //std::cerr<<a<<std::endl;
    R = 4.9;    //4.9;  // mikro-meters
    /*if (2*R/a >= dim){
        std::cerr<< "your synapse is bigger than the lattice"<<std::endl;
    }*/
    Area = R*R/(a*a)*pi; // area of synapse
    SynapseRadius = R/a;
    //AreaCell = 4*R*R/(a*a)*pi;
    AreaWhole = dim*dim;
    //std::cout << "Surfaces = " << Area << '\t' << AreaCell << '\t' << AreaWhole << std::endl;
    std::cout << "Synapse surface = " << Area << '\t' << "Lattice = " << AreaWhole << '\t' << "Syanpse Radius = " << SynapseRadius << std::endl;
    //Tcell.CircleMedium();
    //APC.CircleMedium();
    //APC.Honeycomb();
    //APC.Horizontal();
    //APC.Square();
    //APC.SquareandCircle();
    //APC.CircleandSquare();
    //APC.CircleandSquareandCircle();
//    APC.CircleBig();
//    APC.CircleMedium();
//    Tcell.CircleBig();
//    Tcell.CircleMedium();
    //APC.CircleSmall();
    //APC.Vertical();
    //APC.GridSizesOver25();
    //APC.GridSizesOver15();
    //APC.GridSizesOver10();
    //APC.GridSizesOver5();
    //APC.GridSizesOver3();
    //APC.GridSizesOver2(); // just a cross in the middle...no point in using it..


    //std::cerr << "state before " << APC.access(0,5)->getAgent() << std::endl;
    // for the circular contact region.
    float xc = dim/2; //x center
    float yc = dim/2; //y center

    for (int i=0; i<dim; ++i){
        for (int j=0; j<dim; ++j){
            double distance = sqrt(pow(xc- (float)i, 2.0)+pow(yc- (float)j, 2.0));
            if ((distance >= SynapseRadius) /*&& (distance <= 72)*/){
                Tcell.CreateAgent(border, i, j);
                APC.CreateAgent(border, i, j);
            }
        }
    }

    //Test for -- dr --
    /*for (int i=0; i<dim; ++i){
        for (int j=0; j<dim; ++j){
            double distance = sqrt(pow(3*dim/4- (float)i, 2.0)+pow(3*dim/4- (float)j, 2.0));
            if ((distance >= 1) && (distance <= 1+dr/a)){
                Tcell.CreateAgent(border, i, j);
                APC.CreateAgent(border, i, j);
            }
        }
    }//*/

    for ( int i=0; i<dim;++i){
        Tcell.CreateAgent(border, i, 0);
        Tcell.CreateAgent(border, i, dim-1);
        Tcell.CreateAgent(border, 0, i);
        Tcell.CreateAgent(border, dim-1, i);
        APC.CreateAgent(border, i, 0);
        APC.CreateAgent(border, i, dim-1);
        APC.CreateAgent(border, 0, i);
        APC.CreateAgent(border, dim-1, i);
    }

    /*for (int i=0; i<dim; ++i){
        for (int j=0; j<dim; ++j){
           double distan = sqrt(pow(xc- (float)i, 2.0)+pow(yc- (float)j, 2.0));
           if (distan>2.){
               Tcell.CreateAgent(border,i,j);
            }
           Tcell.CreateAgent(TCR,xc,yc);
           if (distan <= 2.)
               Tcell.CreateAgent(TCR,i,j);
        }
    }*/

    //std::cerr << "state at after " << APC.access(0,5)->getAgent() << std::endl;
    /*for (int i=0; i<=dim-1; ++i){
        Tcell.CreateAgent(border, i, 0);
        Tcell.CreateAgent(border, i, dim-1);
        APC.CreateAgent(border, i, 0);
        APC.CreateAgent(border, i, dim-1);
    }
    for (int i=0; i<dim-1; ++i){
        Tcell.CreateAgent(border, 0, i);
        Tcell.CreateAgent(border, dim-1, i);
        APC.CreateAgent(border, 0, i);
        APC.CreateAgent(border, dim-1, i);
    }*/

    //TesteAll();
    //exit(-1);
    /*for (int i = 50; i < 100; ++i){
        for (int j = 50; j <100; ++j){
            Tcell.CreateAgent(TCR, i, j);
        }
    }*/

    T  = 0.01225; // s WE SHOULD CHOOSE T=0.01225 IF WE WANT P_near_m=1
    std::cerr << "T = " << T <<std::endl;
    chronos = 0.0;
    if (DefaultValues){

#define NoDTDepDiffusion false // true = diffusing at each time step

        Dm  = 0.10;  // mikro-meters^2 s^(-1)
        Dc  = 0.06;  // mikro-meters^2 s^(-1)

        //Wrep = -1.0*std::pow(a/a0,2.);      // repulsion weight
        Wrep = -1.; // repulsion weight
        Lrep = 0.425 / a; //0.4 / a;  // repulsion length
        WrepPD1_PDL1 = -1.; // repulsion weight on PD1-PDL1
        LrepPD1_PDL1 = 0.425 / a; //0.4 / a;  // repulsion length
        //Lrep = 0.12 / a;
        Watt = 0.;     // attraction weight
        Latt = (R) / a;   // attraction length

//----------- CD2-CD48 -----------------------------------------------------------------------------------------------------------------------------------------------
        Watt_CD2CD48 = 1.;
        Att_CD2CD48 = (0.00) / a; // Self attraction between CD2-CD48 complexes

        Watt_CD2TM = 0.0;
        Att_CD2TM = (0.00) / a;

        Watt_TCRTM = 0.0;
        Att_TCRTM = (0.00) / a; // 0.1 -> to get the diagonal neighbors..., 0.2 -> to get two diagonal neighbors...

        Watt_TMTM = 1.0;
        Att_TMTM = (0.00) / a; // 0.1 -> to get the diagonal neighbors..., 0.2 -> to get two diagonal neighbors...

        Watt_CD2_CD28 = 1.;
        Att_CD2_CD28 = (0.35) / a; // attraction between CD2-CD28 complexes

        Watt_CD2_by_TM = 0.0; // attraction between CD2 and TCR-pMHC complexes
        Att_CD2_by_TM = (0.00)/a;

        Rep_byTM_CD2CD48 = (0.00) / a; // Repulsion between CD2-CD48 and TCR-pMHC complexes
        Wrep_byTM_CD2CD48 = -1.;

        Att_byTM_CD2CD48 = (0.00) / a; // Attraction between CD2-CD48 and TCR-pMHC complexes NOT LIKELY
        Watt_byTM_CD2CD48 = 0.2;
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------
//----------- CD28-CD80 ----------------------------------------------------------------------------------------------------------------------------------------------
        Att_CD28CD80 = (0.21) / a; // attraction between CD28-CD80 complexes
        Watt_CD28CD80 = 1.;

        Rep_byTM_CD28CD80 = (0.00) / a; // Repulsion between CD28-CD80 and TCR-pMHC complexes
        Wrep_byTM_CD28CD80 = -1.;
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------- Centripetal forces -------------------------------------------------------------------------------------------------------------------------------------

        CentrVecTM = 1.0;
        CentrVecLI = 0.00;
        CentrVecCD2CD58 = 0.0;
        CentrVecCD28CD80 = 0.20;
        CentrVecPD1_PDL1 = 0.0;

        OutVecCD2CD58 = -0.0;
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------

        N_A  = 6.02214085774 * pow(10.,23.); // M^(-1) Avogadro number

//------------- for CD45 ----------------------------------------------------------------------------------------------------------------------------------------------
        Watt_by_TM = 0.25;
        Att_by_TM  = (0.00) / a; // maybe i need a little smaller attraction or even higher close to the attraction between TMs....
        Wrep_by_TM = -1.;
        Rep_by_TM  = (0.28) / a; // two gridpoints are affected by the force
        Wrep_by_LI = -1.;
        Rep_by_LI  = (0.0) / a; // With (R-3) or (R-2)/a double localization without the help of Attractive force, with (R-4)/a the biggest ammount is in the cSMAC
        //Rep_by_LI  = (0.06) / a;
                                 // We need both froces in order to have CD45 in cSMAC in the case of LFA-1 defficient mice.

        Wrep_by_CC = -1.;
        Rep_by_CC  = (0.28) / a;

        AttCD45  = 0.4/a;
        WattCD45 = 0.2;

        RepCD45 = 0.14/a;
        WrepCD45 = -1.;
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------

        L_TM = 15. * pow(10.,-3.); // micrometers-> Bond length TCR-pMHC
        L_LI = 45. * pow(10.,-3.); // micrometers-> Bond length LFA1-ICAM1
        Length_of_CD45 = 45 * pow(10.,-3.);
        Koff_TM = 0.1;   // s^(-1)    // On- and Off-Probabilities calculation
        Koff_LI = 0.03;  // s^(-1)
        Kon_TM = 2. * pow(10.,4.) * pow(10.,15.); //2. * pow(10.,4.) * pow(10.,15.); // M^(-1) s^(-1)
        Kon_LI = 3. * pow(10.,5.) * pow(10.,15.); // M^(-1) s^(-1)


//------------- For CD2-CD48 Interactions ------------------------------------------//
        L_CD2_CD48 = 15. * pow(10.,-3.); // same as TCR-pMHC complex

        Koff_CD2_CD48 = 4.;  // maybe it is better for CD2-CD48
        Kon_CD2_CD48 =  10. * pow(10.,5) * pow(10.,15.);
        //Koff_CD2_CD48 = 4.; // maybe better for CD2-CD58
        //Kon_CD2_CD48 = 36. * pow(10.,5) * pow(10.,15.); // maybe go for 2*10^5
        //Koff_CD2_CD48 = 0.1;
        //Kon_CD2_CD48 =  2. * pow(10.,4.) * pow(10.,15.);

//----------------------------------------------------------------------------------//
//------------- For CD28-CD80 Interactions ------------------------------------------//
        L_CD28_CD80 = 15. * pow(10.,-3.); // same as TCR-pMHC complex

        Koff_CD28_CD80 = 1.6;  // maybe it is better for CD2-CD48
        Kon_CD28_CD80 =  4 * pow(10.,5) * pow(10.,15.); //normal
        //Kon_CD28_CD80 =  2 * pow(10.,4) * pow(10.,15.);

//----------------------------------------------------------------------------------//
//------------- For PD1-PDL1 Interactions ------------------------------------------//
        L_PD1_PDL1 = 11. * pow(10.,-3.); // same as TCR-pMHC complex

        Koff_PD1_PDL1 = 1.44;  // maybe it is better for CD2-CD48
        Kon_PD1_PDL1 =  1.84 * pow(10.,5) * pow(10.,15.); //normal
//----------------------------------------------------------------------------------//


        init_TCR  = Area*0.12;
        init_LFA1 = 2.5*init_TCR;
        init_pMHC = Area*0.12;
        init_ICAM = 2.5*init_TCR;
        init_CD45 = Area*0.2;
        init_CD2  = Area*0.15; // when TCR 006
        init_CD48 = Area*0.15;
//        init_CD28 = Area*0.075*0.0; // FullAgents -->  init_CD28 = Area*0.30*0.6
//        init_CD80 = Area*0.075*0.0;//*/

        std::cerr<< "Area = "<<Area<<std::endl;

        std::cerr<< "#TCR = " << init_TCR << '\n' << "#LFA-1 = " << init_LFA1 << '\n' << "#CD2 = " << init_CD2 << '\n' << "#CD28 = " << init_CD28
                 << '\n' << "#pMHC = " << init_pMHC << '\n' << "#ICAM = " << init_ICAM << '\n' << "#CD48 = " << init_CD48 << '\n' << "#CD80 = " << init_CD80 <<std::endl;


        Pnucleate     = 0.0005;
        Pdenucleate   = 0.0005;
        Ppolymerize   = 0.01;
        Pdepolymerize = 0.005;
        Pnucl_polym   = 0.5;
        Pdenucl_polym = 0.00005;
    }


    Toff_TM = (1. / Koff_TM);
    Toff_LI = (1. / Koff_LI);
    Ton_TM  = (L_TM * pow(a0, 2.) * N_A) / Kon_TM;
    Ton_LI  = (L_LI * pow(a0, 2.) * N_A) / Kon_LI;

    //------------- For CD2-CD48 Interactions ------------------------------------------//

    Toff_CD2_CD48 = std::pow(a/a0,2.)*(1. / Koff_CD2_CD48);
    Ton_CD2_CD48  = (L_CD2_CD48 * pow(a, 2.) * N_A) / Kon_CD2_CD48;

    Poff_CD2_CD48 = T/Toff_CD2_CD48;
    Pon_CD2_CD48 = T/Ton_CD2_CD48;

    std::cerr<< "Pon_CD2_CD48 = " << Pon_CD2_CD48 << ", Poff_CD2_CD48 = " << Poff_CD2_CD48 << std::endl;


    //--------------------------------------------------------------------------------//

    //------------- For CD28-CD80 Interactions ------------------------------------------//

    Toff_CD28_CD80 = std::pow(a/a0,2.)*(1 / Koff_CD28_CD80);
    Ton_CD28_CD80  = (L_CD28_CD80 * pow(a, 2.) * N_A) / Kon_CD28_CD80;

    Poff_CD28_CD80 = T/Toff_CD28_CD80;
    Pon_CD28_CD80 = T/Ton_CD28_CD80;

    //--------------------------------------------------------------------------------//


    //------------- For PD1-PDL1 Interactions ------------------------------------------//

    Toff_PD1_PDL1 = std::pow(a/a0,2.)*(1 / Koff_PD1_PDL1);
    Ton_PD1_PDL1  = (L_PD1_PDL1 * pow(a, 2.) * N_A) / Kon_PD1_PDL1;

    Poff_PD1_PDL1 = T/Toff_PD1_PDL1;
    Pon_PD1_PDL1 = T/Ton_PD1_PDL1;

    //--------------------------------------------------------------------------------//

    Tinternalize = 0.45; // Tint = 1 ~> Pint = 0.01225 ~> 0.95 -1.15 %
                         // Tint = 1.75 ~> Pint = 0.007 ~> 0.65 - 0.75 %
                         // Tint = 0.75 ~> Pint = 0.016 ~> 1.20 - 1.35 %
    Pinternalize = T/Tinternalize;
    RangeInternalization = 0.050; // diameter is 100 nm, in micrometers,
    GridDistanceInternalization = RangeInternalization / a; // radius of internalization


    Poff_TM = T/Toff_TM; // probability to unbind TM
    Poff_LI = T/Toff_LI; // probability to unbind LI
    Pon_TM  = T/Ton_TM;  // probability to bind TM
    Pon_LI  = T/Ton_LI;  // probability to bind LI   //int N = pi * pow(R/a,2) + 1<=; // circular interface

    PextTCR  = T/TextTCR; // Pext ~> P externalization
    PextLFA1 = T/TextLFA1;
    PextCD45 = T/TextCD45;
    PextCD2  = T/TextCD2;
    PextCD28  = T/TextCD28;
    PextMHC  = T/TextMHC;
    PextICAM = T/TextICAM;

    P_outFluxCD45 = 0.0005;
    P_outFluxTCR  = 0.0005;
    P_outFluxLFA  = 0.0005;
    P_outFluxMHC  = 0.0005;
    P_outFluxICAM = 0.0005;
    P_outFluxCD2  = 0.0005;
    P_outFluxCD48 = 0.0005;

    P_inFluxCD45 = .05;
    P_inFluxTCR  = .05;
    P_inFluxLFA  = .05;
    P_inFluxMHC  = .05;
    P_inFluxICAM = .05;
    P_inFluxCD2  = .05;
    P_inFluxCD48 = .05;

    std::cerr << "PonTM = " << Pon_TM << ", PonLI = " << Pon_LI << ", PoffTM = " << Poff_TM << ", PoffLI = " <<Poff_LI
              << ", PextTCR = " << PextTCR << std::endl;


    P_near_m = (4*T*Dm*sqrt(2)) / ((pow(a,2))*(sqrt(2)+1)); // molecule to nearest neighbor
    P_near_c = (4*T*Dc*sqrt(2)) / ((pow(a,2))*(sqrt(2)+1)); // complex to nearest neighbor
    P_diag_m = P_near_m / (sqrt(2));      // molecule to nearest neighbor
    P_diag_c = P_near_c / (sqrt(2));      // complex to nearest neaighbor, basically it is LI complex
    std::cout<< "P_near_m = " << P_near_m <<std::endl;    std::cout<< "P_diag_m= " << P_diag_m <<std::endl;
    std::cout<< "P_near_c= " << P_near_c <<std::endl;    std::cout<< "P_diag_c= " << P_diag_c <<std::endl;


    int tcr = 0;
    while (tcr < init_TCR){
        int xt = RandGen(); int yt = RandGen();
        float xc = dim/2; //x center
        float yc = dim/2; //y center
        double distance = sqrt(pow(xc- (float)xt, 2.0)+pow(yc- (float)yt, 2.0));
        if ((distance < SynapseRadius)/* && (distance< dim-2)*/){
            if (Tcell.CreateAgent(TCR, xt, yt)){
                ++tcr;
            }
        }
    }
    int lfa = 0;
    while (lfa < init_LFA1){
        int xl = RandGen(); int yl = RandGen();
        float xc = dim/2; //x center
        float yc = dim/2; //y center
        double distance = sqrt(pow(xc- (float)xl, 2.0)+pow(yc- (float)yl, 2.0));
        if (distance < SynapseRadius){
            if (Tcell.CreateAgent(LFA1, xl, yl)){
                ++lfa;
            }
        }
    }
    int mhc = 0;
    while (mhc < init_pMHC){
        int xm = RandGen(); int ym = RandGen();
        float xc = dim/2; //x center
        float yc = dim/2; //y center
        double distance = sqrt(pow(xc- (float)xm, 2.0)+pow(yc- (float)ym, 2.0));
        if (distance <SynapseRadius){
            if (APC.CreateAgent(pMHC, xm, ym)){
                ++mhc;
            }
        }
    }
    int icam = 0;
    while (icam < init_ICAM){
        int xi = RandGen(); int yi = RandGen();
        float xc = dim/2; //x center
        float yc = dim/2; //y center
        double distance = sqrt(pow(xc- (float)xi, 2.0)+pow(yc- (float)yi, 2.0));
        if (distance < SynapseRadius){
            if (APC.CreateAgent(ICAM1, xi, yi)){
                ++icam;
            }
        }
    }
    int cd45 = 0;
    while (cd45 < init_CD45){
        int xcd = RandGen(); int ycd = RandGen();
        float xc = dim/2; //x center
        float yc = dim/2; //y center
        double distance = sqrt(pow(xc- (float)xcd, 2.0)+pow(yc- (float)ycd, 2.0));
        if ((distance <SynapseRadius)/* && (distance >60)*/){
            if (Tcell.CreateAgent(CD45, xcd, ycd)){
                ++cd45;
            }
        }
    }
    int cd2 = 0;
    while (cd2 < init_CD2){
        int xcd2 = RandGen(); int ycd2 = RandGen();
        float xc = dim/2; //x center
        float yc = dim/2; //y center
        double distance = sqrt(pow(xc- (float)xcd2, 2.0)+pow(yc- (float)ycd2, 2.0));
        if (distance < SynapseRadius){
            if (Tcell.CreateAgent(CD2, xcd2, ycd2)){
                ++cd2;
            }
        }
    }
    int cd28 = 0;
    while (cd28 < init_CD28){
        int xcd28 = RandGen(); int ycd28 = RandGen();
        float xc = dim/2; //x center
        float yc = dim/2; //y center
        double distance = sqrt(pow(xc- (float)xcd28, 2.0)+pow(yc- (float)ycd28, 2.0));
        if (distance < SynapseRadius){
            if (Tcell.CreateAgent(CD28, xcd28, ycd28)){
                ++cd28;
            }
        }
    }
    int cd48 = 0;
    while (cd48 < init_CD48){
        int xcd48 = RandGen(); int ycd48 = RandGen();
        float xc = dim/2; //x center
        float yc = dim/2; //y center
        double distance = sqrt(pow(xc- (float)xcd48, 2.0)+pow(yc- (float)ycd48, 2.0));
        if (distance < SynapseRadius){
            if (APC.CreateAgent(CD48, xcd48, ycd48)){
                ++cd48;
            }
        }
    }
    int cd80 = 0;
    while (cd80 < init_CD80){
        int xcd80 = RandGen(); int ycd80 = RandGen();
        float xc = dim/2; //x center
        float yc = dim/2; //y center
        double distance = sqrt(pow(xc- (float)xcd80, 2.0)+pow(yc- (float)ycd80, 2.0));
        if (distance < SynapseRadius){
            if (APC.CreateAgent(CD80, xcd80, ycd80)){
                ++cd80;
            }
        }
    }
    int emptyT = 0;
    int emptyA = 0;
    for (int i=0; i<dim;++i){
        for(int j=0; j<dim;++j){
            if(Tcell.access(i,j)->getAgent() == nothing)
                ++emptyT;

            if(APC.access(i,j)->getAgent() == nothing)
                ++emptyA;
        }
    }
    std::cerr<< "Empty T = " << emptyT << '\t' << "Empty APC = " << emptyA << std::endl;

    NbInD_5  = 0.0, NbInD_10 = 0.0, NbInD_15 = 0.0, NbInD_20 = 0.0, NbInD_25 = 0.0, NbInD_30 = 0.0, NbInD_35 = 0.0, NbInD_40 = 0.0, NbInD_45 = 0.0, NbInD_50 = 0.0,
    NbInD_55 = 0.0, NbInD_60 = 0.0, NbInD_65 = 0.0, NbInD_70 = 0.0;

    PossiblePositionsTotal = 0.0;
// T R A C K I N G -------------------------------------------------------------------------------------------------------------------------
//    APC.TrackRandom(CD48, 10);
    Tcell.TrackRandom(TCR,10); // when i don't put anything it implicitly is equal to 1, so you will get 2 but we will solve it
//    Tcell.TrackRandom(LFA1,5);
    //APC.TrackRandom(pMHC);

    for (int i=0; i<dim;++i){
        for(int j=0; j<dim;++j){
            double distance = sqrt(pow(xc- (float)i, 2.0)+pow(yc- (float)j, 2.0));
            if (distance < (a0/a)*5.){
                if(Tcell.access(i,j)->getAgent() != border){
                    ++NbInD_5;
                }
            }
            else if (distance < (a0/a)*10.){
                if(Tcell.access(i,j)->getAgent() != border){
                    ++NbInD_10;
                }
            }
            else if (distance < (a0/a)*15.){
                if(Tcell.access(i,j)->getAgent() != border){
                    ++NbInD_15;
                }
            }
            else if (distance < (a0/a)*20.){
                if(Tcell.access(i,j)->getAgent() != border){
                    ++NbInD_20;
                }
            }
            else if (distance < (a0/a)*25.){
                if(Tcell.access(i,j)->getAgent() != border){
                    ++NbInD_25;
                }
            }
            else if (distance < (a0/a)*30.){
                if(Tcell.access(i,j)->getAgent() != border){
                    ++NbInD_30;
                }
            }
            else if (distance < (a0/a)*35.){
                if(Tcell.access(i,j)->getAgent() != border){
                    ++NbInD_35;
                }
            }
            else if (distance < (a0/a)*40.){
                if(Tcell.access(i,j)->getAgent() != border){
                    ++NbInD_40;
                }
            }
            else if (distance < (a0/a)*45.){
                if(Tcell.access(i,j)->getAgent() != border){
                    ++NbInD_45;
                }
            }
            else if (distance < (a0/a)*50.){
                if(Tcell.access(i,j)->getAgent() != border){
                    ++NbInD_50;
                }
            }
            else if (distance < (a0/a)*55.){
                if(Tcell.access(i,j)->getAgent() != border){
                    ++NbInD_55;
                }
            }
            else if (distance < (a0/a)*60.){
                if(Tcell.access(i,j)->getAgent() != border){
                    ++NbInD_60;
                }
            }
            else if (distance < (a0/a)*65.){
                if(Tcell.access(i,j)->getAgent() != border){
                    ++NbInD_65;
                }
            }
            else if (distance < (a0/a)*70.){
                if(Tcell.access(i,j)->getAgent() != border){
                    ++NbInD_70;
                }
            }
        }
    }

   PossiblePositionsTotal = NbInD_5 + NbInD_10 + NbInD_15 + NbInD_20 + NbInD_25 + NbInD_30 + NbInD_35 + NbInD_40 + NbInD_45 + NbInD_50 + NbInD_55 + NbInD_60 + NbInD_65 + NbInD_70;

   std::cout<<"Empty gridpoints on T-cell: "<< emptyT<< std::endl;

}
masterObserver ms = masterObserver("");

bool isComplex(grid_states s){
    return ((s == TM) || (s == LI) || (s == CD2_CD48) || (s == CD28_CD80) || (s == PD1_PDL1));
}

void SwapElements(lattice &Tcell,lattice &APC){
    std::map<gridpoint*, gridpoint*>::iterator itT;
    //std::cerr<< "Unmoved -> " << Tcell.Unmoved.size()<<std::endl;
    int cpr= 0;
    for (itT=Tcell.Unmoved.begin(); itT!=Tcell.Unmoved.end(); ++itT) {
    cpr++;

        gridpoint* oldpos = itT->first;
        gridpoint* targetPos = itT->second;
        grid_states typeOldPos = oldpos->getAgent();
        grid_states typeTargetPos = targetPos->getAgent();

        int oX = oldpos->coord_x; // old X coordinate
        int oY = oldpos->coord_y; // old Y coordinate
        int nX = targetPos->coord_x; // new Y coordinate
        int nY = targetPos->coord_y; // new Y coordinate

        gridpoint* oppositeHere = APC.access(oldpos->coord_x, oldpos->coord_y);
        gridpoint* oppositeTarget = APC.access(targetPos->coord_x, targetPos->coord_y);
        grid_states typeOppositeHere = oppositeHere->getAgent();
        grid_states typeOppositeTarget = oppositeTarget->getAgent();

        std::map<gridpoint*, gridpoint*>::iterator targetFromList = Tcell.Unmoved.find(targetPos); // if the target of the targetPos is oldpos

        if((targetFromList == Tcell.Unmoved.end()) && (typeTargetPos != nothing)) {
            //std::cerr << "Continue biatch" << std::endl;
            continue;
            std::cerr << "Continue biatch - After continue" << std::endl;
        }

        if ((typeOldPos == nothing) || (typeOldPos == border)) {
            continue; // check, go to next for
        }
        if (isComplex(typeOldPos)){ // if you are a complex..
            if ((typeTargetPos == nothing) && (typeOppositeTarget == nothing)){
                Tcell.Move(oldpos->getCoord_x(), oldpos->getCoord_y(), targetPos->getCoord_x(), targetPos->getCoord_y());
                APC.Move(oppositeHere->getCoord_x(), oppositeHere->getCoord_y(), oppositeTarget->getCoord_x(), oppositeTarget->getCoord_y());
            }
            else if ((isComplex(typeTargetPos)) && (targetFromList->second == oldpos) && (typeOldPos!=typeTargetPos) && (RandReal() <= .5)){
                grid_states targetAgent = targetPos->getAgent();
                grid_states oldAgent = oldpos->getAgent();
                grid_states oppositeOldAgent = oppositeHere->getAgent();
                if(Tcell.access(nX, nY)->getAgent() != APC.access(nX, nY)->getAgent()){
                    std::cerr<< "Jesus fuck III trinity is back"<<std::endl;
                }
                Tcell.StateChange(targetAgent, oX, oY);
                Tcell.StateChange(oldAgent, nX, nY);
                APC.StateChange(targetAgent, oX, oY);  // here problem
                APC.StateChange(oppositeOldAgent, nX, nY);

                if(Tcell.access(oX, oY)->getAgent() != APC.access(oX, oY)->getAgent()){
                    std::cerr<< "Crap happens again"<<std::endl;
                    exit(-1);
                }
                if(oldAgent != oppositeOldAgent) std::cerr << "Skata foufoune Power" << std::endl;
            }
        }
    }

    std::map<gridpoint*, gridpoint*>::iterator itA; // = Tcell.Unmoved.begin();
    for (itA=APC.Unmoved.begin(); itA!=APC.Unmoved.end(); ++itA){
        //std::cout << it->first << " => " << it->second << '\n';
        gridpoint* oldpos = itA->first;
        gridpoint* targetPos = itA->second;
        grid_states typeOldPos = oldpos->getAgent();
        grid_states typeTargetPos = targetPos->getAgent();
        //if(typeTargetPos == typeOldPos) continue; // save a lot of time
        gridpoint* oppositeHere = Tcell.access(oldpos->coord_x, oldpos->coord_y);
        gridpoint* oppositeTarget = Tcell.access(targetPos->coord_x, targetPos->coord_y);
        grid_states typeOppositeHere = oppositeHere->getAgent();
        grid_states typeOppositeTarget = oppositeTarget->getAgent();
        //grid_states typeTargetPos = oppositeTarget->getAgent();

        int oX = oldpos->coord_x; // old X coordinate
        int oY = oldpos->coord_y; // old Y coordinate
        int nX = targetPos->coord_x; // new Y coordinate
        int nY = targetPos->coord_y; // new Y coordinate

        if ((typeOldPos == nothing) || (typeOldPos == border)) {
            continue; // check, go to next for
        }
        std::map<gridpoint*, gridpoint*>::iterator targetOfOppositeTarget = Tcell.Unmoved.find(oppositeTarget);
        std::map<gridpoint*, gridpoint*>::iterator targetFromList = APC.Unmoved.find(targetPos);
    }
}


void simulation::Simulate(double TimeEvolution){

    std::ofstream DensityTCR;
    DensityTCR.open("DensityTCR.txt", std::ofstream::out | std::ofstream::app);
    std::ofstream DensityMHC;
    DensityMHC.open("DensityMHC.txt", std::ofstream::out | std::ofstream::app);
    std::ofstream DensityLFA1;
    DensityLFA1.open("DensityLFA1.txt", std::ofstream::out | std::ofstream::app);
    std::ofstream DensityCD45;
    DensityCD45.open("DensityCD45.txt", std::ofstream::out | std::ofstream::app);
    std::ofstream DensitypMHC;
    DensitypMHC.open("DensitypMHC.txt", std::ofstream::out | std::ofstream::app);
    std::ofstream DensityICAM1;
    DensityICAM1.open("DensityICAM.txt", std::ofstream::out | std::ofstream::app);
    std::ofstream DensityTM;
    DensityTM.open("DensityTM.txt", std::ofstream::out | std::ofstream::app);
    std::ofstream DensityLI;
    DensityLI.open("DensityLI.txt", std::ofstream::out | std::ofstream::app);
    std::ofstream DensityCC;
    DensityCC.open("DensityCC.txt", std::ofstream::out | std::ofstream::app);
    std::ofstream DensityCD28CD80;
    DensityCD28CD80.open("DensityCD28CD80.txt", std::ofstream::out | std::ofstream::app);//*/

    std::ofstream DensityTM_CD2CD58coloc;
    DensityTM_CD2CD58coloc.open("DensityTM_CD2CD58coloc.txt", std::ofstream::out | std::ofstream::app);//*/

    std::ofstream InitialAmountOfMolecules;
    InitialAmountOfMolecules.open("InitialAmountOfMolecules.txt", std::ofstream::out | std::ofstream::app);//*/

    std::ofstream freeVSboundMHC;
    freeVSboundMHC.open("freeVSboundMHC.txt", std::ofstream::out | std::ofstream::app);

    std::ofstream freeVSboundLFA;
    freeVSboundLFA.open("freeVSboundLFA.txt", std::ofstream::out | std::ofstream::app);

    std::ofstream freeVSboundCD2;
    freeVSboundCD2.open("freeVSboundCD2.txt", std::ofstream::out | std::ofstream::app);

    double DistFromCentre = 0.0;
    int TCR_5 = 0, TCR_10 = 0, TCR_15 = 0, TCR_20 = 0, TCR_25 = 0, TCR_30 = 0, TCR_35 = 0, TCR_40 = 0, TCR_45 = 0, TCR_50 = 0, TCR_55 = 0, TCR_60 = 0, TCR_65 = 0, TCR_70 = 0;
    int MHC_5 = 0, MHC_10 = 0, MHC_15 = 0, MHC_20 = 0, MHC_25 = 0, MHC_30 = 0, MHC_35 = 0, MHC_40 = 0, MHC_45 = 0, MHC_50 = 0, MHC_55 = 0, MHC_60 = 0, MHC_65 = 0, MHC_70 = 0;
    int LFA_5 = 0, LFA_10 = 0, LFA_15 = 0, LFA_20 = 0, LFA_25 = 0, LFA_30 = 0, LFA_35 = 0, LFA_40 = 0, LFA_45 = 0, LFA_50 = 0, LFA_55 = 0, LFA_60 = 0, LFA_65 = 0, LFA_70 = 0;
    int CD_5 = 0, CD_10 = 0, CD_15 = 0, CD_20 = 0, CD_25 = 0, CD_30 = 0, CD_35 = 0, CD_40 = 0, CD_45 = 0, CD_50 = 0, CD_55 = 0, CD_60 = 0, CD_65 = 0, CD_70 = 0;
    //int MHC_5 = 0, MHC_10 = 0, MHC_15 = 0, MHC_20 = 0, MHC_25 = 0, MHC_30 = 0, MHC_35 = 0, MHC_40 = 0, MHC_45 = 0, MHC_50 = 0, MHC_55 = 0, MHC_60 = 0, MHC_65 = 0, MHC_70 = 0;
    //int ICAM_5 = 0, ICAM_10 = 0, ICAM_15 = 0, ICAM_20 = 0, ICAM_25 = 0, ICAM_30 = 0, ICAM_35 = 0, ICAM_40 = 0, ICAM_45 = 0, ICAM_50 = 0, ICAM_55 = 0, ICAM_60 = 0, ICAM_65 = 0, ICAM_70 = 0;
    int TM_5 = 0, TM_10 = 0, TM_15 = 0, TM_20 = 0, TM_25 = 0, TM_30 = 0, TM_35 = 0, TM_40 = 0, TM_45 = 0, TM_50 = 0, TM_55 = 0, TM_60 = 0, TM_65 = 0, TM_70 = 0;
    int LI_5 = 0, LI_10 = 0, LI_15 = 0, LI_20 = 0, LI_25 = 0, LI_30 = 0, LI_35 = 0, LI_40 = 0, LI_45 = 0, LI_50 = 0, LI_55 = 0, LI_60 = 0, LI_65 = 0, LI_70 = 0;
    int CC_5 = 0, CC_10 = 0, CC_15 = 0, CC_20 = 0, CC_25 = 0, CC_30 = 0, CC_35 = 0, CC_40 = 0, CC_45 = 0, CC_50 = 0, CC_55 = 0, CC_60 = 0, CC_65 = 0, CC_70 = 0;
    int TMCC_5 = 0, TMCC_10 = 0, TMCC_15 = 0, TMCC_20 = 0, TMCC_25 = 0, TMCC_30 = 0, TMCC_35 = 0, TMCC_40 = 0, TMCC_45 = 0, TMCC_50 = 0, TMCC_55 = 0, TMCC_60 = 0, TMCC_65 = 0, TMCC_70 = 0;
    int CD28CD80_5 = 0, CD28CD80_10 = 0, CD28CD80_15 = 0, CD28CD80_20 = 0, CD28CD80_25 = 0, CD28CD80_30 = 0, CD28CD80_35 = 0, CD28CD80_40 = 0, CD28CD80_45 = 0, CD28CD80_50 = 0,
        CD28CD80_55 = 0, CD28CD80_60 = 0, CD28CD80_65 = 0, CD28CD80_70 = 0;//*/


    int OneNeighboringTCRs = 0, TwoNeighboringTCRs = 0;
    int OneNeighboringMHCs = 0, TwoNeighboringMHCs = 0;

    for (double TimeStep = 0.; TimeStep < TimeEvolution; TimeStep+=T){
        chronos += T;
        if ((std::fmod(chronos, 1.))<=T){
            for (int i = 0; i < dim; ++i){
                for (int j = 0; j< dim; ++j){
                    int NoTCRpMHC=0, NoCD28CD80=0, NoCD2CD48=0, NoPD1PDL1=0, NoTCR=0, NopMHC=0;
                    double densityDistance = 0.21/a; // four neighbors
                    double circle = 3.141592654 * densityDistance * densityDistance;

                    for (int k = std::max(0,i-7); k < std::min(dim, i+7); ++k){
                        for(int l = std::max(0,j-7); l < std::min(dim, j+7); ++l){
                            double DistFromAgent = sqrt(pow(double(k-i), 2.0)+pow(double(l-j), 2.0));

                            if (DistFromAgent <= 0.245/a){
                                if (Tcell.Pick(k,l)==TM) {NoTCRpMHC++;}
                                if (Tcell.Pick(k,l)==CD28_CD80) {NoCD28CD80++;}
                                if (Tcell.Pick(k,l)==CD2_CD48) {NoCD2CD48++;}
                                if (Tcell.Pick(k,l)==PD1_PDL1) {NoPD1PDL1++;}
                                if (Tcell.Pick(k,l)==TCR) {NoTCR++;}
                                if (APC.Pick(k,l)==pMHC) {NopMHC++;}
                            }
                        }
                    }

//                    if ((std::min(NoTCRpMHC, NoCD28CD80) > 7) && (Tcell.Pick(i,j)!=border))
//                        if (Colocalization.Pick(i,j)==nothing){
//                            Colocalization.CreateAgent(Nucleation, i,j); //Colocalization
//                        }
//                        else{
//                            Colocalization.StateChange(Nucleation, i,j); //Colocalization
//                        }
//                    else{
//                        if (Colocalization.Pick(i,j)==nothing){
//                            Colocalization.CreateAgent(Polymerization, i,j); //Colocalization
//                        }
//                        else{
//                            Colocalization.StateChange(Polymerization, i,j); //Colocalization
//                        }
//                    }

//                    if ((std::min(NoCD2CD48, NoCD28CD80) > 2) && (Tcell.Pick(i,j)!=border))
//                        if (Colocalization.Pick(i,j)==nothing){
//                            Colocalization.CreateAgent(Nucleation, i,j); //Colocalization
//                        }
//                        else{
//                            Colocalization.StateChange(Nucleation, i,j); //Colocalization
//                        }
//                    else{
//                        if (Colocalization.Pick(i,j)==nothing){
//                            Colocalization.CreateAgent(Polymerization, i,j); //Colocalization
//                        }
//                        else{
//                            Colocalization.StateChange(Polymerization, i,j); //Colocalization
//                        }
//                    }

                    if ((std::min(NoTCRpMHC, NoCD2CD48) > 2) && (Tcell.Pick(i,j)!=border))
                        if (Colocalization.Pick(i,j)==nothing){
                            Colocalization.CreateAgent(CD2_CD48, i,j); //Colocalization
                        }
                        else{
                            Colocalization.StateChange(CD2_CD48, i,j); //Colocalization
                        }
                    else{
                        if (Colocalization.Pick(i,j)==nothing){
                            Colocalization.CreateAgent(nothing, i,j); //Colocalization
                        }
                        else{
                            Colocalization.StateChange(nothing, i,j); //Colocalization
                        }
                    }

//                    if ((std::min(std::min(NoCD28CD80, NoCD2CD48), NoPD1PDL1) > 2) && (Tcell.Pick(i,j)!=border)) // Colocalization of CD2, CD28 and PD1
//                        if (Colocalization.Pick(i,j)==nothing){
//                            Colocalization.CreateAgent(NucleationPolymerized, i,j); //Colocalization
//                        }
//                        else{
//                            Colocalization.StateChange(NucleationPolymerized, i,j); //Colocalization
//                        }
//                    else{
//                        if (Colocalization.Pick(i,j)==nothing){
//                            Colocalization.CreateAgent(nothing, i,j); //Colocalization
//                        }
//                        else{
//                            Colocalization.StateChange(nothing, i,j); //Colocalization
//                        }
//                    }

//                    if ((std::min(NoCD2CD48, NoCD28CD80) > 2) && (Tcell.Pick(i,j)!=border)){ // C Y A N
//                        if (Colocalization.Pick(i,j)==nothing){
//                            Colocalization.CreateAgent(Nucleation, i,j); //Colocalization
//                        }
//                        else{
//                            Colocalization.StateChange(Nucleation, i,j); //Colocalization
//                        }
//                    }
//                    else if ((std::min(NoCD2CD48, NoPD1PDL1) > 2) && (Tcell.Pick(i,j)!=border)){ // R E D
//                        if (Colocalization.Pick(i,j)==nothing){
//                            Colocalization.CreateAgent(NucleationPolymerized, i,j); //Colocalization
//                        }
//                        else{
//                            Colocalization.StateChange(NucleationPolymerized, i,j); //Colocalization
//                        }
//                    }
//                    else if ((std::min(NoCD28CD80, NoPD1PDL1) > 1) && (Tcell.Pick(i,j)!=border)){ // Y E L L O W
//                        if (Colocalization.Pick(i,j)==nothing){
//                            Colocalization.CreateAgent(Polymerization, i,j); //Colocalization
//                        }
//                        else{
//                            Colocalization.StateChange(Polymerization, i,j); //Colocalization
//                        }
//                    }
//                    else{
//                        if (Colocalization.Pick(i,j)==nothing){
//                            Colocalization.CreateAgent(nothing, i,j); //Colocalization
//                        }
//                        else{
//                            Colocalization.StateChange(nothing, i,j); //Colocalization
//                        }
//                    }

//                    if ((std::min(NoCD28CD80, NoCD2CD48) > 2) && (Tcell.Pick(i,j)!=border)){
//                        if (Colocalization.Pick(i,j)==nothing){
//                            Colocalization.CreateAgent(CD28_CD80, i,j); //Colocalization
//                        }
//                        else{
//                            Colocalization.StateChange(CD28_CD80, i,j); //Colocalization
//                        }
//                    }
//                    else{
//                        if (Colocalization.Pick(i,j)==nothing){
//                            Colocalization.CreateAgent(nothing, i,j); //Colocalization
//                        }
//                        else{
//                            Colocalization.StateChange(nothing, i,j); //Colocalization
//                        }
//                    }
//                    else if ((std::min(NoPD1PDL1, NoCD2CD48) > 2) && (Tcell.Pick(i,j)!=border)){
//                        if (Colocalization.Pick(i,j)==nothing){
//                            Colocalization.CreateAgent(CD2_CD48, i,j); //Colocalization
//                        }
//                        else{
//                            Colocalization.StateChange(CD2_CD48, i,j); //Colocalization
//                        }
//                    }
//                    else{
//                        if (Colocalization.Pick(i,j)==nothing){
//                            Colocalization.CreateAgent(nothing, i,j); //Colocalization
//                        }
//                        else{
//                            Colocalization.StateChange(nothing, i,j); //Colocalization
//                        }
//                    }
//                    else if ((std::min(NoPD1PDL1, NoCD28CD80) > 2) && (Tcell.Pick(i,j)!=border)){
//                        if (Colocalization.Pick(i,j)==nothing){
//                            Colocalization.CreateAgent(PD1_PDL1, i,j); //Colocalization
//                        }
//                        else{
//                            Colocalization.StateChange(PD1_PDL1, i,j); //Colocalization
//                        }
//                    }
//                    else{
//                        if (Colocalization.Pick(i,j)==nothing){
//                            Colocalization.CreateAgent(nothing, i,j); //Colocalization
//                        }
//                        else{
//                            Colocalization.StateChange(nothing, i,j); //Colocalization
//                        }
//                    }
                }
            }
        }

        if ((std::fmod(chronos, 30))<=T){
            for (int i=0; i<dim; ++i){
                for (int j=0; j<dim; ++j){
                    if (Tcell.access(i,j)->getAgent() == TCR){
                        DistFromCentre = sqrt(pow(dim/2 - i, 2.0)+pow(dim/2 - j, 2.0));
                        if      (DistFromCentre < (a0/a)*5.)   {TCR_5++; }
                        else if (DistFromCentre < (a0/a)*10.)  {TCR_10++;}
                        else if (DistFromCentre < (a0/a)*15.)  {TCR_15++;}
                        else if (DistFromCentre < (a0/a)*20.)  {TCR_20++;}
                        else if (DistFromCentre < (a0/a)*25.)  {TCR_25++;}
                        else if (DistFromCentre < (a0/a)*30.)  {TCR_30++;}
                        else if (DistFromCentre < (a0/a)*35.)  {TCR_35++;}
                        else if (DistFromCentre < (a0/a)*40.)  {TCR_40++;}
                        else if (DistFromCentre < (a0/a)*45.)  {TCR_45++;}
                        else if (DistFromCentre < (a0/a)*50.)  {TCR_50++;}
                        else if (DistFromCentre < (a0/a)*55.)  {TCR_55++;}
                        else if (DistFromCentre < (a0/a)*60.)  {TCR_60++;}
                        else if (DistFromCentre < (a0/a)*65.)  {TCR_65++;}
                        else if (DistFromCentre <= (a0/a)*70.) {TCR_70++;}
                        else if (DistFromCentre > (a0/a)*70.) {std::cerr<< "Wrong Distance"<<std::endl;}
                        //DensityTCR << chronos << '\t' << Tcell.access(i,j)->getID() <<'\t' << DistFromCentre << std::endl;
                    }
                    if (Tcell.access(i,j)->getAgent() == pMHC){
                        DistFromCentre = sqrt(pow(dim/2 - i, 2.0)+pow(dim/2 - j, 2.0));
                        if      (DistFromCentre < (a0/a)*5.)   {MHC_5++; }
                        else if (DistFromCentre < (a0/a)*10.)  {MHC_10++;}
                        else if (DistFromCentre < (a0/a)*15.)  {MHC_15++;}
                        else if (DistFromCentre < (a0/a)*20.)  {MHC_20++;}
                        else if (DistFromCentre < (a0/a)*25.)  {MHC_25++;}
                        else if (DistFromCentre < (a0/a)*30.)  {MHC_30++;}
                        else if (DistFromCentre < (a0/a)*35.)  {MHC_35++;}
                        else if (DistFromCentre < (a0/a)*40.)  {MHC_40++;}
                        else if (DistFromCentre < (a0/a)*45.)  {MHC_45++;}
                        else if (DistFromCentre < (a0/a)*50.)  {MHC_50++;}
                        else if (DistFromCentre < (a0/a)*55.)  {MHC_55++;}
                        else if (DistFromCentre < (a0/a)*60.)  {MHC_60++;}
                        else if (DistFromCentre < (a0/a)*65.)  {MHC_65++;}
                        else if (DistFromCentre <= (a0/a)*70.) {MHC_70++;}
                        else if (DistFromCentre > (a0/a)*70.) {std::cerr<< "Wrong Distance"<<std::endl;}
                        //DensityMHC << chronos << '\t' << Tcell.access(i,j)->getID() <<'\t' << DistFromCentre << std::endl;
                    }
                    else if (Tcell.access(i,j)->getAgent() == LFA1){
                        DistFromCentre = sqrt(pow(dim/2 - i, 2.0)+pow(dim/2 - j, 2.0));
                        if      (DistFromCentre < (a0/a)*5.)   {LFA_5++; }
                        else if (DistFromCentre < (a0/a)*10.)  {LFA_10++;}
                        else if (DistFromCentre < (a0/a)*15.)  {LFA_15++;}
                        else if (DistFromCentre < (a0/a)*20.)  {LFA_20++;}
                        else if (DistFromCentre < (a0/a)*25.)  {LFA_25++;}
                        else if (DistFromCentre < (a0/a)*30.)  {LFA_30++;}
                        else if (DistFromCentre < (a0/a)*35.)  {LFA_35++;}
                        else if (DistFromCentre < (a0/a)*40.)  {LFA_40++;}
                        else if (DistFromCentre < (a0/a)*45.)  {LFA_45++;}
                        else if (DistFromCentre < (a0/a)*50.)  {LFA_50++;}
                        else if (DistFromCentre < (a0/a)*55.)  {LFA_55++;}
                        else if (DistFromCentre < (a0/a)*60.)  {LFA_60++;}
                        else if (DistFromCentre < (a0/a)*65.)  {LFA_65++;}
                        else if (DistFromCentre <= (a0/a)*70.) {LFA_70++;}
                        else if (DistFromCentre > (a0/a)*70.)  {std::cerr<< "Wrong Distance"<<std::endl;}
                        //DensityLFA1 << chronos << '\t' << Tcell.access(i,j)->getID() <<'\t' << DistFromCentre << std::endl;
                    }
                    else if (Tcell.access(i,j)->getAgent() == CD45){
                        DistFromCentre = sqrt(pow(dim/2 - i, 2.0)+pow(dim/2 - j, 2.0));
                        if      (DistFromCentre < (a0/a)*5.)   {CD_5++; }
                        else if (DistFromCentre < (a0/a)*10.)  {CD_10++;}
                        else if (DistFromCentre < (a0/a)*15.)  {CD_15++;}
                        else if (DistFromCentre < (a0/a)*20.)  {CD_20++;}
                        else if (DistFromCentre < (a0/a)*25.)  {CD_25++;}
                        else if (DistFromCentre < (a0/a)*30.)  {CD_30++;}
                        else if (DistFromCentre < (a0/a)*35.)  {CD_35++;}
                        else if (DistFromCentre < (a0/a)*40.)  {CD_40++;}
                        else if (DistFromCentre < (a0/a)*45.)  {CD_45++;}
                        else if (DistFromCentre < (a0/a)*50.)  {CD_50++;}
                        else if (DistFromCentre < (a0/a)*55.)  {CD_55++;}
                        else if (DistFromCentre < (a0/a)*60.)  {CD_60++;}
                        else if (DistFromCentre < (a0/a)*65.)  {CD_65++;}
                        else if (DistFromCentre <= (a0/a)*70.) {CD_70++;}
                        else if (DistFromCentre > (a0/a)*70.)  {std::cerr<< "Wrong Distance"<<std::endl;}
                        //DensityCD45 << chronos << '\t' << Tcell.access(i,j)->getID() <<'\t' << DistFromCentre << std::endl;
                    }
                    else if (Tcell.access(i,j)->getAgent() == TM){
                        DistFromCentre = sqrt(pow(dim/2 - i, 2.0)+pow(dim/2 - j, 2.0));
                        if      (DistFromCentre < (a0/a)*5.)   {TM_5++; }
                        else if (DistFromCentre < (a0/a)*10.)  {TM_10++;}
                        else if (DistFromCentre < (a0/a)*15.)  {TM_15++;}
                        else if (DistFromCentre < (a0/a)*20.)  {TM_20++;}
                        else if (DistFromCentre < (a0/a)*25.)  {TM_25++;}
                        else if (DistFromCentre < (a0/a)*30.)  {TM_30++;}
                        else if (DistFromCentre < (a0/a)*35.)  {TM_35++;}
                        else if (DistFromCentre < (a0/a)*40.)  {TM_40++;}
                        else if (DistFromCentre < (a0/a)*45.)  {TM_45++;}
                        else if (DistFromCentre < (a0/a)*50.)  {TM_50++;}
                        else if (DistFromCentre < (a0/a)*55.)  {TM_55++;}
                        else if (DistFromCentre < (a0/a)*60.)  {TM_60++;}
                        else if (DistFromCentre < (a0/a)*65.)  {TM_65++;}
                        else if (DistFromCentre <= (a0/a)*70.) {TM_70++;}
                        else if (DistFromCentre > (a0/a)*70.)  {std::cerr<< "Wrong Distance"<<std::endl;}
                        //DensityTM << chronos << '\t' << Tcell.access(i,j)->getID() <<'\t' << DistFromCentre << std::endl;
                    }
                    else if (Tcell.access(i,j)->getAgent() == LI){
                        DistFromCentre = sqrt(pow(dim/2 - i, 2.0)+pow(dim/2 - j, 2.0));
                        if      (DistFromCentre < (a0/a)*5.)   {LI_5++; }
                        else if (DistFromCentre < (a0/a)*10.)  {LI_10++;}
                        else if (DistFromCentre < (a0/a)*15.)  {LI_15++;}
                        else if (DistFromCentre < (a0/a)*20.)  {LI_20++;}
                        else if (DistFromCentre < (a0/a)*25.)  {LI_25++;}
                        else if (DistFromCentre < (a0/a)*30.)  {LI_30++;}
                        else if (DistFromCentre < (a0/a)*35.)  {LI_35++;}
                        else if (DistFromCentre < (a0/a)*40.)  {LI_40++;}
                        else if (DistFromCentre < (a0/a)*45.)  {LI_45++;}
                        else if (DistFromCentre < (a0/a)*50.)  {LI_50++;}
                        else if (DistFromCentre < (a0/a)*55.)  {LI_55++;}
                        else if (DistFromCentre < (a0/a)*60.)  {LI_60++;}
                        else if (DistFromCentre < (a0/a)*65.)  {LI_65++;}
                        else if (DistFromCentre <= (a0/a)*70.) {LI_70++;}
                        else if (DistFromCentre > (a0/a)*70.)  {std::cerr<< "Wrong Distance"<<std::endl;}
                        //DensityLI << chronos << '\t' << Tcell.access(i,j)->getID() <<'\t' << DistFromCentre << std::endl;
                    }
                    else if (Tcell.access(i,j)->getAgent() == CD2_CD48){
                        DistFromCentre = sqrt(pow(dim/2 - i, 2.0)+pow(dim/2 - j, 2.0));
                        if      (DistFromCentre < (a0/a)*5.)   {CC_5++; }
                        else if (DistFromCentre < (a0/a)*10.)  {CC_10++;}
                        else if (DistFromCentre < (a0/a)*15.)  {CC_15++;}
                        else if (DistFromCentre < (a0/a)*20.)  {CC_20++;}
                        else if (DistFromCentre < (a0/a)*25.)  {CC_25++;}
                        else if (DistFromCentre < (a0/a)*30.)  {CC_30++;}
                        else if (DistFromCentre < (a0/a)*35.)  {CC_35++;}
                        else if (DistFromCentre < (a0/a)*40.)  {CC_40++;}
                        else if (DistFromCentre < (a0/a)*45.)  {CC_45++;}
                        else if (DistFromCentre < (a0/a)*50.)  {CC_50++;}
                        else if (DistFromCentre < (a0/a)*55.)  {CC_55++;}
                        else if (DistFromCentre < (a0/a)*60.)  {CC_60++;}
                        else if (DistFromCentre < (a0/a)*65.)  {CC_65++;}
                        else if (DistFromCentre <= (a0/a)*70.) {CC_70++;}
                        else if (DistFromCentre > (a0/a)*70.)  {std::cerr<< "Wrong Distance"<<std::endl;}
                        //DensityCC << chronos << '\t' << Tcell.access(i,j)->getID() <<'\t' << DistFromCentre << std::endl;
                    }
                    else if (Colocalization.access(i,j)->getAgent() == Polymerization){
                        DistFromCentre = sqrt(pow(dim/2 - i, 2.0)+pow(dim/2 - j, 2.0));
                        if      (DistFromCentre < (a0/a)*5.)   {TMCC_5++; }
                        else if (DistFromCentre < (a0/a)*10.)  {TMCC_10++;}
                        else if (DistFromCentre < (a0/a)*15.)  {TMCC_15++;}
                        else if (DistFromCentre < (a0/a)*20.)  {TMCC_20++;}
                        else if (DistFromCentre < (a0/a)*25.)  {TMCC_25++;}
                        else if (DistFromCentre < (a0/a)*30.)  {TMCC_30++;}
                        else if (DistFromCentre < (a0/a)*35.)  {TMCC_35++;}
                        else if (DistFromCentre < (a0/a)*40.)  {TMCC_40++;}
                        else if (DistFromCentre < (a0/a)*45.)  {TMCC_45++;}
                        else if (DistFromCentre < (a0/a)*50.)  {TMCC_50++;}
                        else if (DistFromCentre < (a0/a)*55.)  {TMCC_55++;}
                        else if (DistFromCentre < (a0/a)*60.)  {TMCC_60++;}
                        else if (DistFromCentre < (a0/a)*65.)  {TMCC_65++;}
                        else if (DistFromCentre <= (a0/a)*70.) {TMCC_70++;}
                        else if (DistFromCentre > (a0/a)*70.)  {std::cerr<< "Wrong Distance"<<std::endl;}
                        //DensityTMCC << chronos << '\t' << Tcell.access(i,j)->getID() <<'\t' << DistFromCentre << std::endl;
                    }
                    else if (Tcell.access(i,j)->getAgent() == CD28_CD80){
                        DistFromCentre = sqrt(pow(dim/2 - i, 2.0)+pow(dim/2 - j, 2.0));
                        if      (DistFromCentre < (a0/a)*5.)   {CD28CD80_5++; }
                        else if (DistFromCentre < (a0/a)*10.)  {CD28CD80_10++;}
                        else if (DistFromCentre < (a0/a)*15.)  {CD28CD80_15++;}
                        else if (DistFromCentre < (a0/a)*20.)  {CD28CD80_20++;}
                        else if (DistFromCentre < (a0/a)*25.)  {CD28CD80_25++;}
                        else if (DistFromCentre < (a0/a)*30.)  {CD28CD80_30++;}
                        else if (DistFromCentre < (a0/a)*35.)  {CD28CD80_35++;}
                        else if (DistFromCentre < (a0/a)*40.)  {CD28CD80_40++;}
                        else if (DistFromCentre < (a0/a)*45.)  {CD28CD80_45++;}
                        else if (DistFromCentre < (a0/a)*50.)  {CD28CD80_50++;}
                        else if (DistFromCentre < (a0/a)*55.)  {CD28CD80_55++;}
                        else if (DistFromCentre < (a0/a)*60.)  {CD28CD80_60++;}
                        else if (DistFromCentre < (a0/a)*65.)  {CD28CD80_65++;}
                        else if (DistFromCentre <= (a0/a)*70.) {CD28CD80_70++;}
                        else if (DistFromCentre > (a0/a)*70.)  {std::cerr<< "Wrong Distance"<<std::endl;}
                        //DensityCD28CD80 << chronos << '\t' << Tcell.access(i,j)->getID() <<'\t' << DistFromCentre << std::endl;
                    }
                    if (Tcell.access(i,j)->getAgent() == TCR){
                        bool foundOneNeighbor = false;
                        bool foundTwoNeighbors = false;
                        for (int k=-1; k<=1; ++k){
                            for (int l=-1; l<=1; ++l){
                                if (Tcell.access(i+k,j+l)->getAgent() == TM){
                                    if(!foundOneNeighbor) {
                                        OneNeighboringTCRs++;
                                        foundOneNeighbor = true;
                                    }
                                }
                            }
                        }
                        for (int k=-2; k<=2; ++k){
                            for (int l=-2; l<=2; ++l){
                                if (Tcell.access(i+k,j+l)->getAgent() == TM){
                                    if(!foundTwoNeighbors) {
                                        TwoNeighboringTCRs++;
                                        foundTwoNeighbors = true;
                                    }
                                }
                            }
                        }
                    }
                    if (APC.access(i,j)->getAgent() == pMHC){
                        bool foundOneMHCNeighbor = false;
                        bool foundTwoMHCNeighbors = false;
                        for (int k=-1; k<=1; ++k){
                            for (int l=-1; l<=1; ++l){
                                if (Tcell.access(i+k,j+l)->getAgent() == TM){
                                    if(!foundOneMHCNeighbor) {
                                        OneNeighboringMHCs++;
                                        foundOneMHCNeighbor = true;
                                    }
                                }
                            }
                        }
                        for (int k=-2; k<=2; ++k){
                            for (int l=-2; l<=2; ++l){
                                if (Tcell.access(i+k,j+l)->getAgent() == TM){
                                    if(!foundTwoMHCNeighbors) {
                                        TwoNeighboringMHCs++;
                                        foundTwoMHCNeighbors = true;
                                    }
                                }
                            }
                        }
                    }
                }
            }//*/



//            std::cerr << chronos << '\t' <<Tcell.Count(TCR) << '\t' << APC.Count(pMHC) << '\t' << Tcell.Count(TM) << '\t'
//                      << OneNeighboringTCRs << '\t' << Tcell.Count(TM)+OneNeighboringTCRs << '\t' << TwoNeighboringTCRs << '\t'
//                      << Tcell.Count(TM)+TwoNeighboringTCRs  << '\t' << OneNeighboringMHCs << '\t' << Tcell.Count(TM)+OneNeighboringMHCs << '\t'
//                      << TwoNeighboringMHCs << '\t' << Tcell.Count(TM)+TwoNeighboringMHCs << std::endl;

            freeVSboundMHC << chronos << '\t' <<Tcell.Count(TCR) << '\t' << APC.Count(pMHC) << '\t' << Tcell.Count(TM) << std::endl;
            freeVSboundLFA << chronos << '\t' <<Tcell.Count(LFA1) << '\t' << APC.Count(ICAM1) << '\t' << Tcell.Count(LI) << std::endl;
            /* '\t' << Tcell.Count(TM)+OneNeighboringTCRs << '\t' << Tcell.Count(TM)+TwoNeighboringTCRs  << '\t'
                           << Tcell.Count(TM)+OneNeighboringMHCs << '\t' << Tcell.Count(TM)+TwoNeighboringMHCs << std::endl;*/

            freeVSboundCD2 << chronos << '\t' <<Tcell.Count(CD2) /*<< '\t' << APC.Count(CD48)*/ << '\t' << Tcell.Count(CD2_CD48) << std::endl;


            freeVSboundMHC << chronos << '\t' <<Tcell.Count(TCR) << '\t' << APC.Count(pMHC) << '\t' << Tcell.Count(TM) << std::endl;
//                           << '\t' << (double)Tcell.Count(TM)/APC.Count(pMHC) << '\t' << (double)APC.Count(pMHC)/Tcell.Count(TM)
//                           << '\t' << (double)APC.Count(pMHC)/(Tcell.Count(TM)+APC.Count(pMHC)) << '\t' << (double)Tcell.Count(TM)/(Tcell.Count(TM)+APC.Count(pMHC))
//                           << '\t' << cellKD << '\t' << "LFA-1" << '\t' <<  (double)Tcell.Count(LFA1)/(Tcell.Count(LI)+Tcell.Count(LFA1)) << '\t'
//                           << (double)Tcell.Count(LI)/(Tcell.Count(LI)+Tcell.Count(LFA1)) << '\t' << LI_cellKD << std::endl;
//            // APC.Count(pMHC)/Tcell.Count(TM) cannot use this because initially TM=0 and does not like division by zero.....

            DensityTCR << chronos << '\t' << TCR_5/NbInD_5 << '\t' << TCR_10/NbInD_10 << '\t' << TCR_15/NbInD_15 << '\t' << TCR_20/NbInD_20 << '\t' << TCR_25/NbInD_25 << '\t'
                       << TCR_30/NbInD_30 << '\t' << TCR_35/NbInD_35 << '\t' << TCR_40/NbInD_40 << '\t' << TCR_45/NbInD_45  << '\t' << TCR_50/NbInD_50 << '\t' << TCR_55/NbInD_55
                       << '\t' << TCR_60/NbInD_60 << '\t' << TCR_65/NbInD_65 << '\t' << TCR_70/NbInD_70 << std::endl;

            DensityMHC << chronos << '\t' << MHC_5/NbInD_5 << '\t' << MHC_10/NbInD_10 << '\t' << MHC_15/NbInD_15 << '\t' << MHC_20/NbInD_20 << '\t' << MHC_25/NbInD_25 << '\t'
                       << MHC_30/NbInD_30 << '\t' << MHC_35/NbInD_35 << '\t' << MHC_40/NbInD_40 << '\t' << MHC_45/NbInD_45  << '\t' << MHC_50/NbInD_50 << '\t' << MHC_55/NbInD_55
                       << '\t' << MHC_60/NbInD_60 << '\t' << MHC_65/NbInD_65 << '\t' << MHC_70/NbInD_70 << std::endl;

            DensityLFA1 << chronos << '\t' << LFA_5/NbInD_5 << '\t' << LFA_10/NbInD_10 << '\t' << LFA_15/NbInD_15 << '\t' << LFA_20/NbInD_20 << '\t' << LFA_25/NbInD_25 << '\t'
                        << LFA_30/NbInD_30 << '\t' << LFA_35/NbInD_35 << '\t' << LFA_40/NbInD_40 << '\t' << LFA_45/NbInD_45  << '\t' << LFA_50/NbInD_50  << '\t' << LFA_55/NbInD_55
                        << '\t' << LFA_60/NbInD_60 << '\t' << LFA_65/NbInD_65 << '\t' << LFA_70/NbInD_70 << '\t' << std::endl;//*/

            DensityCD45 << chronos << '\t' << CD_5/NbInD_5 << '\t' << CD_10/NbInD_10 << '\t' << CD_15/NbInD_15 << '\t' << CD_20/NbInD_20 << '\t' << CD_25/NbInD_25  << '\t'
                        << CD_30/NbInD_30 << '\t' << CD_35/NbInD_35 << '\t' << CD_40/NbInD_40 << '\t' << CD_45/NbInD_45  << '\t' << CD_50/NbInD_50  << '\t' << CD_55/NbInD_55
                        << '\t' << CD_60/NbInD_60 << '\t' << CD_65/NbInD_65 << '\t' << CD_70/NbInD_70 << '\t' << std::endl;//*/

            DensityTM << chronos << '\t' << TM_5/NbInD_5 << '\t' << TM_10/NbInD_10 << '\t' << TM_15/NbInD_15 << '\t' << TM_20/NbInD_20 << '\t' << TM_25/NbInD_25 << '\t'
                      << TM_30/NbInD_30 << '\t' << TM_35/NbInD_35 << '\t' << TM_40/NbInD_40 << '\t' << TM_45/NbInD_45  << '\t' << TM_50/NbInD_50  << '\t' << TM_55/NbInD_55
                      << '\t' << TM_60/NbInD_60 << '\t' << TM_65/NbInD_65 << '\t' << TM_70/NbInD_70 << '\t' << std::endl;

            DensityLI << chronos << '\t' << LI_5/NbInD_5 << '\t' << LI_10/NbInD_10 << '\t' << LI_15/NbInD_15 << '\t' << LI_20/NbInD_20 << '\t' << LI_25/NbInD_25  << '\t'
                      << LI_30/NbInD_30 << '\t' << LI_35/NbInD_35 << '\t' << LI_40 /NbInD_40<< '\t' << LI_45/NbInD_45  << '\t' << LI_50/NbInD_50  << '\t' << LI_55/NbInD_55
                      << '\t' << LI_60/NbInD_60 << '\t' << LI_65/NbInD_65 << '\t' << LI_70/NbInD_70 << '\t' << std::endl;

            DensityCC << chronos << '\t' << CC_5/NbInD_5 << '\t' << CC_10/NbInD_10 << '\t' << CC_15/NbInD_15 << '\t' << CC_20/NbInD_20 << '\t' << CC_25/NbInD_25  << '\t'
                      << CC_30/NbInD_30 << '\t' << CC_35/NbInD_35 << '\t' << CC_40 /NbInD_40<< '\t' << CC_45/NbInD_45  << '\t' << CC_50/NbInD_50  << '\t' << CC_55/NbInD_55
                      << '\t' << CC_60/NbInD_60 << '\t' << CC_65/NbInD_65 << '\t' << CC_70/NbInD_70 << '\t' << std::endl;

            DensityTM_CD2CD58coloc << chronos << '\t' << TMCC_5/NbInD_5 << '\t' << TMCC_10/NbInD_10 << '\t' << TMCC_15/NbInD_15 << '\t' << TMCC_20/NbInD_20 << '\t' << TMCC_25/NbInD_25  << '\t'
                      << TMCC_30/NbInD_30 << '\t' << TMCC_35/NbInD_35 << '\t' << TMCC_40 /NbInD_40<< '\t' << TMCC_45/NbInD_45  << '\t' << TMCC_50/NbInD_50  << '\t' << TMCC_55/NbInD_55
                      << '\t' << TMCC_60/NbInD_60 << '\t' << TMCC_65/NbInD_65 << '\t' << TMCC_70/NbInD_70 << '\t' << std::endl;


            DensityCD28CD80 << chronos << '\t' << CD28CD80_5/NbInD_5 << '\t' << CD28CD80_10/NbInD_10 << '\t' << CD28CD80_15/NbInD_15 << '\t' << CD28CD80_20/NbInD_20 << '\t'
                            << CD28CD80_25/NbInD_25  << '\t' << CD28CD80_30/NbInD_30 << '\t' << CD28CD80_35/NbInD_35 << '\t' << CD28CD80_40 /NbInD_40<< '\t' << CD28CD80_45/NbInD_45
                            << '\t' << CD28CD80_50/NbInD_50  << '\t' << CD28CD80_55/NbInD_55 << '\t' << CD28CD80_60/NbInD_60 << '\t' << CD28CD80_65/NbInD_65 << '\t'
                            << CD28CD80_70/NbInD_70 << '\t' << std::endl; //*/
        }


//-------------------------------------------------------------------------------------------------------------------------------------------------------------//

        double BindCoef = 0.0; // coefficient that affects the binding rates of TCR-pMHC and LFA-1-ICAM-1
        int coeffTM = 0;
        int coeff_CD2_CD48 = 0;
        //std::cerr<< Prob<< std::endl;

#ifdef USE_SHUFFLE
        std::vector<gridpoint*> toShuffle;
        toShuffle.insert(toShuffle.end(), Tcell.GetListOfAgents(TCR)->begin(), Tcell.GetListOfAgents(TCR)->end());
        toShuffle.insert(toShuffle.end(), Tcell.GetListOfAgents(LFA1)->begin(), Tcell.GetListOfAgents(LFA1)->end());
        toShuffle.insert(toShuffle.end(), Tcell.GetListOfAgents(CD45)->begin(), Tcell.GetListOfAgents(CD45)->end());
        toShuffle.insert(toShuffle.end(), Tcell.GetListOfAgents(TM)->begin(), Tcell.GetListOfAgents(TM)->end());
        toShuffle.insert(toShuffle.end(), Tcell.GetListOfAgents(LI)->begin(), Tcell.GetListOfAgents(LI)->end());
        toShuffle.insert(toShuffle.end(), Tcell.GetListOfAgents(CD2)->begin(), Tcell.GetListOfAgents(CD2)->end());
        toShuffle.insert(toShuffle.end(), Tcell.GetListOfAgents(CD2_CD48)->begin(), Tcell.GetListOfAgents(CD2_CD48)->end());
        toShuffle.insert(toShuffle.end(), Tcell.GetListOfAgents(CD28)->begin(), Tcell.GetListOfAgents(CD28)->end());
        toShuffle.insert(toShuffle.end(), Tcell.GetListOfAgents(CD28_CD80)->begin(), Tcell.GetListOfAgents(CD28_CD80)->end());
        toShuffle.insert(toShuffle.end(), APC.GetListOfAgents(ICAM1)->begin(), APC.GetListOfAgents(ICAM1)->end());
        toShuffle.insert(toShuffle.end(), APC.GetListOfAgents(pMHC)->begin(), APC.GetListOfAgents(pMHC)->end());
        toShuffle.insert(toShuffle.end(), APC.GetListOfAgents(CD48)->begin(), APC.GetListOfAgents(CD48)->end());
        toShuffle.insert(toShuffle.end(), APC.GetListOfAgents(CD80)->begin(), APC.GetListOfAgents(CD80)->end());


        // maybe the best way to randomly shuffle..
        std::random_device rng;
        std::mt19937 urng(rng());
        std::shuffle(toShuffle.begin(), toShuffle.end(), urng);

        // this function will be deprecated in C++14 and removed in C++17
        //std::random_shuffle(toShuffle.begin(), toShuffle.end());

        Tcell.Unmoved.clear();
        APC.Unmoved.clear();

        int totalNbAgentsToCheck = toShuffle.size();
        for(int check = 0; check < totalNbAgentsToCheck; ++check){
            int x = toShuffle[check]->getCoord_x();
            int y = toShuffle[check]->getCoord_y();

#else
        for (int checks=0; checks <  3*(Tcell.Count() + APC.Count()-APC.Count(TM)-APC.Count(LI)); ++checks) {
            int x = RandGen();              int y = RandGen();
            double Prob = ((double)Tcell.Count())/((double)(Tcell.Count()+APC.Count()-APC.Count(TM)-APC.Count(LI)));
#endif


//--------------------B I N D I N G  and  U N B I N D I N G  and M O V E M E N T  and  I N T E R A C T I O N S------------------------------------


#ifdef USE_SHUFFLE
            switch(toShuffle[check]->getAgent()){
#else
            double RandNum = RandReal();
            if (RandNum <= Prob){
                int count = 40000;
                while((count > 0) && (Tcell.Pick(x,y) == nothing) && (Tcell.Pick(x,y) == border)) { // if the site is empty or border pick new X and Y
                    x = RandGen();
                    y = RandGen();
                    std::cerr << "Fuck Fuck Fuck Fuck Fuck " << count << std::endl;
                    count--;
                }
                switch(Tcell.Pick(x,y)){
#endif
                case CD45:{
                    double P_near_CD45 = P_near_m;
                    double P_diag_CD45 = P_diag_m;
                    int ACTION = RandInteraction();
                    if (ACTION == diffusion){
                        double RR = RandReal();
                        if (RR <= P_near_CD45 + P_diag_CD45){
                            int delta = 0;
                            int epsilon = 0;
                            Tcell.getDelta_Epsilon(delta, epsilon);
                            Tcell.Move(x, y, x+delta, y+epsilon);
                        }
                    }
                    if (ACTION == complex_forces){
//                        if (Tcell.CheckNeighbors(x, y) > 7){
//                            break;
//                        }
//                        else if(Tcell.CheckNeighbors(x, y)<=7){
                            Vector2 Sum = Vector2(0,0);
                            for (int i=std::max(0, x-(int(Lrep)+1)); i<std::min(dim, x+(int(Lrep)+1)); ++i){
                                for (int j=std::max(0, y-(int(Lrep)+1)); j< std::min(dim, y+(int(Lrep)+1)); ++j){
                                    if (Tcell.Pick(i, j) == LI){
                                        double dist_CD45LI = sqrt(pow(double(x-i), 2.0)+pow(double(y-j), 2.0));
                                        if (dist_CD45LI < Rep_by_LI){
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Wrep_by_LI));}
                                    }//*/
                                    /*if (Tcell.Pick(i, j) == TM){
                                        double dist_CD45TM = sqrt(pow(double(x-i), 2.0)+pow(double(y-j), 2.0));
                                        if ((dist_CD45TM < Att_by_TM)){ // IN CASE YOU PICK YOUR OWN POSITION
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Watt_by_TM));} // what weight should i add?
                                    }//*/
                                    if (Tcell.Pick(i, j) == TM){
                                        double dist_CD45TM = sqrt(pow(double(x-i), 2.0)+pow(double(y-j), 2.0));
                                        if ((dist_CD45TM < Rep_by_TM)){
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Wrep_by_TM));} // what weight should i add?
                                    }//*/
                                    if (Tcell.Pick(i, j) == CD2_CD48){
                                        double dist_CD45CC = sqrt(pow(double(x-i), 2.0)+pow(double(y-j), 2.0));
                                        if (dist_CD45CC < Rep_by_CC){
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Wrep_by_CC));}
                                    }//*/
                                    /*if (Tcell.Pick(i, j) == CD45){
                                        double dist_CD45 = sqrt(pow(double(x-i), 2.0)+pow(double(y-j), 2.0));
                                        if ((dist_CD45 < RepCD45)){
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(WrepCD45));} // what weight should i add?
                                    }//*/
                                    /*if (Tcell.Pick(i, j) == CD45){
                                        double dist_CD45 = sqrt(pow(double(x-i), 2.0)+pow(double(y-j), 2.0));
                                        if ((dist_CD45 < AttCD45)){
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Tcell.Lennard_Jones(AttCD45,WattCD45,dist_CD45)));} // what weight should i add?
                                    }//*/
                                }
                            }
                            double RandRealNumber = RandReal();
                            if ((Sum.Length() > 1e-12) && (RandRealNumber <= (P_diag_CD45 + P_near_CD45)*(a/a0)*Sum.Length())){
                                int delta = 0;
                                int epsilon = 0;
                                Tcell.getDelta_Epsilon_from_Angle(delta, epsilon, Sum);
                                if (Tcell.Pick(x+delta, y+epsilon) == nothing){
                                    Tcell.Move(x, y, x+delta, y+epsilon);
                                }
//                                else {
//                                    Tcell.Unmoved.insert(std::pair<gridpoint*, gridpoint*>(Tcell.access(x,y), Tcell.access(x+delta, y+epsilon)));
//                                    APC.Unmoved.insert(std::pair<gridpoint*, gridpoint*>(APC.access(x,y), APC.access(x+delta, y+epsilon)));
//                                }
//                            }
                        }
                    } // end else*/
                    break;}
                case CD2:{
                    double T_SameNeighbor = 0.0;
                    double TCRmediated = 0.0;
                    double CoeffBinding = 1.0;
                    int AdhesionSites = 1;
                    int SameNeighbors = 2;
                    for (int i=-AdhesionSites; i<=AdhesionSites; ++i){
                        for (int j=-AdhesionSites; j<=AdhesionSites; ++j){
                            if((Tcell.Pick(x+i, y+j) == TM) ){
                                ++TCRmediated;
                            }
                        }
                    }
////-------------------------------------------------------------------------------------------------------
           if (TCRmediated>=SameNeighbors){
               CoeffBinding = 1000*(1.+((double)TCRmediated))*(1.+((double)TCRmediated))*(1.+((double)TCRmediated))*(1.+((double)TCRmediated))*(1.+((double)TCRmediated));
           }
           else if (TCRmediated<SameNeighbors){
               CoeffBinding = 1.;
           }
//                    if (T_SameNeighbor>=SameNeighbors){
//                        CoeffBinding = 0.00001*(1.+((double)T_SameNeighbor))*(1.+((double)T_SameNeighbor))*(1.+((double)T_SameNeighbor))*(1.+((double)T_SameNeighbor));
//                    }
//                    else if (T_SameNeighbor<SameNeighbors){
//                        CoeffBinding = 1.;
//                    }
//-------------------------------------------------------------------------------------------------------
//                    if (T_SameNeighbor<SameNeighbors){
//                        f_Nnn = 1.;
//                    }
//                    else if (T_SameNeighbor>=SameNeighbors){
//                        f_Nnn = 1./(((double)T_SameNeighbor-2)*((double)T_SameNeighbor-2));
//                    }
//-------------------------------------------------------------------------------------------------------
//                    f_Nnn = 1./(1+(double)T_SameNeighbor*T_SameNeighbor);
//                    std::cerr<< "T_SameNeighbor" << '\t' << T_SameNeighbor << '\t' << "f_Nnn"<< '\t' << f_Nnn
//                             << '\t' << "P" << '\t' <<(P_near_c + P_diag_c)*f_Nnn <<std::endl;
//-------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------
//                    coeff_CD2_CD48 = 1;
//                    int negativeRegulation = 0;
//                    double Binding = 0.;
//                    int NeighoringSites = 1; // to define the neighboring sites faster for the two 'for' loops/*
//                    for (int i=-NeighoringSites; i<=NeighoringSites;++i){
//                        for (int j=-NeighoringSites; j<=NeighoringSites;++j){
//                            if (Tcell.Pick(x+i,y+j) == CD2_CD48){
//                                coeff_CD2_CD48++;
//                            }
//                            if (Tcell.Pick(x+i,y+j) == LI){
//                                negativeRegulation++;
//                            }
//                        }
//                    }
//                    if(coeff_CD2_CD48 <=3){
//                        Binding = 1.;
//                    }
//                    else if (coeff_CD2_CD48 >3){
//                        Binding = 0.000001*(double)coeff_CD2_CD48;
//                    }
//                    if (coeff_CD2_CD48>3)
//                        std::cerr<< "coeff_CD2_CD48 = "<< coeff_CD2_CD48<< '\t' << "Binding = " << Binding << std::endl;
//                    double BindingProbability = Pon_CD2_CD48*Binding;
//-------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------
                    BindCoef = 1.;
                    int ACTION = RandInteraction();

//                    if (BindingProbability>Pon_CD2_CD48)
//                        std::cerr<<"BindingProbability = "<<BindingProbability<<std::endl;
                    if (/*Action < 0.3*/ ACTION == binding_kinetics){
//                            if (Pon_CD2_CD48 * (double)coeff_CD2_CD48 <= 1.0){
                        if((APC.Pick(x, y) == CD48) && (RandReal() <  /*BindingProbability*/ Pon_CD2_CD48/**CoeffBinding*/ /* Binding*/)){ // check the opposite lattice in the same position
//                            if (Pon_CD2_CD48*CoeffBinding > Pon_CD2_CD48)
//                                std::cerr << "Pon_CD2_CD48" << '\t' << Pon_CD2_CD48<< '\t' << Pon_CD2_CD48 * (double)CoeffBinding <<std::endl;
                            Tcell.StateChange(CD2_CD48,x,y);
                            APC.StateChange(CD2_CD48,x,y);
                        }
                        //                            }
                    }
                    else if (ACTION == diffusion){
                        double RR = RandReal();
                       if ((NoDTDepDiffusion) || (RR<P_near_m+P_diag_m)){
                            int delta = 0;
                            int epsilon = 0;
                            Tcell.getDelta_Epsilon(delta, epsilon);
                            Tcell.Move(x, y,x+delta, y+epsilon);
                        }
                    }
                    break;}
                case CD28:{
                    BindCoef = 1.;
                    int ACTION = RandInteraction();
                    if (ACTION == binding_kinetics){
                        if((APC.Pick(x, y) == CD80) && (RandReal() <  Pon_CD28_CD80*BindCoef)){ // check the opposite lattice in the same position
                            Tcell.StateChange(CD28_CD80,x,y);
                            APC.StateChange(CD28_CD80,x,y);
                        }
                    }
                    if (ACTION == diffusion){
                        double RR = RandReal();
                       if ((NoDTDepDiffusion) || (RR<P_near_m+P_diag_m)){
                            int delta = 0;
                            int epsilon = 0;
                            Tcell.getDelta_Epsilon(delta, epsilon);
                            Tcell.Move(x, y,x+delta, y+epsilon);
                        }
                    }
                    break;}
                case TCR:{
                    BindCoef = 1.;
                    coeffTM  = 1;//Tcell.Count(TM);
//                    for (int i=-2; i<=2;++i){
//                        for (int j=-2; j<=2;++j){
//                            if (Tcell.Pick(x+i,y+j) == TM){
//                                coeffTM++;
//                            }
//                        }
//                    }
//                    if (coeffTM>1){
//                        std::cerr << "coeffTM = " << coeffTM << '\t' <<(double)coeffTM <<std::endl;
//                    }
                    int ACTION = RandInteraction();
                    //double Action = RandReal();
//                    if (ACTION == binding_kinetics){
//                        if((APC.Pick(x, y) == pMHC) && (RandReal() < Pon_TM*BindCoef)){ // check the opposite lattice in the same position
//                            Tcell.StateChange(TM,x,y);
//                            APC.StateChange(TM,x,y);
//                        }
//                    }
                    if (ACTION == binding_kinetics){
                        if((APC.Pick(x, y) == pMHC) && (RandReal() < Pon_TM * (double)coeffTM)){ // check the opposite lattice in the same position
                            Tcell.StateChange(TM,x,y);
                            APC.StateChange(TM,x,y);
                        }
                    }
                    else if (ACTION == diffusion){
                        double RR = RandReal();
                       if ((NoDTDepDiffusion) || (RR<P_near_m+P_diag_m)){
                            int delta = 0;
                            int epsilon = 0;
                            Tcell.getDelta_Epsilon(delta, epsilon);
                            Tcell.Move(x, y,x+delta, y+epsilon);
                            //Tracker.StateChange(TCR,x,y);
                        }
                    }//*/
                    break;}
                case LFA1:{
                    BindCoef = 1.;
                    int ACTION = RandInteraction();
                    //double Action = RandReal();
                    if (/*Action < 0.3*/ ACTION == binding_kinetics){
                        if((APC.Pick(x, y) == ICAM1) && (RandReal() <  Pon_LI*BindCoef)){
                            Tcell.StateChange(LI,x,y);
                            APC.StateChange(LI,x,y);;
                        }
                    }
                    else if (/*(Action > 0.3) && (Action < 0.8)*/ ACTION == diffusion){
                        double RR = RandReal();
                        if ((NoDTDepDiffusion) || (RR<P_near_m+P_diag_m)){
                            int delta = 0;
                            int epsilon = 0;
                            Tcell.getDelta_Epsilon(delta, epsilon);
                            Tcell.Move(x, y,x+delta, y+epsilon);
                        }
                    }//*/
                    break;}
                case TM:{
                    if(APC.Pick(x,y) != TM) std::cerr << "No TM" << std::endl;
                    double T_SameNeighbor = 0.0;
                    for (int i=-1; i<=1; ++i){
                        for (int j=-1; j<=1; ++j){
                            if( /*((((i!=0) && (j==0)) || ((i==0) &&(j!=0))) && !((i==0) && (j==0)) ) &&*/ (Tcell.Pick(x+i, y+j) == TM) ){ // (i*j !=0)
                                ++T_SameNeighbor;
                                //std::cerr<<T_SameNeighbor<<std::endl;
                            }
                        }
                    }
                    //double f_Nnn = 1/(1+T_SameNeighbor);
                    double P_near_TM = ((4*sqrt(2)*T*Dc) / (pow(a,2)*(sqrt(2)+1)));//*(f_Nnn);
                    double P_diag_TM = P_near_TM / sqrt(2);
                    //if (P_near_TM < 0.017 ) std::cerr<<"P_near_TM = "<<P_near_TM <<std::endl;
                    //std::cerr<<"Ptm="<<a*T*(P_diag_TM + P_near_TM)<<std::endl;
                    int ACTION = RandInteraction();
                    if (ACTION == binding_kinetics){
                        if( RandReal() < Poff_TM){ // Unbind
                            Tcell.StateChange(TCR,x,y);
                            APC.StateChange(pMHC,x,y);
                        }
                    }
                    else if (ACTION == diffusion){
                        double RR = RandReal();
                        if ( (RR <= P_near_TM + P_diag_TM) ){
                            int delta = 0;
                            int epsilon = 0;
                            Tcell.getDelta_Epsilon(delta, epsilon);
                            if ( (Tcell.Pick(x+delta, y+epsilon) == nothing) && (APC.Pick(x+delta, y+epsilon) == nothing) ){
                                Tcell.Move(x, y,x+delta, y+epsilon);
                                APC.Move(x, y,x+delta, y+epsilon);
                            }
                            else {
                                Tcell.Unmoved.insert(std::pair<gridpoint*, gridpoint*>(Tcell.access(x,y), Tcell.access(x+delta, y+epsilon)));
                                APC.Unmoved.insert(std::pair<gridpoint*, gridpoint*>(APC.access(x,y), APC.access(x+delta, y+epsilon)));
                            }
                        }
                    }//*/
                    else if (ACTION == complex_forces){
                            Vector2 Sum = Vector2(0,0);
                            //Sum.Add(Vector2(x, y, x+delta, y+epsilon));
                            for (int i=std::max(0, x-(int(Lrep)+1)); i<std::min(dim, x+(int(Lrep)+1)); ++i){
                                for (int j=std::max(0, y-(int(Lrep)+1)); j< std::min(dim, y+(int(Lrep)+1)); ++j){
                                    /*if (Tcell.Pick(i, j) == TM){ // TM - PD1 Attraction
                                        double dist_TM_TCR = sqrt(pow(double(x-i), 2.0)+pow(double(y-j), 2.0));
                                        if ((dist_TM_TCR > 0) && (dist_TM_TCR < Att_TMTM)){
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Watt_TMTM));}
                                    }//*/
                                    if (Tcell.Pick(i, j) == LI){
                                        double dist_TMLI = sqrt(pow(double(x-i), 2.0)+pow(double(y-j), 2.0));
                                        if (dist_TMLI < Lrep){
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Wrep));}
                                    }//*/
                                    /*if (Tcell.Pick(i, j) == CD2_CD48){ // TM - PD1 Attraction
                                        double dist_TM_CD2 = sqrt(pow(double(x-i), 2.0)+pow(double(y-j), 2.0));
                                        if ((dist_TM_CD2 > 0) && (dist_TM_CD2 < Att_CD2TM)){
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Watt_CD2TM));}
                                    }//*/
                                    /*if (Tcell.Pick(i, j) == TM){
                                        double dist_TMTM = sqrt(pow(double(x-i), 2.0)+pow(double(y-j), 2.0));
                                        if ((dist_TMTM > 0) && (dist_TMTM < Latt)){ // IN CASE YOU PICK YOUR OWN POSITION
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Watt));}
                                    }//*/
                                    /*if (Tcell.Pick(i, j) == CD45){
                                        double dist_CD45TM = sqrt(pow(double(x-i), 2.0)+pow(double(y-j), 2.0));
                                        if ((dist_CD45TM < Att_by_TM)){
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Watt_by_TM));} // what weight should i add?
                                    }//*/
                                    if (Tcell.Pick(i, j) == CD45){
                                        double dist_CD45TM = sqrt(pow(double(x-i), 2.0)+pow(double(y-j), 2.0));
                                        if ((dist_CD45TM < Rep_by_TM)){
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Wrep_by_TM));} // what weight should i add?
                                    }//*/
                                    /*if (Tcell.Pick(i, j) == CD2){
                                        double dist_CD2TM = sqrt(pow(double(x-i), 2.0)+pow(double(y-j), 2.0));
                                        if ((dist_CD2TM > 0) && (dist_CD2TM < Att_CD2_by_TM)){ // IN CASE YOU PICK YOUR OWN POSITION
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Watt_CD2_by_TM));} // what weight should i add?
                                    }//*/
                                    /*if (Tcell.Pick(i, j) == CD2_CD48){ // Repulsion
                                        double dist_CD2CD48TM = sqrt(pow(double(x-i), 2.0)+pow(double(y-j), 2.0));
                                        if (dist_CD2CD48TM < Rep_byTM_CD2CD48){
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Wrep_byTM_CD2CD48));}
                                    }//*/
                                    /*if (Tcell.Pick(i, j) == CD2_CD48){ // Attraction
                                        double dist_CD2CD48TM = sqrt(pow(double(x-i), 2.0)+pow(double(y-j), 2.0));
                                        if (dist_CD2CD48TM < Att_byTM_CD2CD48){
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Watt_byTM_CD2CD48));}
                                    }*/
                                    /*if (Tcell.Pick(i, j) == CD28_CD80){
                                        double dist_CD28CD80TM = sqrt(pow(double(x-i), 2.0)+pow(double(y-j), 2.0));
                                        if (dist_CD28CD80TM < Rep_byTM_CD28CD80){
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Wrep_byTM_CD28CD80));}
                                    }//*/
                                }
                            }
//                            float xc = dim/2; //x center
//                            float yc = dim/2; //y center
//                            float cSMAC_radius = (1.5)/a;
//                            double dist_TM_center = sqrt(pow(double(x-xc), 2.0)+pow(double(y-yc), 2.0));
//                            if (chronos>10.){
////                                std::cerr<< "dist_TM_center = " << dist_TM_center << '\t' << "cSMAC_radius = " << cSMAC_radius << std::endl;
//                                if (dist_TM_center > cSMAC_radius){
//                                    std::cerr<< "dist_TM_center - " << dist_TM_center << '\t' << "cSMAC_radius = " << cSMAC_radius << std::endl;
////                                std::cerr<< "I enter in here?"<<std::endl;
//                                    CentrVecTM=-1.0;
//                                }
//                            }
//                            if (chronos>1200.){
//                                    CentrVecTM=0.5;
//                            }
                            Sum.Add(Vector2(x, y, dim/2, dim/2).getNormalize().getMultiplied(CentrVecTM)); // 4
                            double RandRealNumber = RandReal();
                            static int cpWrong = 0;
                            static int cpTotal = 0;
                            static double cpsum = 0;
                            static double cpMax = 0.;

                            cpTotal++;
                            cpsum += (P_diag_TM + P_near_TM)*Sum.Length();
                            cpMax = std::max(cpMax, ((P_diag_TM + P_near_TM)*Sum.Length()));
                            if(((P_diag_TM + P_near_TM)*Sum.Length()) > 1.0) {
                                cpWrong++;
                                //std::cerr << (double) cpWrong / (double) cpTotal << " tot=" << cpTotal << " Max=" << cpMax << " mean=" << cpsum / cpTotal << std::endl;
                            }



                            static int printCpt = 0;
                            printCpt++;
                            //if ((std::fmod(chronos, 60.))<=T) std::cerr << "                                                    t=" << T << " Length " << Sum.Length() << std::endl;


                            //std::cerr<<(P_diag_TM + P_near_TM)*Sum.Length()<<std::endl;
                            if ((Sum.Length() > 1e-12) && (RandRealNumber <= (P_diag_TM + P_near_TM)*(a/a0)*Sum.Length())){
                                int delta = 0;
                                int epsilon = 0;
                                Tcell.getDelta_Epsilon_from_Angle(delta, epsilon, Sum);
                                if ( (Tcell.Pick(x+delta, y+epsilon) == nothing) && (APC.Pick(x+delta, y+epsilon) == nothing) ){
                                    Tcell.Move(x, y, x+delta, y+epsilon);
                                    APC.Move(x, y, x+delta, y+epsilon);
                                }
                                else {
                                    Tcell.Unmoved.insert(std::pair<gridpoint*, gridpoint*>(Tcell.access(x,y), Tcell.access(x+delta, y+epsilon)));
                                    APC.Unmoved.insert(std::pair<gridpoint*, gridpoint*>(APC.access(x,y), APC.access(x+delta, y+epsilon)));
                                }
                            //}//*/
                        }
                    } // end else if (comlpex_forces)*/
                    break;}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                case CD2_CD48:{
                    int ACTION = RandInteraction();
                    if (ACTION == binding_kinetics){
                        if( RandReal() < Poff_CD2_CD48){ // Unbind
                            Tcell.StateChange(CD2,x,y);
                            APC.StateChange(CD48,x,y);
                        }
                    }
                    else if (ACTION == diffusion){
//                        double T_SameNeighbor = 0.0;
//                        double f_Nnn = 1.0;
//                        int AdhesionSites = 1;
//                        int SameNeighbors = 4;
//                        for (int i=-AdhesionSites; i<=AdhesionSites; ++i){
//                            for (int j=-AdhesionSites; j<=AdhesionSites; ++j){
//                                if((Tcell.Pick(x+i, y+j) == CD2_CD48) ){
//                                    ++T_SameNeighbor;
//                                }
//                            }
//                        }
//-------------------------------------------------------------------------------------------------------
//                        if (T_SameNeighbor<=SameNeighbors){
//                            f_Nnn = 1./((1.+((double)T_SameNeighbor))*(1.+((double)T_SameNeighbor))*(1.+((double)T_SameNeighbor))*(1.+((double)T_SameNeighbor)));
//                        }
//                        else if (T_SameNeighbor>SameNeighbors){
//                            f_Nnn = 0.;
//                        }
//-------------------------------------------------------------------------------------------------------
//                        if (T_SameNeighbor<SameNeighbors){
//                            f_Nnn = 1.;
//                        }
//                        else if (T_SameNeighbor>=SameNeighbors){
//                            f_Nnn = 1./(((double)T_SameNeighbor-2)*((double)T_SameNeighbor-2));
//                        }
//-------------------------------------------------------------------------------------------------------
//                        f_Nnn = 1./(1+(double)T_SameNeighbor*T_SameNeighbor);
//                        std::cerr<< "T_SameNeighbor" << '\t' << T_SameNeighbor << '\t' << "f_Nnn"<< '\t' << f_Nnn
//                                 << '\t' << "P" << '\t' <<(P_near_c + P_diag_c)*f_Nnn <<std::endl;
                        double RR = RandReal();
                        if ((RR <= (P_near_c + P_diag_c)/**f_Nnn)*/)){ // P_near_c*f_Nnn -> adhesion
//                            if ((P_near_c + P_diag_c)*f_Nnn < 0.6){
//                            std::cerr<< "RR" << '\t' << RR << '\t' <<"T_SameNeighbor" << '\t' << T_SameNeighbor << '\t' << "f_Nnn"<< '\t' << f_Nnn
//                                     << '\t' << "P" << '\t' <<(P_near_c + P_diag_c)*f_Nnn <<std::endl;
//                            }
//                            if ((P_near_c + P_diag_c)*f_Nnn < 0.6){
//                                std::cerr<< "P = " << (P_near_c + P_diag_c) << ", P*f_Nnn = " << (P_near_c + P_diag_c)*f_Nnn << std::endl;
//                            }
                            int delta = 0;
                            int epsilon = 0;
                            Tcell.getDelta_Epsilon(delta, epsilon);
                            if ((Tcell.Pick(x+delta, y+epsilon) == nothing) && (APC.Pick(x+delta, y+epsilon) == nothing)){
                                Tcell.Move(x, y,x+delta, y+epsilon);
                                APC.Move(x, y,x+delta, y+epsilon);
                            }
                            else {
                                Tcell.Unmoved.insert(std::pair<gridpoint*, gridpoint*>(Tcell.access(x,y), Tcell.access(x+delta, y+epsilon)));
                                APC.Unmoved.insert(std::pair<gridpoint*, gridpoint*>(APC.access(x,y), APC.access(x+delta, y+epsilon)));
                            }
                        }
                    }
                    else if (ACTION == complex_forces){
//                        double T_SameNeighbor = 0.0;
//                        double f_Nnn = 0.0;
//                        int AdhesionSites = 1;
//                        for (int i=-AdhesionSites; i<=AdhesionSites; ++i){
//                            for (int j=-AdhesionSites; j<=AdhesionSites; ++j){
//                                if((Tcell.Pick(x+i, y+j) == CD2_CD48) ){
//                                    ++T_SameNeighbor;
//                                }
//                            }
//                        }
//                        if (T_SameNeighbor<3){
//                            f_Nnn = 1.;
//                        }
//                        else if (T_SameNeighbor>=3){
//                            f_Nnn = 1./T_SameNeighbor*T_SameNeighbor*T_SameNeighbor*T_SameNeighbor*T_SameNeighbor*T_SameNeighbor;
//                        }
//                        std::cerr<< "T_SameNeighbor" << '\t' << T_SameNeighbor << '\t' << "f_Nnn"<< '\t' << f_Nnn <<std::endl;
                            Vector2 Sum = Vector2(0,0);
                            for (int i=std::max(0, x-(int(Lrep)+1)); i<std::min(dim, x+(int(Lrep)+1)); ++i){
                                for (int j=std::max(0, y-(int(Lrep)+1)); j< std::min(dim, y+(int(Lrep)+1)); ++j){
                                    if (Tcell.Pick(i, j) == LI){
                                        double dist_CD2CD48LI = sqrt(pow(double(x-i), 2.0)+pow(double(y-j), 2.0));
                                        if (dist_CD2CD48LI < Lrep){
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Wrep));}
                                    }//*/
                                    /*if (Tcell.Pick(i, j) == CD2_CD48){ // Self Attraction
                                        double dist_CD2CD48 = sqrt(pow(double(x-i), 2.0)+pow(double(y-j), 2.0));
                                        if ((dist_CD2CD48 > 0) && (dist_CD2CD48 < Att_CD2CD48)){
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Watt_CD2CD48));}
                                    }//*/
                                    /*if (Tcell.Pick(i, j) == TM){
                                        double dist_TM_CD2 = sqrt(pow(double(x-i), 2.0)+pow(double(y-j), 2.0));
                                        if ((dist_TM_CD2 > 0) && (dist_TM_CD2 < Att_CD2TM)){
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Watt_CD2TM));}
                                    }//*/
                                    /*if (Tcell.Pick(i, j) == CD28_CD80){ // CD28 - CD2 Attraction
                                        double dist_CD28_CD80 = sqrt(pow(double(x-i), 2.0)+pow(double(y-j), 2.0));
                                        if ((dist_CD28_CD80 > 0) && (dist_CD28_CD80 < Att_CD2_CD28)){
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Watt_CD2_CD28));}
                                    }//*/
                                    /*if (Tcell.Pick(i, j) == TM){
                                        double dist_CD2CD48TM = sqrt(pow(double(x-i), 2.0)+pow(double(y-j), 2.0));
                                        if (dist_CD2CD48TM < Rep_byTM_CD2CD48){
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Wrep_byTM_CD2CD48));}
                                    }//*/
                                    if (Tcell.Pick(i, j) == CD45){
                                        double dist_CD45CC = sqrt(pow(double(x-i), 2.0)+pow(double(y-j), 2.0));
                                        if (dist_CD45CC < Rep_by_CC){
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Wrep_by_CC));}
                                    }//*/
                                    /*if (Tcell.Pick(i, j) == TM){
                                        double dist_CD2CD48TM = sqrt(pow(double(x-i), 2.0)+pow(double(y-j), 2.0));
                                        if (dist_CD2CD48TM < Att_byTM_CD2CD48){
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Watt_byTM_CD2CD48));}
                                    }*/
                                }
                            }
                            Sum.Add(Vector2(x, y, dim/2, dim/2).getNormalize().getMultiplied(CentrVecCD2CD58));
                            double RandRealNumber = RandReal();
                            if ((Sum.Length() > 1e-12) && (RandRealNumber <= ((P_diag_c + P_near_c)/**f_Nnn*/)*(a/a0)*Sum.Length()) ){
                                Sum = Sum.getNormalize();
                                int delta = 0;
                                int epsilon = 0;
                                Tcell.getDelta_Epsilon_from_Angle(delta, epsilon, Sum);
                                if ( (Tcell.Pick(x+delta, y+epsilon) == nothing) && (APC.Pick(x+delta, y+epsilon) == nothing) ){
                                    Tcell.Move(x, y, x+delta, y+epsilon);
                                    APC.Move(x, y, x+delta, y+epsilon);
                                }
                                else {
                                    Tcell.Unmoved.insert(std::pair<gridpoint*, gridpoint*>(Tcell.access(x,y), Tcell.access(x+delta, y+epsilon)));
                                    APC.Unmoved.insert(std::pair<gridpoint*, gridpoint*>(APC.access(x,y), APC.access(x+delta, y+epsilon)));
                                }
                            }
                        //}
                    } // end else*/
                    break;}
                case CD28_CD80:{
                    if(APC.Pick(x,y) != CD28_CD80) std::cerr << "No CD28_CD80" << std::endl;
                    double CC_SameNeighbor = 0.0;
                    for (int i=-1; i<=1; ++i){
                        for (int j=-1; j<=1; ++j){
                            if( ((((i!=0) && (j==0)) || ((i==0) &&(j!=0))) && !((i==0) && (j==0)) ) && (Tcell.Pick(x+i, y+j) == CD28_CD80) ){ // (i*j !=0)
                                ++CC_SameNeighbor;
                                //std::cerr<<T_SameNeighbor<<std::endl;
                            }
                        }
                    }
                    double f_Nnn = 1/(1+CC_SameNeighbor);
                    double P_near_CC = P_near_c;//*(f_Nnn);
                    double P_diag_CC = P_near_CC / sqrt(2);
                    int ACTION = RandInteraction();
                    if (ACTION == binding_kinetics){
                        if( RandReal() < Poff_CD28_CD80){ // Unbind
                            Tcell.StateChange(CD28,x,y);
                            APC.StateChange(CD80,x,y);
                        }
                    }
                    else if (ACTION == diffusion){
                        double RR = RandReal();
                        if ((RR <= P_near_CC + P_diag_CC)){ // P_near_c*f_Nnn -> adhesion
                            int delta = 0;
                            int epsilon = 0;
                            Tcell.getDelta_Epsilon(delta, epsilon);
                            if ((Tcell.Pick(x+delta, y+epsilon) == nothing) && (APC.Pick(x+delta, y+epsilon) == nothing)){
                                Tcell.Move(x, y,x+delta, y+epsilon);
                                APC.Move(x, y,x+delta, y+epsilon);
                            }
                            else {
                                Tcell.Unmoved.insert(std::pair<gridpoint*, gridpoint*>(Tcell.access(x,y), Tcell.access(x+delta, y+epsilon)));
                                APC.Unmoved.insert(std::pair<gridpoint*, gridpoint*>(APC.access(x,y), APC.access(x+delta, y+epsilon)));
                            }
                        }
                    }//*/
                    else if (ACTION == complex_forces){
//                        int delta = 0;
//                        int epsilon = 0;
//                        Tcell.getDelta_Epsilon(delta, epsilon);
                        /*if (Tcell.CheckNeighbors(x,y) > 7){
                            break;
                        }
                        else if(Tcell.CheckNeighbors(x,y)<=7){*/
                            Vector2 Sum = Vector2(0,0);
                            //Sum.Add(Vector2(x, y, x+delta, y+epsilon));
                            for (int i=std::max(0, x-(int(Lrep)+1)); i<std::min(dim, x+(int(Lrep)+1)); ++i){
                                for (int j=std::max(0, y-(int(Lrep)+1)); j< std::min(dim, y+(int(Lrep)+1)); ++j){
                                    if (Tcell.Pick(i, j) == LI){
                                        double dist_CD28CD80LI = sqrt(pow(double(x-i), 2.0)+pow(double(y-j), 2.0));
                                        if (dist_CD28CD80LI < Lrep){
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Wrep));}
                                    }//*/
                                    /*if (Tcell.Pick(i, j) == CD2_CD48){ // Self Attraction
                                        double dist_CD28_CD80 = sqrt(pow(double(x-i), 2.0)+pow(double(y-j), 2.0));
                                        if ((dist_CD28_CD80 > 0) && (dist_CD28_CD80 < Att_CD2_CD28)){
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Watt_CD2_CD28));}
                                    }//*/
                                    /*if (Tcell.Pick(i, j) == TM){
                                        double dist_CD28CD80TM = sqrt(pow(double(x-i), 2.0)+pow(double(y-j), 2.0));
                                        if (dist_CD28CD80TM < Rep_byTM_CD28CD80){
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Wrep_byTM_CD28CD80));}
                                    }//*/
                                    /*if (Tcell.Pick(i, j) == CD28_CD80){
                                        double dist_CD28CD80 = sqrt(pow(double(x-i), 2.0)+pow(double(y-j), 2.0));
                                        if ((dist_CD28CD80 > 0) && (dist_CD28CD80 < Att_CD28CD80)){
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Watt_CD28CD80));}
                                    }//*/
                                }
                            }
                            Sum.Add(Vector2(x, y, dim/2, dim/2).getNormalize().getMultiplied(CentrVecCD28CD80));

                            static int printCptCD = 0;
                            printCptCD++;
                            //if ((std::fmod(chronos, 60.))<=T) std::cerr << "t=" << T << " CD28 Length " << Sum.Length() << std::endl;

                            double RandRealNumber = RandReal();
                            if ((Sum.Length() > 1e-12) && (RandRealNumber <= (P_diag_CC + P_near_CC)*(a/a0)*Sum.Length()) ){
                                Sum = Sum.getNormalize();
                                int delta = 0;
                                int epsilon = 0;
                                Tcell.getDelta_Epsilon_from_Angle(delta, epsilon, Sum);
                                if ( (Tcell.Pick(x+delta, y+epsilon) == nothing) && (APC.Pick(x+delta, y+epsilon) == nothing) ){
                                    Tcell.Move(x, y, x+delta, y+epsilon);
                                    APC.Move(x, y, x+delta, y+epsilon);
                                }
                                else {
                                    Tcell.Unmoved.insert(std::pair<gridpoint*, gridpoint*>(Tcell.access(x,y), Tcell.access(x+delta, y+epsilon)));
                                    APC.Unmoved.insert(std::pair<gridpoint*, gridpoint*>(APC.access(x,y), APC.access(x+delta, y+epsilon)));
                                }
                            }
                        //}
                    } // end else*/
                    break;}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                case LI:{
                    double T_SameNeighbor = 0.0;
                    for (int i=-1; i<=1; ++i){
                        for (int j=-1; j<=1; ++j){
                            if( ((((i!=0) && (j==0)) || ((i==0) &&(j!=0))) && !((i==0) && (j==0)) ) && (Tcell.Pick(x+i, y+j) == LI) ){ // (i*j !=0)
                                ++T_SameNeighbor;
                                //std::cerr<<T_SameNeighbor<<std::endl;
                            }
                        }
                    }
                    double P_near_LI = ((4*sqrt(2)*T*Dc) / (pow(a,2)*(sqrt(2)+1)));//(f_Nnn);
                    double P_diag_LI = P_near_LI / sqrt(2); //*/
                    if(APC.Pick(x,y) != LI) std::cerr << "No LI" << std::endl;
                    int ACTION = RandInteraction();
                    if (ACTION == binding_kinetics){
                        if(RandReal() <= Poff_LI){ // Unbind
                            Tcell.StateChange(LFA1,x,y);
                            APC.StateChange(ICAM1,x,y);
                        }
                    }
                    else if (ACTION ==diffusion){
                        double RR = RandReal();
                        if ( (RR <= P_near_LI + P_diag_LI) ){
                            int delta = 0;
                            int epsilon = 0;
                            Tcell.getDelta_Epsilon(delta, epsilon);
                            if ( (Tcell.Pick(x+delta, y+epsilon) == nothing) && (APC.Pick(x+delta, y+epsilon) == nothing) ){
                                Tcell.Move(x, y, x+delta, y+epsilon);
                                APC.Move(x, y, x+delta, y+epsilon);
                            }
                            else {
                                Tcell.Unmoved.insert(std::pair<gridpoint*, gridpoint*>(Tcell.access(x,y), Tcell.access(x+delta, y+epsilon)));
                                APC.Unmoved.insert(std::pair<gridpoint*, gridpoint*>(APC.access(x,y), APC.access(x+delta, y+epsilon)));
                            }
                        }
                    }//*/
                    else if(ACTION == complex_forces){
//                        int delta = 0;
//                        int epsilon = 0;
//                        Tcell.getDelta_Epsilon(delta, epsilon);
                        /*if (Tcell.CheckNeighbors(x,y) > 7){
                            break;
                        }
                        else if(Tcell.CheckNeighbors(x,y)<=7){*/
                            int nbInteractors = 0;
                            Vector2 Sum = Vector2(0,0);
                            //Sum.Add(Vector2(x, y, x+delta, y+epsilon));
                            for (int i=std::max(0, x-(int(Lrep)+1)); i<std::min(dim, x+(int(Lrep)+1)); ++i){
                                for (int j=std::max(0, y-(int(Lrep)+1)); j< std::min(dim, y+(int(Lrep)+1)); ++j){
                                    if (Tcell.Pick(i, j) == TM){
                                        double dist_LITM = sqrt(pow((double) x-i, 2.0)+pow((double) y-j, 2.0)); // it should work the function Distance too
                                        if ((dist_LITM < Lrep)){
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Wrep));
                                            nbInteractors++;}
                                    }//*/
                                    if (Tcell.Pick(i, j) == CD2_CD48){
                                        double dist_LICD2CD48 = sqrt(pow((double) x-i, 2.0)+pow((double) y-j, 2.0)); // it should work the function Distance too
                                        if ((dist_LICD2CD48 < Lrep)){
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Wrep));
                                            nbInteractors++;}
                                    }//
                                    /*if (Tcell.Pick(i, j) == CD28_CD80){
                                        double dist_LICD2CD48 = sqrt(pow((double) x-i, 2.0)+pow((double) y-j, 2.0)); // it should work the function Distance too
                                        if ((dist_LICD2CD48 < Lrep)){
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Wrep));
                                            nbInteractors++;}
                                    }//*/
                                    if (Tcell.Pick(i, j) == CD45){
                                        double dist_CD45LI = sqrt(pow(double(x-i), 2.0)+pow(double(y-j), 2.0));
                                        if (dist_CD45LI < Rep_by_LI){
                                            Sum.Add(Vector2(x, y, i, j).getNormalize().getMultiplied(Wrep_by_LI));}
                                    }//*/
                                }
                            }
//                            if (chronos>1200.){
//                                CentrVecLI=0.0;
//                            }
                            Sum.Add(Vector2(x, y, dim/2, dim/2).getNormalize().getMultiplied(CentrVecLI));
                            double RandRealNumber = RandReal();
                            if ((Sum.Length() > 1e-12) &&  (RandRealNumber <= (P_diag_LI + P_near_LI)*(a/a0)*Sum.Length())){
                                //Sum = Sum.getNormalize();
                                int delta = 0; int epsilon = 0;
                                Tcell.getDelta_Epsilon_from_Angle(delta, epsilon, Sum);
                                if ( (Tcell.Pick(x+delta, y+epsilon) == nothing) && (APC.Pick(x+delta, y+epsilon) == nothing) ){
                                    Tcell.Move(x, y, x+delta, y+epsilon);
                                    APC.Move(x, y, x+delta, y+epsilon);
                                }
                                else {
                                    Tcell.Unmoved.insert(std::pair<gridpoint*, gridpoint*>(Tcell.access(x,y), Tcell.access(x+delta, y+epsilon)));
                                    APC.Unmoved.insert(std::pair<gridpoint*, gridpoint*>(APC.access(x,y), APC.access(x+delta, y+epsilon)));
                                }
                            }
                        //}
                    } //end else if complex_forces*/
                    break;}
#ifdef USE_SHUFFLE
#else
                    case border:{break;}
                default:{break;}
                } // end switch
            }
//-----------From the APC side -------------------------------------------------------------------------------------------------------------
            else{
                int count2 = 20000;
                while((count2 >0) && (APC.Pick(x,y) == nothing) && (APC.Pick(x,y) == border)){ // if the site is empty or border
                    x = RandGen();
                    y = RandGen();
                    count2--;
                }
                switch (APC.Pick(x, y)){
#endif
                case CD48:{
//                    double APC_SameNeighbor = 0.0;
//                    double CoeffBinding = 1.0;
//                    int AdhesionSites = 1;
//                    int SameNeighbors = 2;
//                    for (int i=-AdhesionSites; i<=AdhesionSites; ++i){
//                        for (int j=-AdhesionSites; j<=AdhesionSites; ++j){
//                            if((APC.Pick(x+i, y+j) == CD2_CD48) ){
//                                ++APC_SameNeighbor;
//                            }
//                        }
//                    }
////-------------------------------------------------------------------------------------------------------
//                    if (APC_SameNeighbor>=SameNeighbors){
//                        CoeffBinding = 0.00001*(1.+((double)APC_SameNeighbor))*(1.+((double)APC_SameNeighbor))*(1.+((double)APC_SameNeighbor))*(1.+((double)APC_SameNeighbor));
//                    }
//                    else if (APC_SameNeighbor<SameNeighbors){
//                        CoeffBinding = 1.;
//                    }
//-------------------------------------------------------------------------------------------------------
//                    coeff_CD2_CD48 = 1;
//                    int negativeRegulation = 0;
//                    double Binding = 0.;
//                    int NeighoringSites = 1; // to define the neighboring sites faster for the two 'for' loops/*
//                    for (int i=-NeighoringSites; i<=NeighoringSites;++i){
//                        for (int j=-NeighoringSites; j<=NeighoringSites;++j){
//                            if (Tcell.Pick(x+i,y+j) == CD2_CD48){
//                                coeff_CD2_CD48++;
//                            }
//                            if (Tcell.Pick(x+i,y+j) == LI){
//                                negativeRegulation++;
//                            }
//                        }
//                    }
//                    if(coeff_CD2_CD48 <=3){
//                        Binding = 1.;
//                    }
//                    else if (coeff_CD2_CD48 >3){
//                        Binding = 0.000001*(double)coeff_CD2_CD48;
//                    }
//                    std::cerr<< "coeff_CD2_CD48 = "<< coeff_CD2_CD48<< '\t' << "Binding = " << Binding << std::endl;
//                    BindCoef = 1.;
//                    double BindingProbability = Pon_CD2_CD48*Binding;
                    int ACTION = RandInteraction();
                    if (ACTION == binding_kinetics){
                        if((Tcell.Pick(x, y) == CD2) && (RandReal() < /*BindingProbability*/Pon_CD2_CD48/**CoeffBinding*//*Binding*/)){
                            Tcell.StateChange(CD2_CD48,x,y);
                            APC.StateChange(CD2_CD48,x,y);
                        }
                    }
                    else if (ACTION == diffusion){
                        double RR = RandReal();
                       if ((NoDTDepDiffusion) || (RR<P_near_m+P_diag_m)){
                            int delta = 0;
                            int epsilon = 0;
                            APC.getDelta_Epsilon(delta, epsilon);
                            APC.Move(x, y,x+delta, y+epsilon);
                        }
                    }
                    break;}
                case CD80:{
                    BindCoef = 1.;
                    int ACTION = RandInteraction();
                    if (ACTION == binding_kinetics){
                        if((Tcell.Pick(x, y) == CD28) && (RandReal() < Pon_CD28_CD80*BindCoef)){
                            Tcell.StateChange(CD28_CD80,x,y);
                            APC.StateChange(CD28_CD80,x,y);
                        }
                    }
                    else if (ACTION == diffusion){
                        double RR = RandReal();
                       if ((NoDTDepDiffusion) || (RR<P_near_m+P_diag_m)){
                            int delta = 0;
                            int epsilon = 0;
                            APC.getDelta_Epsilon(delta, epsilon);
                            APC.Move(x, y,x+delta, y+epsilon);
                        }
                    }
                    break;}
                case pMHC:{
                    BindCoef = 1.;
                    int ACTION = RandInteraction();
                    if (ACTION == binding_kinetics){
                        if((Tcell.Pick(x, y) == TCR) && (RandReal() < Pon_TM*BindCoef)){
                            Tcell.StateChange(TM,x,y);
                            APC.StateChange(TM,x,y);
                            //std::cerr<< "pMHC binding" << std::endl;
                        }
                    }
                    else if (ACTION == diffusion){
                        double RR = RandReal();
                        if ((NoDTDepDiffusion) || (RR<P_near_m+P_diag_m)){
                            int delta = 0;
                            int epsilon = 0;
                            APC.getDelta_Epsilon(delta, epsilon);
                            APC.Move(x, y,x+delta, y+epsilon);
                        }
                    }
                    break;}
                case ICAM1:{
                    BindCoef = 1.;
                    int ACTION = RandInteraction();
                    if (ACTION == binding_kinetics){
                        if((Tcell.Pick(x, y) == LFA1) && (RandReal() < Pon_LI*BindCoef)){
                            Tcell.StateChange(LI,x,y);
                            APC.StateChange(LI,x,y);
                            //std::cerr<< "ICAM binding" << std::endl;
                        }
                    }
                    else if (ACTION == diffusion){
                        double RR = RandReal();
                       if ((NoDTDepDiffusion) || (RR<P_near_m+P_diag_m)){
                            int delta = 0;
                            int epsilon = 0;
                            APC.getDelta_Epsilon(delta, epsilon);
                            APC.Move(x, y,x+delta, y+epsilon);
                        }
                    }
                    break;}
                case border:{break;}
                default:{break;}

#ifdef USE_SHUFFLE
            } // end switch
#else


                } // end switch APC
            } // end else
#endif
        } // checks loop ends here !
        //SwapElements(Tcell, APC, Actin);
    } // Time loop end here

//    DensityTCR.close();
//    DensityMHC.close();
//    DensityLFA1.close();
//    DensityCD45.close();
    DensityTM.close();
    DensityLI.close();
    //DensitypMHC.close();
    //DensityICAM1.close();
    DensityCC.close();
    DensityTM_CD2CD58coloc.close();
//    DensityCD28CD80.close();
//    DensityPD1_PDL1.close();

    InitialAmountOfMolecules.close();

    freeVSboundMHC.close();
    freeVSboundLFA.close();
    freeVSboundCD2.close();



    ms.writeFiles();

}

//    Sanity Checks

//    int IDposition(int x, int y){ return dim*x+ y;}
//    int getX(int ID){return ID/dim;}
//    int getY(int ID){return ID - dim*getX(ID);}


//void TesteAll(){
//    std::ofstream writer("Theta_delta_epsilon.txt");
//    float PI = 3.141592654;
//    for (double theta=-PI; theta<PI;theta+=0.3){
//        //std::cerr<<theta<<std::endl;
//        int delta =0;
//        int epsilon =0;

//    if ( (theta >= -PI) && (theta <= -3.*PI/4.) ){
//        if (RandReal() <= 1.){
//            delta = -1.;
//        }
//        if (RandReal() < tan(theta)){
//            epsilon = -1.;
//        }
//    }
//    else if ( (theta > -3.*PI/4.) && (theta <= -PI/2.) ){
//        if (RandReal() < -tan(theta+PI/2.)){
//            delta = -1.;
//        }
//        if (RandReal() <= 1.){
//            epsilon = -1.;
//        }
//    }
//    else if ( (theta > -PI/2.) && (theta <= -PI/4.) ){
//        if (RandReal() < -tan(theta+PI/2.)){
//            delta = 1.;
//        }
//        if (RandReal() <= 1.){
//            epsilon = -1.;
//        }
//    }
//    else if ( (theta > -PI/4.) && (theta <= 0.0) ){
//        if (RandReal() <= 1.){
//            delta = 1.;
//        }
//        if (RandReal() < -tan(theta)){
//            epsilon = -1;
//        }
//    }
//    else if ( (theta > 0.0) && (theta <= PI/4.) ){
//        if (RandReal() <= 1.){
//            delta = 1.;
//        }
//        if (RandReal() < tan(theta)){
//            epsilon = 1.;
//        }
//    }
//    else if ( (theta > PI/4.) && (theta <= PI/2.) ){
//        if (RandReal() < -tan(theta+PI/2.)){
//            delta = 1.;
//        }
//        if (RandReal() <= 1.){
//            epsilon = 1.;
//        }
//    }
//    else if ( (theta > PI/2.) && (theta <= 3*PI/4.) ){
//        if (RandReal() < -tan(theta+PI/2.)){
//            delta = -1.;
//        }
//        if (RandReal() <= 1.){
//            epsilon = 1.;
//        }
//    }
//    else if ( (theta > 3*PI/4) && (theta <= PI) ){
//        if (RandReal() <= 1){
//            delta = -1.;
//        }
//        if (RandReal() < -tan(theta)){
//            epsilon = 1.;
//        }
//    }
//    writer<< theta*180/PI<<'\t'<<delta<<'\t'<<epsilon<<std::endl;
//    //std::cerr<< "theta = "<<theta*180/PI<< "  delta = " <<delta<< "   epslilon  " << epsilon << std::endl;
//    }


//    {Vector2 V2(1,1);
//    V2 = V2.getNormalize();
//    int delta = -10, epsilon = -10;
//    Vector2 V1(0, 0);
//    for(int i = 0; i <= 1000; ++i){
//        lattice::getDelta_Epsilon_from_Angle(delta, epsilon, V2);
//        std::cerr << "Vector : " << V1.X << " " <<V1.Y << "\t Direction  " << delta << "   " << epsilon << std::endl;
//        V1.Add(Vector2(delta,epsilon));
//        //if (delta*epsilon!=0) std::cerr<<"Diag"<<std::endl;
//        //else std::cerr<<"Near"<<std::endl;
//    }
//    std::cerr<< "Before  "<< V1.Length() <<"  "<< V1.X << "  "  << V1.Y  << "  " <<std::endl;
//    V1 = V1.getNormalize();
//    std::cerr<< V1.Length() <<"  "<< V1.X << "  "  << V1.Y  << "  " << "  V2x = " << V2.X << "   V2y = "<<  V2.Y <<std::endl;

//    }


//if(1){
//    for(int i = -10; i <= 10; ++i){
//        Vector2 V1(i, 1);
//        int delta = -10, epsilon = -10;
//        lattice::getDelta_Epsilon_from_Angle(delta, epsilon, V1);
//        std::cerr << "Vector : " << V1.X << " " <<V1.Y << "\t Direction  " << delta << " " << epsilon << std::endl;
//    }
//    for(int i = -10; i <= 10; ++i){
//        Vector2 V1(1, i);
//        int delta = -10, epsilon = -10;
//        lattice::getDelta_Epsilon_from_Angle(delta, epsilon, V1);
//        std::cerr << "Vector : " << V1.X << " " <<V1.Y << "\t Direction  " << delta << " " << epsilon << std::endl;
//    }
//    for(int i = -10; i <= 10; ++i){
//        Vector2 V1(i, -1);
//        int delta = -10, epsilon = -10;
//        lattice::getDelta_Epsilon_from_Angle(delta, epsilon, V1);
//        std::cerr << "Vector : " << V1.X << " " <<V1.Y << "\t Direction  " << delta << " " << epsilon << std::endl;
//    }
//    for(int i = -10; i <= 10; ++i){
//        Vector2 V1(-1, i);
//        int delta = -10, epsilon = -10;
//        lattice::getDelta_Epsilon_from_Angle(delta, epsilon, V1);
//        std::cerr << "Vector : " << V1.X << " " <<V1.Y << "\t Direction  " << delta << " " << epsilon << std::endl;
//    }
//    {
//        Vector2 V1(-1, .65);
//        int delta = -10, epsilon = -10;
//        lattice::getDelta_Epsilon_from_Angle(delta, epsilon, V1);
//        std::cerr << "Vector : " << V1.X << " " <<V1.Y << "\t Direction  " << delta << " " << epsilon << std::endl;
//    }
//    {
//        Vector2 V1(0.0, 0.0);
//        int delta = -10, epsilon = -10;
//        lattice::getDelta_Epsilon_from_Angle(delta, epsilon, V1);
//        std::cerr << "Vector : " << V1.X << " " <<V1.Y << "\t Direction  " << delta << " " << epsilon << std::endl;
//    }

//    {
//        Vector2 V1(-1.0, -.65);
//        int delta = -10, epsilon = -10;
//        lattice::getDelta_Epsilon_from_Angle(delta, epsilon, V1);
//        std::cerr << "Vector : " << V1.X << " " <<V1.Y << "\t Direction  " << delta << " " << epsilon << std::endl;
//    }
//}
//}
