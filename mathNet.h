#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <unistd.h>
#include <list>

#include "params.h"

using namespace std;

extern const int ADD;
extern const int SUBTRACT;
extern const int MULTIPLY;
extern const int POW;
extern const int DIVIDE;
extern const int NUM_FUNCTIONS;

class Genome;
class Population;

double inline Mod(double x, double y)
{
    if (0 == y)
        return x;

    return x - y * floor(x/y);
}


inline double randfloat() {return (rand())/(RAND_MAX+1.0);}

//Returns a normally distributed deviate with 0 mean and unit variance
//Algorithm is from Numerical Recipes in C, Second Edition
inline double gaussrand() {
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;

  if (iset==0) {
    do {
      v1=2.0*((rand())/(RAND_MAX+1.0))-1.0;
      v2=2.0*((rand())/(RAND_MAX+1.0))-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq>=1.0 || rsq==0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  }
  else {
    iset=0;
    return gset;
  }
}

class gridEntry {

    public:

    Genome* myGenome;
    double myDist;

    gridEntry(Genome* genomePtr) {
        myGenome = genomePtr;
    };

};

class TreeID {

    public:

    int islandID;
    int myID;

    TreeID(int isl, int id) {
        islandID = isl;
        myID = id;
    };

    TreeID(TreeID *copyMe) {
        islandID = copyMe->islandID;
        myID = copyMe->myID;
    };

    bool isEqual(TreeID *compareID);
};

class Tree {

    public:

    Tree* leftTree;
    Tree* rightTree;

    TreeID* tree_id;

    int isTerminal;
    int isConstant;
    double myVal;
    int myDepth;
    int numConstants;
    int numTerminals;
    int numNodes;

    Tree() {};

    Tree(bool makeTerminal, int currentDepth, int maxDepth, bool isGrow, int my_num_outputs, int my_isl, int my_id) {
        myDepth = currentDepth;
        tree_id = new TreeID(my_isl, my_id);

        if (currentDepth == 0) {
            numConstants = 0;
            numTerminals = 0;
            numNodes = 0;
        }

        if (makeTerminal) {
            isTerminal = true;
            numTerminals = 1;
            leftTree = NULL;
            rightTree = NULL;
            if (randfloat() < PERCENT_CONSTANT) {
                isConstant = true;
                myVal = (2.0*(double)MAX_CONSTANT*randfloat() - (double)MAX_CONSTANT);
                numConstants = 1;
                numNodes = 1;
            } else {
                isConstant = false;
                myVal = (double)((int)((NUM_INPUTS + my_num_outputs) * randfloat()));
                numConstants = 0;
                numNodes = 1;
            }
        } else {
            isTerminal = false;
            isConstant = false;
            myVal = (int)(NUM_FUNCTIONS * randfloat());

            if (currentDepth == (maxDepth-1)) {
                leftTree = new Tree(true, (myDepth+1), maxDepth, isGrow, my_num_outputs, my_isl, my_id);
                rightTree = new Tree(true, (myDepth+1), maxDepth, isGrow, my_num_outputs, my_isl, my_id);
            } else if (isGrow) {
                bool nextIsTerminal = false;
                if (randfloat() < PERCENT_TERMINAL) {
                    nextIsTerminal = true;
                }
                leftTree = new Tree(nextIsTerminal, (myDepth+1), maxDepth, isGrow, my_num_outputs, my_isl, my_id);
                nextIsTerminal = false;
                if (randfloat() < PERCENT_TERMINAL) {
                    nextIsTerminal = true;
                }
                rightTree = new Tree(nextIsTerminal, (myDepth+1), maxDepth, isGrow, my_num_outputs, my_isl, my_id);
            } else {
                leftTree = new Tree(false, (myDepth+1), maxDepth, isGrow, my_num_outputs, my_isl, my_id);
                rightTree = new Tree(false, (myDepth+1), maxDepth, isGrow, my_num_outputs, my_isl, my_id);
            }
            numConstants = leftTree->numConstants + rightTree->numConstants;
            numTerminals = leftTree->numTerminals + rightTree->numTerminals;
            numNodes = 1 + leftTree->numNodes + rightTree->numNodes;
        }
    }

    Tree (Tree *copyMe) {
        leftTree = NULL;
        rightTree = NULL;
        if (copyMe->leftTree != NULL)
            leftTree = new Tree(copyMe->leftTree);
        if (copyMe->rightTree != NULL)
            rightTree = new Tree(copyMe->rightTree);

        isTerminal = copyMe->isTerminal;
        isConstant = copyMe->isConstant;
        myVal = copyMe->myVal;
        myDepth = copyMe->myDepth;
        numConstants = copyMe->numConstants;
        numTerminals = copyMe->numTerminals;
        numNodes = copyMe->numNodes;
        tree_id = new TreeID(copyMe->tree_id);
    }

    ~Tree() {
        delete tree_id;
        if (!isTerminal) {
            if (leftTree != NULL)
                delete leftTree;
            if (rightTree != NULL)
                delete rightTree;
        }
    }

    void reduceEqn();
    bool isEqual(Tree *compareTree);
    double evalEquation(double currentData[], int current_num_trees);
    void loadTree(FILE *myFile, int islID, int myID);
    Tree* getNonTerminalTree(int &i);
    bool xOverNonTerminalTree(int &i, Tree *xOverTree);
    Tree* getTerminalTree(int &i);
    bool xOverTerminalTree(int &i, Tree *xOverTree);
    Tree* getConstantTree(int &i);
    Tree* getSubTree(int &i);
    bool xOverSubTree(int &i, Tree *xOverTree);
    string produceMyEqn(int total_num_trees);
    void recountNodes(int curDepth);
    string getTreeString();
    void refCheck(list<int> &treeList, int treeNum, TreeID* treeNumID, int &myCounter, bool isChecked[]);
    void swapRefs(const list<int> &oldIDs, const list<int> &newIDs);
    void updateID(int newIslandID, int newID);
};

class Genome {

    public:

    list<Tree*> trees; //The equation networks of this genome
    list<TreeID*> tree_id_lookup; //ID info for equation networks

    int pop_id;
    int num_trees;
    int fitness;
    int numOffspringProduced;
    double age;
    int mainGridLocation[2];
    Genome* closestInd;
    bool is_alive;
    double currentComm;
    double myPosition[2];
    double *myData, *initialVal;
    bool isNan;
    int distID;

    Genome() {}

    //Generate numInputs initial networks
    Genome(int myID, int maxDepth, bool isGrow, int isl_id) {
        Tree *new_tree;
        int i;
        pop_id = myID;
        num_trees = 0;
        fitness = 0;
        age = 0.0;
        isNan = false;

        is_alive = true;

        generateRandomPosition();
        mainGridLocation[0] = (int)((myPosition[0] + ((double)ISLAND_SIZE / 2.0)) / NN_BLOCK_SIZE);
        mainGridLocation[1] = (int)((myPosition[1] + ((double)ISLAND_SIZE / 2.0)) / NN_BLOCK_SIZE);

        for (i = 0; i < NUM_OUTPUTS; i++) {
            new_tree = new Tree(false, 0, maxDepth, isGrow, NUM_OUTPUTS, isl_id, i);

            trees.push_back(new_tree);
            tree_id_lookup.push_back(new TreeID(new_tree->tree_id));
            num_trees++;
        }

        myData = new double[num_trees];
        initialVal = new double[num_trees];
        for (i = 0; i < NUM_OUTPUTS; i++) {
            initialVal[i] = (2.0 * randfloat()) - 1.0;
        }
    }

    //Constructor that copies the contents of another Genome
    Genome(Genome *copyMe) {
        int i;
        list<Tree*>::iterator curTree;
        Tree *new_tree;
        pop_id = copyMe->pop_id;
        num_trees = copyMe->num_trees;
        myData = new double[num_trees];
        initialVal = new double[num_trees];
        fitness = copyMe->fitness;
        closestInd = copyMe->closestInd;
        is_alive = copyMe->is_alive;
        age = copyMe->age;
        mainGridLocation[0] = copyMe->mainGridLocation[0];
        mainGridLocation[1] = copyMe->mainGridLocation[1];
        isNan = copyMe->isNan;

        for (i = 0; i < 2; i++) {
            myPosition[i] = copyMe->myPosition[i];
        }

        i = 0;
        for (curTree = copyMe->trees.begin(); curTree != copyMe->trees.end(); ++curTree) {
//            cout << "Copying eq net...";
            initialVal[i] = copyMe->initialVal[i];
            //myData[i] = copyMe->myData[i];
            new_tree = new Tree((*curTree));
            trees.push_back(new_tree);
            tree_id_lookup.push_back(new TreeID((*curTree)->tree_id));
//            cout << "done copying eq net." << endl;
            i++;
        }
    }

    // destructor
    ~Genome() {
        list<Tree*>::iterator curTree;
        list<TreeID*>::iterator curTreeID;

        delete [] myData;
        delete [] initialVal;

        for (curTree = trees.begin(); curTree != trees.end(); ++curTree)
            delete (*curTree);
        for (curTreeID = tree_id_lookup.begin(); curTreeID != tree_id_lookup.end(); ++curTreeID)
            delete (*curTreeID);

        trees.clear();
        tree_id_lookup.clear();
    }

    void save_EMM_txt(int myNum);
    void evalEMM(double currentData[]);
    int numRefsToTree(int treeNum);
    void updateIDrefs(const list<int> &oldIDs, const list<int> &newIDs);
    Genome* removeUnusedEquations();
    int countTotalNumNodes();
    Genome* reproduce(Genome* parent2, int my_isl_id, int &max_tree_id, double beta);
    void setNearestNeighbor(list<gridEntry*> *fullGrid, bool &isWithinReproductionDistance, Population *myPop);
    void setNearestNeighbor(list<gridEntry*> *fullGrid, Population *myPop);
    void generateRandomPosition();
    void initializeGenome();
    void addToGrid(list<gridEntry*> *fullGrid);
    void removeFromGrid(list<gridEntry*> *fullGrid);
    void updateOnGrid(list<gridEntry*> *fullGrid, int oldGridPosition[]);
};

class Population {

    public:

    list<Genome*> genomes; //The genomes (organisms) in the population
    int numInputs;
    int pop_size;
    int gen;
    int my_isl_id;
    int max_tree_id;
    double bestFitness;
    double meanFitness;
    double bestFitSoFar;
    double beta;

    Population() {
        numInputs = NUM_INPUTS;
        pop_size = 0;
        bestFitSoFar = DBL_MAX;
        gen = 0;
        my_isl_id = -1;
        max_tree_id = -1;
        beta = 1.0;
    }

    Population(int isl_id) {
        numInputs = NUM_INPUTS;
        pop_size = ISLAND_POP;
        bestFitSoFar = DBL_MAX;
        gen = 0;
        my_isl_id = isl_id;
        max_tree_id = NUM_OUTPUTS - 1;
        beta = 1.0;

        int counter = 0;
        Genome *new_genome;

        for (int i = 0; i < pop_size; i++) {
            // ramped half-and-half (depth 1 to 2, half grow half full)
            if (i < (pop_size/4)) {
                new_genome = new Genome(i, 1, false, isl_id);
            } else if ((i >= (pop_size/4)) && (i < (pop_size/2))) {
                new_genome = new Genome(i, 2, false, isl_id);
            } else if ((i >= (pop_size/2)) && (i < (3*pop_size/4))) {
                new_genome = new Genome(i, 1, true, isl_id);
            } else {
                new_genome = new Genome(i, 2, true, isl_id);
            }
            new_genome->pop_id = counter;
            counter++;

            genomes.push_back(new_genome);
        }
    }

    // destructor
    ~Population() {
        list<Genome*>::iterator curGenome;
        for (curGenome = genomes.begin(); curGenome != genomes.end(); ++curGenome)
            delete (*curGenome);

        genomes.clear();
    }

    //Constructor to create a blank population
//    Population() {}

    void savePop(string myFilename);
    void loadPop(string myFilename);
    void printStats();
    void saveFitness(string myFilename);
    void updateDensityGrid(int densityGrid[], int currentDensityTotals[], int &currDensityGridPointer, int &prevDensityGridPointer);
};
