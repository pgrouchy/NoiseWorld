#include <sstream>
#include <limits>

#include "mathNet.h"

const int ADD = 0;
const int SUBTRACT = 1;
const int MULTIPLY = 2;
const int DIVIDE = 3;
const int NUM_FUNCTIONS = 4;

// Method to save an EMM as a human readable text file
void Genome::save_EMM_txt(int myNum) {
    list<Tree*>::iterator curTree;

    int counter;
    ofstream myFile;
    stringstream tempNum;
    tempNum << myNum;
    string myName = "myEMM_ind" + tempNum.str();
    string myFilename = myName + ".emm";

    do {
        myFile.open(myFilename.c_str(), ios::out | ios::trunc);
        if (!myFile)
            sleep(0.1);
    } while (!myFile);

    counter = 0;
    for (curTree = trees.begin(); curTree != trees.end(); ++curTree) {
        if (counter == 0) {
            myFile << "theta^(t+1) = ";
        } else if (counter == 1) {
            myFile << "omega_out^(t+1) = ";
        } else {
            myFile << "v_" << (counter - 1) << "^(t+1) = ";
        }
        myFile << (*curTree)->produceMyEqn(this->num_trees) << "\n";
        counter++;
    }
    myFile << "\n";
    myFile.close();
}

// This method generates a human readable equation for a given tree in an EMM
string Tree::produceMyEqn(int total_num_trees) {
    stringstream ss;
    string myString;

    if (this->isTerminal) {
        if (this->isConstant) {
            ss << myVal;
        } else {
            int counter = (int)myVal;
            if (counter == 0) {
                ss << "x^t";
            } else if (counter == 1) {
                ss << "y^t";
            } else if (counter == 2) {
                ss << "w_in^t";
            } else if (counter == 3) {
                ss << "theta^t";
            } else if (counter == 4) {
                ss << "omega_out^t";
            } else {
                ss << "v_" << (counter - 4) << "^t";
            }
        }
    } else {
        switch ((int)myVal) {
            case ADD: {
                ss << "(" << this->leftTree->produceMyEqn(total_num_trees) << ") + (" << this->rightTree->produceMyEqn(total_num_trees) << ")";
                break;
            }
            case SUBTRACT: {
                ss << "(" << this->leftTree->produceMyEqn(total_num_trees) << ") - (" << this->rightTree->produceMyEqn(total_num_trees) << ")";
                break;
            }
            case MULTIPLY: {
                ss << "(" << this->leftTree->produceMyEqn(total_num_trees) << ") * (" << this->rightTree->produceMyEqn(total_num_trees) << ")";
                break;
            }
            case DIVIDE: {
                ss << "(" << this->leftTree->produceMyEqn(total_num_trees) << ") / (" << this->rightTree->produceMyEqn(total_num_trees) << ")";
                break;
            }
        }
    }

    return ss.str();
}

// This method evaluates an EMM for one timestep using input and previous output data provided in "currentData"
void Genome::evalEMM(double currentData[]) {
    int i;
    list<Tree*>::iterator curTree;

    curTree = this->trees.begin();
    for (i = 0; i < this->num_trees; i++) {
        // Evaluate an equation tree in the EMM
        this->myData[i] = (*curTree)->evalEquation(currentData, this->num_trees);

        if (this->myData[i] > DBL_MAX) {
            this->myData[i] = DBL_MAX;
        } else if (this->myData[i] < -DBL_MAX) {
            this->myData[i] = -DBL_MAX;
        }

        ++curTree;
    }
}

// This method evaluates a given equation tree
double Tree::evalEquation(double currentData[], int current_num_trees) {
    double myResult = 0.0;

    if (this->isTerminal) {
        if (this->isConstant) {
            myResult = this->myVal;
        } else {
            myResult = currentData[(int)this->myVal];
        }
    } else {
        switch ((int)this->myVal) {
            case ADD: {
                myResult = this->leftTree->evalEquation(currentData, current_num_trees) + this->rightTree->evalEquation(currentData, current_num_trees);
                break;
            }
            case SUBTRACT: {
                myResult = this->leftTree->evalEquation(currentData, current_num_trees) - this->rightTree->evalEquation(currentData, current_num_trees);
                break;
            }
            case MULTIPLY: {
                myResult = this->leftTree->evalEquation(currentData, current_num_trees) * this->rightTree->evalEquation(currentData, current_num_trees);
                break;
            }
            case DIVIDE: {
                myResult = this->leftTree->evalEquation(currentData, current_num_trees) / this->rightTree->evalEquation(currentData, current_num_trees);
                break;
            }
        }
    }

    return myResult;
}

// This method performs some simple (recursive) equation reductions
void Tree::reduceEqn() {
    if (!this->isTerminal) {
        switch ((int)this->myVal) {
            case ADD: {
                this->leftTree->reduceEqn();
                this->rightTree->reduceEqn();
                if ((this->leftTree->isConstant) && (this->rightTree->isConstant)) { // reduce the addition of two constants to a single constant
                    this->myVal = this->leftTree->myVal + this->rightTree->myVal;
                    this->isConstant = true;
                    this->isTerminal = true;
                    delete this->leftTree;
                    delete this->rightTree;
                    this->leftTree = NULL;
                    this->rightTree = NULL;
                } else if (this->leftTree->isEqual(this->rightTree)) { // reduce the addition of two equal trees to 2x a single copy of that tree
                    this->myVal = MULTIPLY;
                    this->rightTree->myVal = 2.0;
                    this->rightTree->isConstant = true;
                    this->rightTree->isTerminal = true;
                    delete this->rightTree->leftTree;
                    this->rightTree->leftTree = NULL;
                    delete this->rightTree->rightTree;
                    this->rightTree->rightTree = NULL;
                }
                break;
            }
            case SUBTRACT: {
                this->leftTree->reduceEqn();
                this->rightTree->reduceEqn();
                if ((this->leftTree->isConstant) && (this->rightTree->isConstant)) { // reduce the subtraction of two constants to a single constant
                    this->myVal = this->leftTree->myVal - this->rightTree->myVal;
                    this->isConstant = true;
                    this->isTerminal = true;
                    delete this->leftTree;
                    delete this->rightTree;
                    this->leftTree = NULL;
                    this->rightTree = NULL;
                } else if (this->leftTree->isEqual(this->rightTree)) { // reduce the subtraction of two equal trees to 0
                    this->myVal = 0.0;
                    this->isConstant = true;
                    this->isTerminal = true;
                    delete this->leftTree;
                    delete this->rightTree;
                    this->leftTree = NULL;
                    this->rightTree = NULL;
                }
                break;
            }
            case MULTIPLY: {
                this->leftTree->reduceEqn();
                this->rightTree->reduceEqn();
                if ((this->leftTree->isConstant) && (this->rightTree->isConstant)) { // reduce the multiplication of two constants to a single constant
                    this->myVal = this->leftTree->myVal * this->rightTree->myVal;
                    this->isConstant = true;
                    this->isTerminal = true;
                    delete this->leftTree;
                    delete this->rightTree;
                    this->leftTree = NULL;
                    this->rightTree = NULL;
                } else if ((this->leftTree->isConstant) && ((int)this->leftTree->myVal == 0)) { // reduce the multiplication of a subtree by 0 to 0 (i.e., delete the subtree)
                    this->myVal = 0.0;
                    this->isConstant = true;
                    this->isTerminal = true;
                    delete this->leftTree;
                    delete this->rightTree;
                    this->leftTree = NULL;
                    this->rightTree = NULL;
                } else if ((this->rightTree->isConstant) && ((int)this->rightTree->myVal == 0)) { // reduce the multiplication of a subtree by 0 to 0 (i.e., delete the subtree)
                    this->myVal = 0.0;
                    this->isConstant = true;
                    this->isTerminal = true;
                    delete this->leftTree;
                    delete this->rightTree;
                    this->leftTree = NULL;
                    this->rightTree = NULL;
                }
                break;
            }
            case DIVIDE: {
                this->leftTree->reduceEqn();
                this->rightTree->reduceEqn();
                if ((this->leftTree->isConstant) && (this->rightTree->isConstant)) { // reduce the division of two constants to a single constant
                    if ((int)this->rightTree->myVal != 0) { // skip any division-by-zero operations, such agents will be selected to die after initial evaluation of EMM (see island.cpp)
                        this->myVal = this->leftTree->myVal / this->rightTree->myVal;
                        this->isConstant = true;
                        this->isTerminal = true;
                        delete this->leftTree;
                        delete this->rightTree;
                        this->leftTree = NULL;
                        this->rightTree = NULL;
                    }
                } else if ((this->leftTree->isConstant) && ((int)this->leftTree->myVal == 0)) { // reduce the division of 0 by a subtree to 0
                    this->myVal = 0.0;
                    this->isConstant = true;
                    this->isTerminal = true;
                    delete this->leftTree;
                    delete this->rightTree;
                    this->leftTree = NULL;
                    this->rightTree = NULL;
                }
                break;
            }
        }
    }
}

// Check (recursively) to see if the contents of two equation trees are identical
bool Tree::isEqual(Tree *compareTree) {
    bool treesAreEqual = false;

    if ((!this->isTerminal) && (!compareTree->isTerminal)) {
        if (this->myVal == compareTree->myVal) {
            if ((this->leftTree->isEqual(compareTree->leftTree)) && (this->rightTree->isEqual(compareTree->rightTree))) {
                treesAreEqual = true;
            }
        }
    } else if ((this->isTerminal) && (compareTree->isTerminal)) {
        if (this->isConstant == compareTree->isConstant) {
            if (this->myVal == compareTree->myVal) {
                treesAreEqual = true;
            }
        }
    }

    return treesAreEqual;
}

// Save a population snapshot to disk
void Population::savePop(string myFilename) {
    ofstream myFile;
    int i;
    string treeString;
    list<Genome*>::iterator curGenome;
    list<Tree*>::iterator curTree;

    do {
        myFile.open(myFilename.c_str(), ios::out | ios::trunc);
        if (!myFile)
            sleep(0.1);
    } while (!myFile);

    // Population data
    myFile << "Timestep: " << this->gen << endl;
    myFile << "Isl: " << this->my_isl_id << endl;
    myFile << "MaxID: " << this->max_tree_id << endl;
    myFile << showpoint << setprecision(numeric_limits<double>::digits10 + 2) << "beta: " << this->beta << endl;

    // Agent & agent genome data
    for (curGenome = this->genomes.begin(); curGenome != this->genomes.end(); ++curGenome) {
        myFile << "Pop: " << (*curGenome)->pop_id << endl;
        if ((*curGenome)->is_alive) {
            myFile << "Alive: 1" << endl;
        } else {
            myFile << "Alive: 0" << endl;
        }
        myFile << "NumTrees: " << (*curGenome)->num_trees << endl;
        myFile << "Fit: " << (*curGenome)->fitness << endl;
        myFile << "Age: " << (*curGenome)->age << endl;
        myFile << "Grid: " << (*curGenome)->mainGridLocation[0] << "," << (*curGenome)->mainGridLocation[1] << endl;
        myFile << showpoint << setprecision(numeric_limits<double>::digits10 + 2) << "InitVals: " << (*curGenome)->initialVal[0];
        for (i = 1; i < (*curGenome)->num_trees; i++) {
            myFile << showpoint << setprecision(numeric_limits<double>::digits10 + 2) << "," << (*curGenome)->initialVal[i];
        }
        myFile << endl;
        myFile << showpoint << setprecision(numeric_limits<double>::digits10 + 2) << "Pos: " << (*curGenome)->myPosition[0] << "," << (*curGenome)->myPosition[1] << endl;

        // tree data
        for (curTree = (*curGenome)->trees.begin(); curTree != (*curGenome)->trees.end(); ++curTree) {
            myFile << "Tree: " << (*curTree)->tree_id->islandID << "," << (*curTree)->tree_id->myID << endl;
            treeString = (*curTree)->getTreeString();
            myFile << treeString;
        }
        myFile << endl;
    }

    myFile.close();
}

// Convert (recursively) a tree into a string for saving to disk
string Tree::getTreeString() {
    stringstream ss;
    string tempString;
    if (this->isConstant) {
        ss << 1;
    } else {
        ss << 0;
    }
    if (this->isTerminal) {
        ss << "," << 1;
    } else {
        ss << "," << 0;
    }
    ss << "," << this->myDepth << ",";
    ss << showpoint << setprecision(numeric_limits<double>::digits10 + 2) << this->myVal;
    ss << "," << this->numConstants << "," << this->numNodes << "," << this->numTerminals << endl;
    if (this->leftTree != NULL) {
        tempString = this->leftTree->getTreeString();
        ss << tempString;
    }
    if (this->rightTree != NULL) {
        tempString = this->rightTree->getTreeString();
        ss << tempString;
    }
    return ss.str();
}

// Load a population snapshot from disk (based on the encoding used in savePop)
void Population::loadPop(string myFilename) {
    Genome *myGenome;
    Tree *myTree;
    FILE *myFile;
    int myInt, myInt2, i;
    double myDouble, myDouble2;

    this->pop_size = 0;

    do {
        myFile = fopen(myFilename.c_str(), "r");
        if (myFile == NULL)
            sleep(0.1);
    } while (myFile == NULL);

    // Population details
    fscanf(myFile, "Timestep: %d\n", &myInt);
    this->gen = myInt;
    fscanf(myFile, "Isl: %d\n", &myInt);
    this->my_isl_id = myInt;
    fscanf(myFile, "MaxID: %d\n", &myInt);
    this->max_tree_id = myInt;
    fscanf(myFile, "beta: %lf\n", &myDouble);
    this->beta = myDouble;

    // Agent and agent genome details
    while (fscanf(myFile, "Pop: %d\n", &myInt) == 1) {
        myGenome = new Genome();
        myGenome->pop_id = myInt;
        fscanf(myFile, "Alive: %d\n", &myInt);
        if (myInt == 1) {
            myGenome->is_alive = true;
        } else {
            myGenome->is_alive = false;
        }
        fscanf(myFile, "NumTrees: %d\n", &myInt);
        myGenome->num_trees = myInt;
        fscanf(myFile, "Fit: %d\n", &myInt);
        myGenome->fitness = myInt;
        fscanf(myFile, "Age: %lf\n", &myDouble);
        myGenome->age = myDouble;
        fscanf(myFile, "Grid: %d,%d\n", &myInt, &myInt2);
        myGenome->mainGridLocation[0] = myInt;
        myGenome->mainGridLocation[1] = myInt2;

        myGenome->initialVal = new double[myGenome->num_trees];

        fscanf(myFile, "InitVals: %lf", &myDouble);
        myGenome->initialVal[0] = myDouble;
        for (i = 1; i < myGenome->num_trees; i++) {
            fscanf(myFile, ",%lf", &myDouble);
            myGenome->initialVal[i] = myDouble;
        }
        fscanf(myFile,"\n");

        myGenome->myData = new double[myGenome->num_trees];

        fscanf(myFile, "Pos: %lf,%lf\n", &myDouble, &myDouble2);
        myGenome->myPosition[0] = myDouble;
        myGenome->myPosition[1] = myDouble2;
        myGenome->mainGridLocation[0] = (int)((myGenome->myPosition[0] + ((double)ISLAND_SIZE / 2.0)) / NN_BLOCK_SIZE);
        myGenome->mainGridLocation[1] = (int)((myGenome->myPosition[1] + ((double)ISLAND_SIZE / 2.0)) / NN_BLOCK_SIZE);

        // Tree details
        for (i = 0; i < myGenome->num_trees; i++) {
            fscanf(myFile, "Tree: %d,%d\n", &myInt, &myInt2);
            myTree = new Tree();
            myTree->loadTree(myFile, myInt, myInt2);
            myTree->recountNodes(0);
            myGenome->trees.push_back(myTree);
            myGenome->tree_id_lookup.push_back(new TreeID(myTree->tree_id));
        }

        fscanf(myFile,"\n");

        this->genomes.push_back(myGenome);
        this->pop_size++;
    }

    fclose(myFile);
}

// Load (recursively) a tree from disk
void Tree::loadTree(FILE *myFile, int islID, int myID) {
    int myInt, myInt2, myInt3, myInt4, myInt5, myInt6;
    double myDouble;

    this->tree_id = new TreeID(islID, myID);

    fscanf(myFile, "%d,%d,%d,%lf,%d,%d,%d\n", &myInt, &myInt2, &myInt3, &myDouble, &myInt4, &myInt5, &myInt6);

    if (myInt == 0)
        this->isConstant = false;
    else
        this->isConstant = true;

    if (myInt2 == 0)
        this->isTerminal = false;
    else
        this->isTerminal = true;

    this->myDepth = myInt3;

    this->myVal = myDouble;

    this->numConstants = myInt4;
    this->numNodes = myInt5;
    this->numTerminals = myInt6;

    if (!this->isTerminal) {
        this->leftTree = new Tree();
        this->leftTree->loadTree(myFile, islID, myID);
        this->rightTree = new Tree();
        this->rightTree->loadTree(myFile, islID, myID);
    } else {
        this->leftTree = NULL;
        this->rightTree = NULL;
    }
}

// Returns the i-th  (counting from 0) nonterminal tree (recursive)
// ASSUMES i < num of nonterminal trees in given subtree
// ASSUMES that the number of nodes in the tree is > 1
Tree* Tree::getNonTerminalTree(int &i) {
    if (i == 0) {
        return this;
    } else {
        Tree *myTree;
        myTree = NULL;
        if (!this->leftTree->isTerminal) {
            i--;
            myTree = this->leftTree->getNonTerminalTree(i);
        }
        if (myTree != NULL) {
            return myTree;
        } else {
            if (!this->rightTree->isTerminal) {
                i--;
                myTree = this->rightTree->getNonTerminalTree(i);
                if (myTree != NULL) {
                    return myTree;
                }
            }
        }
        
        return NULL;
    }
}

// Splices xOverTree on to the i-th (counting from 0) nonterminal tree (recursive)
// ASSUMES i < num of nonterminal trees in given subtree
// ASSUMES that the number of nodes in the tree is > 1
// ASSUMES that the initial call was with i > 0
bool Tree::xOverNonTerminalTree(int &i, Tree *xOverTree) {
    if (i == 0) {
        return true;
    } else {
        bool isMyPointer = false;
        if (!this->leftTree->isTerminal) {
            i--;
            isMyPointer = this->leftTree->xOverNonTerminalTree(i, xOverTree);
        }
        if (isMyPointer) {
            delete this->leftTree;
            this->leftTree = new Tree(xOverTree);
        } else if (this->rightTree != NULL) {
            if (!this->rightTree->isTerminal) {
                i--;
                if (this->rightTree->xOverNonTerminalTree(i, xOverTree)) {
                    delete this->rightTree;
                    this->rightTree = new Tree(xOverTree);
                }
            }
        }
    }
    return false;
}

// Returns the i-th (counting from 0) terminal tree (recursive)
// ASSUMES that the number of nodes in the tree is > 1 (i.e., the intiai node will be a nonterminal node)
Tree* Tree::getTerminalTree(int &i) {
    if ((this->isTerminal) && (i == -1)) {
        return this;
    } else if (!this->isTerminal) {
        Tree *myTree;
        myTree = NULL;
        if (this->leftTree->isTerminal) {
            i--;
        }
        myTree = this->leftTree->getTerminalTree(i);

        if (myTree != NULL) {
            return myTree;
        } else {
            if (this->rightTree->isTerminal) {
                i--;
            }
            myTree =  this->rightTree->getTerminalTree(i);
            if (myTree != NULL) {
                return myTree;
            }
        }
    }
    return NULL;
}

// Splices xOverTree on to the i-th (counting from 0) terminal tree (recursive)
// ASSUMES that the number of nodes in the tree is > 1 (i.e., the intiai node will be a nonterminal node)
bool Tree::xOverTerminalTree(int &i, Tree *xOverTree) {
    if ((this->isTerminal) && (i == -1)) {
        return true;
    } else if (!this->isTerminal) {
        bool isMyPointer = false;
        if (this->leftTree->isTerminal) {
            i--;
        }
        isMyPointer = this->leftTree->xOverTerminalTree(i, xOverTree);

        if (isMyPointer) {
            delete this->leftTree;
            this->leftTree = new Tree(xOverTree);
        } else if (this->rightTree != NULL) {
            if (this->rightTree->isTerminal) {
                i--;
            }
            isMyPointer = this->rightTree->xOverTerminalTree(i, xOverTree);
            if (isMyPointer) {
                delete this->rightTree;
                this->rightTree = new Tree(xOverTree);
            }
        }
    }
    return false;
}

// Returns the i-th (counting from 0) constant node (recursive)
// ASSUMES that the number of nodes in the tree is > 1 (i.e., the intiai node will be a nonterminal node)
// ASSUMES that the number of constants in the tree is > 0
// ASSUMES the i < the numbre of constants in the tree
Tree* Tree::getConstantTree(int &i) {
    if ((this->isConstant) && (i == -1)) {
        return this;
    } else if (!this->isTerminal) {
        Tree *myTree;
        myTree = NULL;
        if (this->leftTree->isConstant) {
            i--;
        }
        myTree = this->leftTree->getConstantTree(i);
        if (myTree != NULL) {
            return myTree;
        } else {
            if (this->rightTree->isConstant) {
                i--;
            }
            myTree = this->rightTree->getConstantTree(i);
            if (myTree != NULL) {
                return myTree;
            }
        }
    }
    return NULL;
}

// Return a pointer to the i-th subtree (counting from 0)
Tree* Tree::getSubTree(int &i) {
    if (i == 0) {
        return this;
    } else if (!this->isTerminal) {
        Tree *myTree = NULL;
        i--;
        myTree = this->leftTree->getSubTree(i);
        if (myTree != NULL) {
            return myTree;
        } else {
            i--;
            myTree = this->rightTree->getSubTree(i);
            if (myTree != NULL) {
                return myTree;
            }
        }
    }
    return NULL;
}

// Splices xOverTree on to the i-th (counting from 0) subtree (recursive)
// ASSUMES that the number of nodes in the tree is > 1 (i.e., the intiai node will be a nonterminal node)
// ASSUMES that i < number of nodes in the tree
// ASSUMES that i > 0
bool Tree::xOverSubTree(int &i, Tree *xOverTree) {
    if (i == 0) {
        return true;
    } else if (!this->isTerminal) {
        bool isMyPointer = false;
        i--;
        isMyPointer = this->leftTree->xOverSubTree(i, xOverTree);
        if (isMyPointer) {
            delete this->leftTree;
            this->leftTree = new Tree(xOverTree);
        } else if (this->rightTree != NULL) {
            i--;
            if (this->rightTree->xOverSubTree(i, xOverTree)) {
                delete this->rightTree;
                this->rightTree = new Tree(xOverTree);
            }
        }
    }
    return false;
}

// Recount the nodes and node types in a tree and update node count values
void Tree::recountNodes(int curDepth) {
    this->myDepth = curDepth;
    curDepth++;
    this->numNodes = 1;
    this->numConstants = 0;
    this->numTerminals = 0;
    if (!this->isTerminal) {
        this->leftTree->recountNodes(curDepth);
        this->numNodes += this->leftTree->numNodes;
        this->numConstants += this->leftTree->numConstants;
        this->numTerminals += this->leftTree->numTerminals;
        if (this->rightTree != NULL) {
            this->rightTree->recountNodes(curDepth);
            this->numNodes += this->rightTree->numNodes;
            this->numConstants += this->rightTree->numConstants;
            this->numTerminals += this->rightTree->numTerminals;
        }
    } else {
        this->numTerminals = 1;
        if (this->isConstant) {
            this->numConstants = 1;
        }
    }
}

// return the number of times a given equation tree (# treeNum in the genome) is referenced by OTHER equation trees in the same genome
//NOTE: This method DOES NOT count self references
int Genome::numRefsToTree(int treeNum) {
    list<Tree*>::iterator curTree;
    list<TreeID*>::iterator curTreeID;
    TreeID *treeNumID;
    list<int> treeList;
    int i, currentTreeNum;
    bool isChecked[this->num_trees];
    int counter = 0;
    curTreeID = this->tree_id_lookup.begin();
    for (i = 0; i < this->num_trees; i++) {
        isChecked[i] = false;
        if (i == treeNum) {
            treeNumID = (*curTreeID);
        }
        ++curTreeID;
    }
    // Only check trees that are directly or indirectly referenced by the output equations
    for (i = 0; i < NUM_OUTPUTS; i++) {
        treeList.push_back(i);
    }
    while (!treeList.empty()) {
        currentTreeNum = treeList.front();
        if (!isChecked[currentTreeNum]) {
            curTree = this->trees.begin();
            for (i = 0; i < currentTreeNum; i++) { // get the pointer to the currentTreeNum tree
                ++curTree;
            }
            (*curTree)->refCheck(treeList, treeNum, treeNumID, counter, isChecked);
            isChecked[currentTreeNum] = true;
        }
        treeList.pop_front();
    }
    return counter;
}

// Add the number of times that treeNumID appears in this tree to the total number of times treeNumID has been seen (myCounter)
// Also adds to the list of other trees to be checked
void Tree::refCheck(list<int> &treeList, int treeNum, TreeID* treeNumID, int &myCounter, bool isChecked[]) {
    int refTreeId;
    if ((this->isTerminal) && (!this->isConstant)) {
        refTreeId = this->myVal - NUM_INPUTS;
        if (refTreeId >= 0) { // terminal is not an input variable
            if (!isChecked[refTreeId]) {
                treeList.push_back(refTreeId);
            }
            if ((refTreeId == treeNum) && (!this->tree_id->isEqual(treeNumID))) {
                myCounter++;
            }
        }
    } else if (!this->isTerminal) {
        this->leftTree->refCheck(treeList, treeNum, treeNumID, myCounter, isChecked);
        this->rightTree->refCheck(treeList, treeNum, treeNumID, myCounter, isChecked);
    }
}

// Changes all references to tree numbers inoldIDs to their corresponding new position in the genome (found in newIDs)
void Genome::updateIDrefs(const list<int> &oldIDs, const list<int> &newIDs) {
    list<Tree*>::iterator curTree;
    for (curTree = this->trees.begin(); curTree != this->trees.end(); ++curTree) { // run through all trees
        (*curTree)->swapRefs(oldIDs, newIDs);
    }
}

// Change all references in a tree from internal variable values listed in oldIDs to their corresponding new values in newIDs (recursive)
void Tree::swapRefs(const list<int> &oldIDs, const list<int> &newIDs) {
    int refTreeId;
    list<int>::const_iterator oldIterator, newIterator;
    if ((this->isTerminal) && (!this->isConstant)) { // this node is an internal variable
        refTreeId = this->myVal - (double)NUM_INPUTS; // get this tree's number
        newIterator = newIDs.begin();
        for (oldIterator = oldIDs.begin(); oldIterator != oldIDs.end(); ++oldIterator) {
            if (refTreeId == (*oldIterator)) { // if this node's internal variable is in the oldIDs list...
                this->myVal = (double)((*newIterator) + NUM_INPUTS); // ...change it to the new value 
                break;
            }
            ++newIterator;
        }
    } else if (!this->isTerminal) {
        this->leftTree->swapRefs(oldIDs, newIDs);
        this->rightTree->swapRefs(oldIDs, newIDs);
    }
}

// Recursively update all nodes' tree IDs with a new ID
void Tree::updateID(int newIslandID, int newID) {
    this->tree_id->islandID = newIslandID;
    this->tree_id->myID = newID;
    if (!this->isTerminal) {
        this->leftTree->updateID(newIslandID, newID);
        this->rightTree->updateID(newIslandID, newID);
    }
}

// Produce a "baby" offspring using genetic material from this genome and parent2 (potentially with mutation) 
Genome* Genome::reproduce(Genome* parent2, int my_isl_id, int &max_tree_id, double beta) {
    list<Genome*>::iterator curGenome;
    list<Tree*>::iterator curTree1, curTree2, curTree3;
    list<TreeID*>::iterator curTreeID1;
    list<TreeID*> commonTreeIDs;
    list<int> oldIDs, newIDs;
    Tree *xOverTree, *tempTree;

    Genome *tempBaby;

    int i, j, counter, counter2, counter3, num_eqns_in_common, tempIslID, tempID, numXOver;
    int xOverPoint1, xOverPoint2;

    bool *isCopied;
    double myChoice, p_m;
    double *tempInitValArray;
    bool isGrow, didXOver;

    // prep parents
    Genome *parent1, *tempParent, *parent1_reorder, *parent2_reorder;

    // prep offspring population
    Genome *baby;

    parent1 = this; // just to make everything a bit clearer

    // re-order parents so that common equations are lined up

    // find and count common trees
    num_eqns_in_common = NUM_OUTPUTS; // always treat the initial NUM_OUTPUTS output equations as the same across all genomes
    commonTreeIDs.clear();
    curTree1 = parent1->trees.begin();
    curTree3 = parent2->trees.begin();
    for (i = 0; i < NUM_OUTPUTS; i++) {
        ++curTree1;
        ++curTree3;
    }

    while (curTree1 != parent1->trees.end()) {
        curTree2 = curTree3;
        while (curTree2 != parent2->trees.end()) {
            if ((*curTree1)->tree_id->isEqual((*curTree2)->tree_id)) { // found an equation in common
                num_eqns_in_common++;
                commonTreeIDs.push_back((*curTree2)->tree_id);
            }
            ++curTree2;
        }
        ++curTree1;
    }

    // re-order parent1
    oldIDs.clear();
    newIDs.clear();
    isCopied = new bool[parent1->num_trees];
    for (i = 0; i < parent1->num_trees; i++) {
        isCopied[i] = false;
    }
    parent1_reorder = new Genome(parent1);

    // delete all tree ids
    for (curTreeID1 = parent1_reorder->tree_id_lookup.begin(); curTreeID1 != parent1_reorder->tree_id_lookup.end(); ++curTreeID1) {
        delete (*curTreeID1);
    }
    parent1_reorder->tree_id_lookup.clear();
    // set tree ids in new order
    curTreeID1 = parent1->tree_id_lookup.begin();
    curTree1 = parent1_reorder->trees.begin();
    for (i = 0; i < NUM_OUTPUTS; i++) { // initial out equations remain unchanged
        parent1_reorder->tree_id_lookup.push_back(new TreeID((*curTreeID1)));
        isCopied[i] = true;
        ++curTreeID1;
        ++curTree1;
    }

    counter = NUM_OUTPUTS;
    for (curTreeID1 = commonTreeIDs.begin(); curTreeID1 != commonTreeIDs.end(); ++curTreeID1) {
        counter2 = 0;
        for (curTree2 = parent1->trees.begin(); curTree2 != parent1->trees.end(); ++curTree2) {
            if ((*curTreeID1)->isEqual((*curTree2)->tree_id)) { // found one of the trees that is in common with parent2
                delete (*curTree1); // delete whichever tree is currently at this spot in the genome
                (*curTree1) = new Tree((*curTree2)); // copy the found common tree to this spot
                parent1_reorder->tree_id_lookup.push_back(new TreeID((*curTreeID1))); // copy the found comm tree's ID to this spot
                ++curTree1; // move to the next tree in the genome
                isCopied[counter2] = true;
                oldIDs.push_back(counter2); // keep track of tree movement so that references internal variables can be updated
                newIDs.push_back(counter); // keep track of tree movement so that references internal variables can be updated
                parent1_reorder->initialVal[counter] = parent1->initialVal[counter2];
                counter++;
            }
            counter2++;
        }
    }

    // now copy the rest of the trees (i.e., the ones with no counterpart in the other parent)
    counter2 = 0;
    for (curTree2 = parent1->trees.begin(); curTree2 != parent1->trees.end(); ++curTree2) {
        if (!isCopied[counter2]) {
            delete (*curTree1);
            (*curTree1) = new Tree((*curTree2));
            parent1_reorder->tree_id_lookup.push_back(new TreeID((*curTree2)->tree_id));
            ++curTree1;
            isCopied[counter2] = true;
            oldIDs.push_back(counter2);
            newIDs.push_back(counter);
            parent1_reorder->initialVal[counter] = parent1->initialVal[counter2];
            counter++;
        }
        counter2++;
    }

    parent1_reorder->updateIDrefs(oldIDs, newIDs); // update references to internal variables, as trees may have been moved around

    // re-order parent2
    // (similar approach as above)
    oldIDs.clear();
    newIDs.clear();
    delete [] isCopied;
    isCopied = new bool[parent2->num_trees];

    for (i = 0; i < parent2->num_trees; i++) {
        isCopied[i] = false;
    }
    parent2_reorder = new Genome(parent2);
    for (curTreeID1 = parent2_reorder->tree_id_lookup.begin(); curTreeID1 != parent2_reorder->tree_id_lookup.end(); ++curTreeID1) {
        delete (*curTreeID1);
    }
    parent2_reorder->tree_id_lookup.clear();
    curTreeID1 = parent2->tree_id_lookup.begin();
    curTree1 = parent2_reorder->trees.begin();
    for (i = 0; i < NUM_OUTPUTS; i++) {
        parent2_reorder->tree_id_lookup.push_back(new TreeID((*curTreeID1)));
        isCopied[i] = true;
        ++curTreeID1;
        ++curTree1;
    }

    counter = NUM_OUTPUTS;
    for (curTreeID1 = commonTreeIDs.begin(); curTreeID1 != commonTreeIDs.end(); ++curTreeID1) {
        counter2 = 0;
        for (curTree2 = parent2->trees.begin(); curTree2 != parent2->trees.end(); ++curTree2) {
            if ((*curTreeID1)->isEqual((*curTree2)->tree_id)) {
                delete (*curTree1);
                (*curTree1) = new Tree((*curTree2));
                parent2_reorder->tree_id_lookup.push_back(new TreeID((*curTreeID1)));
                ++curTree1;
                isCopied[counter2] = true;
                oldIDs.push_back(counter2);
                newIDs.push_back(counter);
                parent2_reorder->initialVal[counter] = parent2->initialVal[counter2];
                counter++;
            }
            counter2++;
        }
    }

    // need to re-id these equations, as they are getting copied to the very end of parent1's equation tree list
    counter = parent1->num_trees;

    counter2 = 0;
    counter3 = num_eqns_in_common;
    for (curTree2 = parent2->trees.begin(); curTree2 != parent2->trees.end(); ++curTree2) {
        if (!isCopied[counter2]) {
            delete (*curTree1);
            (*curTree1) = new Tree((*curTree2));
            parent2_reorder->tree_id_lookup.push_back(new TreeID((*curTree2)->tree_id));
            ++curTree1;
            isCopied[counter2] = true;
            oldIDs.push_back(counter2);
            newIDs.push_back(counter);
            parent2_reorder->initialVal[counter3] = parent2->initialVal[counter2];
            counter++;
            counter3++;
        }
        counter2++;
    }

    parent2_reorder->updateIDrefs(oldIDs, newIDs);
    delete [] isCopied;

    // create baby out of the reordered first parent
    baby = new Genome(parent1_reorder);

    // Copy extra equations of parent2 over to baby, so that they can be referenced by trees from parent2 that might be spliced on to the baby genome
    curTree1 = parent2_reorder->trees.begin();
    curTreeID1 = parent2_reorder->tree_id_lookup.begin();
    for (i = 0; i < num_eqns_in_common; i++) {
        ++curTree1;
        ++curTreeID1;
    }
    while (curTree1 != parent2_reorder->trees.end()) {
        baby->trees.push_back(new Tree((*curTree1)));
        baby->num_trees++;
        baby->tree_id_lookup.push_back(new TreeID((*curTreeID1)));
        ++curTree1;
        ++curTreeID1;
    }

    // recreate the initial value array as the number of trees in the genome may have changed
    delete [] baby->initialVal;
    baby->initialVal = new double[baby->num_trees];
    for (i = 0; i < num_eqns_in_common; i++) {
        baby->initialVal[i] = parent1_reorder->initialVal[i];
    }
    j = num_eqns_in_common;
    for (i = num_eqns_in_common; i < baby->num_trees; i++) {
        if (i < parent1_reorder->num_trees) {
            baby->initialVal[i] = parent1_reorder->initialVal[i];
        } else {
            baby->initialVal[i] = parent2_reorder->initialVal[j];
            j++;
        }
    }

    // myData gets set up properly later, but this needs to be done to avoid seg faults in subsequent operations
    delete [] baby->myData;
    baby->myData = new double[baby->num_trees];

    numXOver = 0;
    tempBaby = new Genome(baby);

    p_m = MUT_RATE * (1.0 / (double)baby->num_trees) * beta;

    while ((numXOver == 0) || (numXOver == num_eqns_in_common)) { // an offspring genome must contain at least one equation from each parent
        numXOver = 0;
        counter = 0;
        curTree2 = parent2_reorder->trees.begin();
        for (curTree1 = baby->trees.begin(); curTree1 != baby->trees.end(); ++curTree1) {
            didXOver = false;
            if ((counter < num_eqns_in_common) && (randfloat() < 0.5)) {
                delete (*curTree1);
                (*curTree1) = new Tree((*curTree2));
                baby->initialVal[counter] = parent2_reorder->initialVal[counter];
                didXOver = true;
                numXOver++;
            }
            if (randfloat() < p_m) { // perform a genetic splice operation
                tempIslID = (*curTree1)->tree_id->islandID;
                tempID = (*curTree1)->tree_id->myID;

                if (!didXOver) { // take subtree from other parent's genome (i.e., not from the parent where this tree originated)
                    j = (int)((parent2->num_trees) * randfloat());
                    i = 0;
                    for (curTree3 = parent2_reorder->trees.begin(); curTree3 != parent2_reorder->trees.end(); ++curTree3) {
                        i++;
                        if (i > j) {
                            break;
                        }
                    }
                } else {
                    j = (int)((parent1_reorder->num_trees) * randfloat());
                    i = 0;
                    for (curTree3 = parent1_reorder->trees.begin(); curTree3 != parent1_reorder->trees.end(); ++curTree3) {
                        i++;
                        if (i > j) {
                            break;
                        }
                    }
                }
                // Get a pointer (xOverTree) to the subtree that will be spliced on to the offspring's current tree (*curTree1)
                if ((randfloat() < PROB_NON_TERM_XOVER) && ((*curTree3)->numNodes != 1)) { // select a nonterminal subtree
                    xOverPoint2 = (int)(((*curTree3)->numNodes - (*curTree3)->numTerminals) * randfloat());
                    xOverTree = (*curTree3)->getNonTerminalTree(xOverPoint2);
                } else { // select a terminal subtree
                    if ((*curTree3)->numNodes == 1) {
                        xOverTree = (*curTree3);
                    } else {
                        xOverPoint2 = (int)((*curTree3)->numTerminals * randfloat());
                        xOverTree = (*curTree3)->getTerminalTree(xOverPoint2);
                    }
                }
                // Splice a copy of the selected subtree on to a spot on the offspring's current tree (*curTree1)
                if ((randfloat() < PROB_NON_TERM_XOVER) && ((*curTree1)->numNodes != 1)) { // select a nonterminal subtree
                    xOverPoint1 = (int)(((*curTree1)->numNodes - (*curTree1)->numTerminals) * randfloat());
                    if (xOverPoint1 == 0) {
                        delete (*curTree1);
                        (*curTree1) = new Tree(xOverTree);
                    } else {
                        (*curTree1)->xOverNonTerminalTree(xOverPoint1, xOverTree);
                    }
                } else { // selected a terminal subtree
                    if ((*curTree1)->numNodes == 1) {
                        delete (*curTree1);
                        (*curTree1) = new Tree(xOverTree);
                    } else {
                        xOverPoint1 = (int)((*curTree1)->numTerminals * randfloat());
                        (*curTree1)->xOverTerminalTree(xOverPoint1, xOverTree);
                    }
                }
                (*curTree1)->updateID(tempIslID, tempID);
                (*curTree1)->recountNodes(0);
            }
            counter++;
            if (counter < num_eqns_in_common) {
                ++curTree2;
            }
        }
        if ((numXOver == 0) || (numXOver == num_eqns_in_common)) { // an offspring genome must contain at least one equation from each parent
            delete baby;
            baby = new Genome(tempBaby); // revert back to initial baby genome and try again
        }
    }

    delete tempBaby;

    delete parent1_reorder;
    delete parent2_reorder;

    // apply mutations
    for (curTree1 = baby->trees.begin(); curTree1 != baby->trees.end(); ++curTree1) {
        if (randfloat() < (0.5 * p_m)) { // add a new equation tree

            baby->num_trees++;
            max_tree_id++;
            // create a new random tree using ramped half and half with a max depth of 1 or 2
            myChoice = randfloat();
            if (myChoice < 0.25) {
                tempTree = new Tree(false, 0, 1, false, baby->num_trees, my_isl_id,  max_tree_id);
            } else if (myChoice < 0.5) {
                tempTree = new Tree(false, 0, 2, false, baby->num_trees, my_isl_id,  max_tree_id);
            } else if (myChoice < 0.75) {
                tempTree = new Tree(false, 0, 1, true, baby->num_trees, my_isl_id,  max_tree_id);
            } else {
                tempTree = new Tree(false, 0, 2, true, baby->num_trees, my_isl_id,  max_tree_id);
            }
            baby->trees.push_back(tempTree);
            baby->tree_id_lookup.push_back(new TreeID(tempTree->tree_id));

            // resize the intial value array and generate a random initial value for this new equation
            tempInitValArray = new double[baby->num_trees - 1];
            for (i = 0; i < (baby->num_trees - 1); i++) {
                tempInitValArray[i] = baby->initialVal[i];
            }
            delete [] baby->initialVal;
            baby->initialVal = new double[baby->num_trees];
            for (i = 0; i < (baby->num_trees - 1); i++) {
                baby->initialVal[i] = tempInitValArray[i];
            }
            baby->initialVal[baby->num_trees - 1] = (2.0 * randfloat()) - 1.0;
            delete [] tempInitValArray;

            // resize myData array, as the number of trees has changed
            delete [] baby->myData;
            baby->myData = new double[baby->num_trees];

            // add a simple connection to this new tree somewhere within the existing trees
            if ((*curTree1)->numNodes == 1) {
                tempTree = (*curTree1);
            } else {
                j = (int)((*curTree1)->numTerminals * randfloat());
                tempTree = (*curTree1)->getTerminalTree(j);
            }
            tempTree->isConstant = false;
            tempTree->myVal = (double)(NUM_INPUTS + baby->num_trees - 1);

            (*curTree1)->recountNodes(0);
        }
    }

    // initial value mutation
    curTree1 = baby->trees.begin();
    for (i = 0; i < baby->num_trees; i++) {
        if (randfloat() < p_m) {
            if (randfloat() < 0.5) {
                baby->initialVal[i] = (2.0 * randfloat()) - 1.0;
            } else {
                baby->initialVal[i] += (MUT_INITIAL_COND_STD * gaussrand());
            }
        }
        ++curTree1;
    }

    // tree mutation
    for (curTree1 = baby->trees.begin(); curTree1 != baby->trees.end(); ++curTree1) {

        if (randfloat() < p_m) { // do a tree mutation
            if (randfloat() < PROB_POINT_MUT) { // do a point mutation
                if (((*curTree1)->numConstants > 0) && (randfloat() < PROB_MUT_CONST)) { // do a perturbation of a constant
                    if ((*curTree1)->isConstant) {
                        xOverTree = (*curTree1);
                    } else {
                        i = (int)((*curTree1)->numConstants * randfloat());
                        xOverTree = (*curTree1)->getConstantTree(i);
                    }
                    xOverTree->myVal += POINT_MUT_STD * gaussrand();
                } else { // mutate a node
                    i = (int)((*curTree1)->numNodes * randfloat());
                    xOverTree = (*curTree1)->getSubTree(i);
                    myChoice = xOverTree->myVal;
                    if (xOverTree->isTerminal) { // picked a terminal node to mutate
                        if (xOverTree->isConstant) { // picked a constant to mutate
                            if (randfloat() < PROB_MUT_TERM_TYPE) { // change from a constant to a variable
                                xOverTree->isConstant = false;
                                myChoice = (int)((NUM_INPUTS + baby->num_trees) * randfloat());
                            } else { // otherwise, set to a new random constant
                                while (myChoice == xOverTree->myVal) {
                                    myChoice = (2.0*(double)MAX_CONSTANT*randfloat() - (double)MAX_CONSTANT);
                                }
                            }
                        } else { // picked a variable to mutate
                            if (randfloat() < PROB_MUT_TERM_TYPE) { // change it to a constant
                                xOverTree->isConstant = true;
                                myChoice = (2.0*(double)MAX_CONSTANT*randfloat() - (double)MAX_CONSTANT);
                            } else { // otherwise change it to a different variable
                                while (myChoice == xOverTree->myVal) {
                                    myChoice = (double)((int)((NUM_INPUTS + baby->num_trees) * randfloat()));
                                }
                            }
                        }
                    } else { // picked a nonterminal node, mutate it to a different nonterminal node
                        while (myChoice == xOverTree->myVal) {
                            myChoice = (double)((int)(NUM_FUNCTIONS * randfloat()));
                        }
                    }
                    xOverTree->myVal = myChoice;
                }
            } else { // subtree mutation

                // generate a new random subtree
                isGrow = false;
                if (randfloat() < 0.5)
                    isGrow = true;
                i = 1;
                if (randfloat() < 0.5)
                    i = 2;

                xOverTree = new Tree(false, 0, i, isGrow, baby->num_trees, (*curTree1)->tree_id->islandID, (*curTree1)->tree_id->myID);

                // tree swap (make the original tree the new random subtree, and then splice the original tree on to the new subtree)
                if (randfloat() < 0.05) {
                    tempTree = xOverTree;
                    xOverTree = (*curTree1);
                    (*curTree1) = tempTree;
                }

                // pick a random node
                i = (int)((*curTree1)->numNodes * randfloat());
                // splice the xOverTree on to that node
                if (i == 0) {
                    delete (*curTree1);
                    (*curTree1) = new Tree(xOverTree);
                } else {
                    (*curTree1)->xOverSubTree(i, xOverTree);
                }

                delete xOverTree;
            }
            (*curTree1)->recountNodes(0);
        }
    }

    // Strip away ALL unused equations
    tempParent = baby->removeUnusedEquations();
    delete baby;
    baby = tempParent;

    // resize myData array, as the number of trees may have changed
    delete [] baby->myData;
    baby->myData = new double[baby->num_trees];

    baby->fitness = 0;
    baby->is_alive = true;
    baby->age = 0.0;
    if (baby->countTotalNumNodes() > MAX_TOTAL_NODES) { // offspring born with more than MAX_TOTAL_NODES will be selected to die at the next island reproduction event
        baby->is_alive = false;
    }
    baby->isNan = false;

    return baby;
}

// Check to see if two tree ids are equal
bool TreeID::isEqual(TreeID *compareID) {
    bool result = false;
    if ((this->myID == compareID->myID) && (this->islandID == compareID->islandID)) {
        result = true;
    }
    return result;
}

// Return a genome that is the same as this one, but with extra equations that are not directly or indirectly referenced by the original NUM_OUTPUTS output equations removed
Genome* Genome::removeUnusedEquations() {
    Genome *tempParent;
    Tree *tempTree;
    list<Tree*>::iterator curTree1;
    list<TreeID*>::iterator curTreeID;
    list<int> oldIDs, newIDs;
    int i, j, numNotConnected = 0;
    bool isConnected[this->num_trees];

    tempParent = new Genome(this);

    for (i = 0; i < NUM_OUTPUTS; i++) { // we will never remove the original output equations
        isConnected[i] = true;
    }
    for (i = NUM_OUTPUTS; i < this->num_trees; i++) {
        isConnected[i] = true;
            
        if (this->numRefsToTree(i) == 0) { // this tree is not directly or indirectly referenced by an original output equation
            numNotConnected++;
            isConnected[i] = false;
        }
    }
    if (numNotConnected > 0) { // remove unconnected trees

        // delete all tree ids, will need to reorder
        curTreeID = tempParent->tree_id_lookup.begin();
        for (curTree1 = tempParent->trees.begin(); curTree1 != tempParent->trees.end(); ++curTree1) {
            delete (*curTree1);
            delete (*curTreeID);
            ++curTreeID;
        }

        // delete all trees, will need to reorder
        tempParent->trees.clear();
        tempParent->tree_id_lookup.clear();
        tempParent->num_trees = 0;
        curTree1 = this->trees.begin();
        curTreeID = this->tree_id_lookup.begin();
        oldIDs.clear();
        newIDs.clear();
        for (i = 0; i < this->num_trees; i++) {
            if (i < NUM_OUTPUTS) { // put original equations back
                tempParent->trees.push_back(new Tree((*curTree1)));
                tempParent->tree_id_lookup.push_back(new TreeID((*curTreeID)));
                tempParent->num_trees++;
            } else if (isConnected[i]) { // put connected equations back (leaving out unconnected equation tree)
                tempTree = new Tree((*curTree1));
                oldIDs.push_back(i); // keep track of order changes
                newIDs.push_back(tempParent->num_trees); // keep track of order changes
                tempParent->trees.push_back(tempTree);
                tempParent->tree_id_lookup.push_back(new TreeID((*curTreeID)));
                tempParent->num_trees++;

            }
            ++curTree1;
            ++curTreeID;
        }

        // resize initial value array as the number of trees has changed
        delete [] tempParent->initialVal;
        tempParent->initialVal = new double[tempParent->num_trees];
        j = 0;
        for (i = 0; i < this->num_trees; i++) {
            if (isConnected[i]) {
                tempParent->initialVal[j] = this->initialVal[i];
                j++;
            }
        }

        // reorder references to internal variables as trees have been reordered
        tempParent->updateIDrefs(oldIDs, newIDs);
    }
    return tempParent;
}

// Return the total number of nodes in a genome
int Genome::countTotalNumNodes() {
    int numNodes = 0;
    list<Tree*>::iterator curTree;
    for (curTree = this->trees.begin(); curTree != this->trees.end(); ++curTree) {
        numNodes += (*curTree)->numNodes;
    }
    return numNodes;
}

// Finds an agent's nearest neighbor by using an expanding search of a discretized position space
// Also sets the boolean isWithinReproductionDistance based on whether the nearest neighbor is within the specified REPRODUCTION_DIST
void Genome::setNearestNeighbor(list<gridEntry*> *fullGrid, bool &isWithinReproductionDistance, Population *myPop) {
    list<Genome*>::iterator curGenome, curGenome2; 
    list<gridEntry*> innerEntriesToCheck;
    list<gridEntry*> outerEntriesToCheck;
    list<gridEntry*>::iterator curGridEntry;
    double minDist, myDist, closestIndDist = DBL_MAX;
    int i, j, xGridBlock, yGridBlock;
    int counter = 1;
    bool neighbourFound = false;

    // first check the grid location that is also occupied by the agent
    for (curGridEntry = fullGrid[this->mainGridLocation[1]*(int)(ISLAND_SIZE/NN_BLOCK_SIZE) + this->mainGridLocation[0]].begin(); curGridEntry != fullGrid[this->mainGridLocation[1]*(int)(ISLAND_SIZE/NN_BLOCK_SIZE) + this->mainGridLocation[0]].end(); ++curGridEntry) {
        if (((*curGridEntry)->myGenome != this)) {
            (*curGridEntry)->myDist = sqrt(pow((this->myPosition[0] - (*curGridEntry)->myGenome->myPosition[0]),2.0) + pow((this->myPosition[1] - (*curGridEntry)->myGenome->myPosition[1]),2.0));
            innerEntriesToCheck.push_back((*curGridEntry));
        }
    }

    do { // check next group of grid sqaures
        for (i = (-1*counter); i <= counter; i++) {
            for (j = (-1*counter); j <= counter; j++) {
                if ((i == (-1*counter)) || (j == (-1*counter)) || (i == counter) || (j == counter)) {
                    xGridBlock = this->mainGridLocation[0] + i;
                    if (xGridBlock >= (int)(ISLAND_SIZE/NN_BLOCK_SIZE)) {
                        continue;
                    } else if (xGridBlock < 0) {
                        continue;
                    }
                    yGridBlock = this->mainGridLocation[1] + j;
                    if (yGridBlock >= (int)(ISLAND_SIZE/NN_BLOCK_SIZE)) {
                        continue;
                    } else if (yGridBlock < 0) {
                        continue;
                    }
                    for (curGridEntry = fullGrid[yGridBlock*(int)(ISLAND_SIZE/NN_BLOCK_SIZE) + xGridBlock].begin(); curGridEntry != fullGrid[yGridBlock*(int)(ISLAND_SIZE/NN_BLOCK_SIZE) + xGridBlock].end(); ++curGridEntry) {
                        (*curGridEntry)->myDist = sqrt(pow((this->myPosition[0] - (*curGridEntry)->myGenome->myPosition[0]),2.0) + pow((this->myPosition[1] - (*curGridEntry)->myGenome->myPosition[1]),2.0));
                        outerEntriesToCheck.push_back((*curGridEntry));
                    }
                }
            }
        }

        // first check the inner grid locations for a neighbor
        for (curGridEntry = innerEntriesToCheck.begin(); curGridEntry != innerEntriesToCheck.end(); ++curGridEntry) {
            if ((*curGridEntry)->myDist < closestIndDist) {
                neighbourFound = true;
                closestIndDist = (*curGridEntry)->myDist;
                this->closestInd = (*curGridEntry)->myGenome;
            }
        }

        if (neighbourFound) { // if we find a neighbor, we still need to search the next layer of grid locations as a closer neighbor may be located there
            for (curGridEntry = outerEntriesToCheck.begin(); curGridEntry != outerEntriesToCheck.end(); ++curGridEntry) {
                if ((*curGridEntry)->myDist < closestIndDist) {
                    closestIndDist = (*curGridEntry)->myDist;
                    this->closestInd = (*curGridEntry)->myGenome;
                }
            }
        } else { // get ready to check the next layer
            innerEntriesToCheck.clear();
            innerEntriesToCheck.splice(innerEntriesToCheck.begin(), outerEntriesToCheck);
            outerEntriesToCheck.clear();
            counter++;
        }
    } while ((!neighbourFound) && (counter <= 2));

    if (!neighbourFound) { // if after 5 layers no neighbor is found, do a search of the entire population
        minDist = DBL_MAX;
        for (curGenome = myPop->genomes.begin(); curGenome != myPop->genomes.end(); ++curGenome) {
            if ((*curGenome) != this) {
                myDist = sqrt(pow((this->myPosition[0] - (*curGenome)->myPosition[0]), 2.0) + pow((this->myPosition[1] - (*curGenome)->myPosition[1]), 2.0));
                if (myDist < minDist) {
                    minDist = myDist;
                    this->closestInd = (*curGenome);
                    closestIndDist = minDist;
                }
            }
        }
    }

    isWithinReproductionDistance = false;
    if (closestIndDist < REPRODUCTION_DIST) {
        isWithinReproductionDistance = true;
    }
}

void Genome::setNearestNeighbor(list<gridEntry*> *fullGrid, Population *myPop) {
    bool temp_bool;
    this->setNearestNeighbor(fullGrid, temp_bool, myPop);
}

// Generate a random position for an agent within the START_RADIUS
void Genome::generateRandomPosition() {
    double r = START_RADIUS * sqrt(randfloat());
    double theta = 2.0 * M_PI * randfloat();
    this->myPosition[0] = r * cos(theta);
    this->myPosition[1] = r * sin(theta);

    this->mainGridLocation[0] = (int)((this->myPosition[0] + ((double)ISLAND_SIZE / 2.0)) / NN_BLOCK_SIZE);
    this->mainGridLocation[1] = (int)((this->myPosition[1] + ((double)ISLAND_SIZE / 2.0)) / NN_BLOCK_SIZE);
}

// Set all outputs and internal variables to their initial values
void Genome::initializeGenome() {
    int i;
    for (i = 0; i < this->num_trees; i++) {
        this->myData[i] = this->initialVal[i];
    }
    this->currentComm = 0.0;
}

// add an agent to the grid used for nearest neighbor search
void Genome::addToGrid(list<gridEntry*> *fullGrid) {
    gridEntry *newGridEntry = new gridEntry(this);
    fullGrid[(this->mainGridLocation[1]*(int)(ISLAND_SIZE/NN_BLOCK_SIZE) + this->mainGridLocation[0])].push_back(newGridEntry);
}

// remove an agent to the grid used for nearest neighbor search
void Genome::removeFromGrid(list<gridEntry*> *fullGrid) {
    list<gridEntry*>::iterator curGridEntry;
    for (curGridEntry = fullGrid[(this->mainGridLocation[1]*(int)(ISLAND_SIZE/NN_BLOCK_SIZE) + this->mainGridLocation[0])].begin(); curGridEntry != fullGrid[(this->mainGridLocation[1]*(int)(ISLAND_SIZE/NN_BLOCK_SIZE) + this->mainGridLocation[0])].end(); ++curGridEntry) {
        if ((*curGridEntry)->myGenome == this) {
            delete (*curGridEntry);
            fullGrid[(this->mainGridLocation[1]*(int)(ISLAND_SIZE/NN_BLOCK_SIZE) + this->mainGridLocation[0])].erase(curGridEntry);
            break;
        }
    }
}

// update an agent's position on the grid used for nearest neighbor search
void Genome::updateOnGrid(list<gridEntry*> *fullGrid, int oldGridPosition[]) {
    // no update needed if genome is still in the same grid square as before
    if ((oldGridPosition[0] != this->mainGridLocation[0]) || (oldGridPosition[1] != this->mainGridLocation[1])) {
        int tempGridPosition[2];
        tempGridPosition[0] = this->mainGridLocation[0];
        tempGridPosition[1] = this->mainGridLocation[1];
        this->mainGridLocation[0] = oldGridPosition[0];
        this->mainGridLocation[1] = oldGridPosition[1];
        this->removeFromGrid(fullGrid);
        this->mainGridLocation[0] = tempGridPosition[0];
        this->mainGridLocation[1] = tempGridPosition[1];
        this->addToGrid(fullGrid);
    }
}
