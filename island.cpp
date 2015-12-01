#include <sstream>
#include <time.h>
#include <sys/stat.h>
#include <sys/file.h>

#include <fcntl.h>

#include "mathNet.h"
#include "params.h"

int main(int argc, char* argv[]) {

    list<Genome*>::iterator curGenome, curGenome2, curGenome3;
    list<Genome*> moveList, newInds, deathList;
    list<Tree*>::iterator curTree;
    Genome *offspring;
    int i, j, k, fd, wrapX, wrapY, counter, migIsland, intStat, deathCounter;
    int myID = atoi(argv[1]);
    int oldGridPosition[2];
    double my_current_speed, my_current_theta, runningTotal, rouletteWheelSpin, rouletteWheelTotal;
    ofstream myFile;
    string myFilename;
    bool withinReproductionDistance, hitWall;
    bool *deleteList;

    int numIslsPerSlave = NUM_ISLANDS / NUM_SLAVES;
    int myX[numIslsPerSlave];
    int myY[numIslsPerSlave];
    int islIDs[numIslsPerSlave];
    int numDeaths[numIslsPerSlave];
    int numBirths[numIslsPerSlave];
    int migIn[numIslsPerSlave];
    int prevNumBirths[numIslsPerSlave];
    int hitWall_counter[numIslsPerSlave];
    int migInCounter[numIslsPerSlave];
    int migInTotal[numIslsPerSlave];
    double totalAge[numIslsPerSlave];
    double migIncrement[numIslsPerSlave];
    double migIncrementCounter[numIslsPerSlave];

    for (i = 0; i < numIslsPerSlave; i++) {
        islIDs[i] = (myID * numIslsPerSlave) + i;
        hitWall_counter[i] = 0;
        migInCounter[i] = 0;
        migInTotal[i] = 0;

        totalAge[i] = 0.0;
        migIncrement[i] = 0.0;
        migIncrementCounter[i] = 0.0;
    }

    double agent_move[2];

    struct stat stFileInfo;

    string command;

    ofstream resultFile;

    struct flock fl = {F_WRLCK, SEEK_SET,   0,      0,     0 };
    fl.l_pid = getpid();

    // this is used to determine neighboring islands for migration
    int migInds[NUM_ISLANDS];
    counter = 0;
    for (i = 0; i < (int)sqrt(NUM_ISLANDS); i++) {
        for (j = 0; j < (int)sqrt(NUM_ISLANDS); j++) {
            migInds[(j * (int)sqrt(NUM_ISLANDS)) + i] = counter;
            for (k = 0; k < numIslsPerSlave; k++) {
                if (counter == islIDs[k]) {
                    myX[k] = i;
                    myY[k] = j;
                }
            }
            counter++;
        }
    }

    // even though agents exist in a continuous world, a grid is used to speed up nearest neighbor calculations
    list<gridEntry*>::iterator curGridEntry;
    list<gridEntry*> *fullGrid[numIslsPerSlave];
    for (k = 0; k < numIslsPerSlave; k++) {
        fullGrid[k] = new list<gridEntry*>[(int)((double)ISLAND_SIZE / NN_BLOCK_SIZE)*(int)((double)ISLAND_SIZE / NN_BLOCK_SIZE)];
    }

    double *currentData;

    Population *myTempPop;
    Population *myPop[numIslsPerSlave];
    Population *myMigPop[(NUM_ISLANDS + 1) * numIslsPerSlave];

    // Initialize randomizer so that different runs produce different results
    srand( (myID * 10000) + (atoi(argv[2]) * 10000000));

    if (!LOAD_POP) {
        // if we are starting a run from the beginning, load the initial population provided by the master process
        for (k = 0; k < numIslsPerSlave; k++) {
            myPop[k] = new Population();
            stringstream load_ss;
            load_ss << "/dev/shm/initPop_" << islIDs[k] << ".pop";
            command = load_ss.str();
            myPop[k]->loadPop(command.c_str());
            // Initialize individuals
            for (curGenome = myPop[k]->genomes.begin(); curGenome != myPop[k]->genomes.end(); ++curGenome) {
                (*curGenome)->initializeGenome();
                (*curGenome)->generateRandomPosition();
                (*curGenome)->addToGrid(fullGrid[k]);
                (*curGenome)->isNan = false;
            }
            // make sure that no two agents are within the reproduction distance of each other
            for (curGenome = myPop[k]->genomes.begin(); curGenome != myPop[k]->genomes.end(); ++curGenome) {
                oldGridPosition[0] = (*curGenome)->mainGridLocation[0];
                oldGridPosition[1] = (*curGenome)->mainGridLocation[1];
                (*curGenome)->setNearestNeighbor(fullGrid[k], withinReproductionDistance, myPop[k]);
                while (withinReproductionDistance) {
                    (*curGenome)->generateRandomPosition();
                    (*curGenome)->setNearestNeighbor(fullGrid[k], withinReproductionDistance, myPop[k]);

                }
                (*curGenome)->updateOnGrid(fullGrid[k], oldGridPosition);
            }
        }
    } /*else {
        TODO: Load populations from snapshots
    }*/

    for (k = 0; k < numIslsPerSlave; k++) {
        numDeaths[k] = 0;
        numBirths[k] = 0;
        prevNumBirths[k] = 0;
        migIn[k] = 0;

        for (i = 0; i < (NUM_ISLANDS + 1); i++) {
            myMigPop[(k * (NUM_ISLANDS + 1)) + i] = new Population();
        }

        myPop[k]->beta = 1.0;
    }

    // we run the simulation until it is killed by the system or user
    while (true) {

        if (((myPop[0]->gen % SYNC_INTERVAL) == 0) && (myPop[0]->gen > 0)) {

            for (k = 0; k < numIslsPerSlave; k++) {

                // recalculate beta at every sync
                myPop[k]->beta = 500.0 / (numBirths[k] - prevNumBirths[k]);
                if (myPop[k]->beta > 100.0) {
                    myPop[k]->beta = 100.0;
                }

                prevNumBirths[k] = numBirths[k];

                // save population snapshots once per era
                if ((myPop[0]->gen % ERA) == 0) {
                    stringstream ss9;
                    ss9 << "/dev/shm/popSnapshot-" << myPop[0]->gen << "_" << islIDs[k] << ".pop";
                    string command = ss9.str();
                    myPop[k]->savePop(command.c_str());
                }

                // this is used to signal to the master process that it's time to zip and save the population snapshots (here once every 10 eras)
                if ((islIDs[k] == 0) && ((myPop[0]->gen % (10 * ERA)) == 0)) {
                    stringstream do_zip;
                    do_zip << "touch /dev/shm/zip_time_flag.txt";
                    string command = do_zip.str();
                    system(command.c_str());
                }

                if ((myPop[0]->gen % ERA) == 0) {

                    stringstream sync_ss1;
                    sync_ss1 << "/dev/shm/sync_dir_" << islIDs[k] << "/folderSync.txt";
                    command = sync_ss1.str();
                    // set file lock
                    fl.l_type = F_WRLCK;
                    fd = open(command.c_str(), O_RDWR);
                    // waits until file lock is free (if it is already set), then sets file lock
                    fcntl(fd, F_SETLKW, &fl);

                    // send updates to master process
                    stringstream ss4;
                    ss4 << "/dev/shm/currentStatus_" << islIDs[k] << ".txt";
                    myFilename = ss4.str();
                    resultFile.open(myFilename.c_str(), ios::out | ios::trunc);
                    resultFile << myPop[k]->pop_size << "," << myPop[0]->gen << "," << numDeaths[k] << "," << numBirths[k] << "," << migIn[k] << "," << totalAge[k] / (double)(numDeaths[k]) << "," << myPop[k]->beta << "," << hitWall_counter[k] << endl;
                    resultFile.close();

                    totalAge[k] = 0.0;

                    numDeaths[k] = 0;
                    numBirths[k] = 0;
                    prevNumBirths[k] = 0;
                    migIn[k] = 0;
                    hitWall_counter[k] = 0;

                    // this flag lets the master process know that this island has updated its statistics
                    if (k == (numIslsPerSlave - 1)) {
                        stringstream ss5;
                        ss5 << "touch /dev/shm/fitFlag_" << myID << ".txt";
                        myFilename = ss5.str();
                        system(myFilename.c_str());
                    }

                    // release lock file
                    fl.l_type = F_UNLCK;
                    fcntl(fd, F_SETLK, &fl);
                    close(fd);
                }
            }

            // create a sync flag file
            // while this flag is present, the island will wait (the master process will remove this flag once all islands have reached this point)
            stringstream isl_sync_1;
            isl_sync_1 << "touch /dev/shm/syncFlag_" << myID << ".txt";
            command = isl_sync_1.str();
            system(command.c_str());

            // sync (wait until the master process removes this island's sync flag)
            stringstream isl_sync_2;
            isl_sync_2 << "/dev/shm/syncFlag_" << myID << ".txt";
            command = isl_sync_2.str();
            do {
                usleep(500000);
                intStat = stat(command.c_str(),&stFileInfo);
            } while (intStat == 0);
        }

        for (k = 0; k < numIslsPerSlave; k++) {

            // prepare agents for the next timestep
            for (curGenome = myPop[k]->genomes.begin(); curGenome != myPop[k]->genomes.end(); ++curGenome) {
                (*curGenome)->currentComm = (*curGenome)->myData[1] + (NOISE_STD * gaussrand());
                (*curGenome)->numOffspringProduced = 0;
            }

            // prepare to delete agents (1 death per birth in the previous timestep)
            deleteList = new bool[myPop[k]->pop_size];
            for (i = 0; i < myPop[k]->pop_size; i++) {
                deleteList[i] = false;
            }

            deathList.clear();
            deathCounter = myPop[k]->pop_size - ISLAND_POP; // check if deaths are needed
            if (deathCounter > 0) {
                for (curGenome = myPop[k]->genomes.begin(); curGenome != myPop[k]->genomes.end(); ++curGenome) {
                    // first select agents that have hit a wall or exceeded numerical limits
                    if ((*curGenome)->isNan) {
                        (*curGenome)->is_alive = false;
                        numDeaths[k]++;
                        totalAge[k] += (*curGenome)->age;
                        deathList.push_back((*curGenome));
                        deathCounter--;
                        if (deathCounter == 0) {
                            break;
                        }
                    }
                }
            }

            // if deaths are still required, select agents randomly
            if (deathCounter > 0) {
                rouletteWheelTotal = 0.0;
                for (curGenome = myPop[k]->genomes.begin(); curGenome != myPop[k]->genomes.end(); ++curGenome) {
                    // newborn agents are not eligible to die
                    if ((*curGenome)->age > 0) {
                        rouletteWheelTotal += 1.0;
                    }
                }

                while (deathCounter > 0) {
                    rouletteWheelSpin = rouletteWheelTotal * randfloat();
                    runningTotal = 0.0;
                    for (curGenome = myPop[k]->genomes.begin(); curGenome != myPop[k]->genomes.end(); ++curGenome) {
                        if ((*curGenome)->age > 0) {
                            runningTotal += 1.0;
                            if (rouletteWheelSpin <= runningTotal)
                                break;
                        }
                    }
                    if ((*curGenome)->is_alive) {
                        (*curGenome)->is_alive = false;
                        numDeaths[k]++;
                        totalAge[k] += (*curGenome)->age;
                        deathList.push_back((*curGenome));
                        deathCounter--;
                    }
                }
            }

            counter = 0;
            // evaluate each genome
            for (curGenome = myPop[k]->genomes.begin(); curGenome != myPop[k]->genomes.end(); ++curGenome) {

                if ((*curGenome)->is_alive) {

                    (*curGenome)->age++;

                    currentData = new double[NUM_INPUTS + (*curGenome)->num_trees];

                    // set current agent inputs             
                    currentData[0] = (*curGenome)->myPosition[0] + (NOISE_STD * gaussrand());
                    currentData[1] = (*curGenome)->myPosition[1] + (NOISE_STD * gaussrand());
                    currentData[2] = (*curGenome)->closestInd->currentComm + (NOISE_STD * gaussrand());

                    // set inputs from previous agent outputs and internal variables
                    for (i = 0; i < (*curGenome)->num_trees; i++) {
                        currentData[i + NUM_INPUTS] = (*curGenome)->myData[i] + (NOISE_STD * gaussrand());
                    }

                    // evaluate the agent's EMM
                    (*curGenome)->evalEMM(currentData);

                    delete [] currentData;

                    // check for NaN or out of bounds outputs, which will cause the agent to die at the next island reproduction event
                    for (i = 0; i < (*curGenome)->num_trees; i++) {
                        if ((isnan((*curGenome)->myData[i])) || ((*curGenome)->myData[i] == DBL_MAX) || ((*curGenome)->myData[i] == -DBL_MAX)) {
                            (*curGenome)->isNan = true;
                            (*curGenome)->myData[i] = (DBL_MAX * randfloat()) - (DBL_MAX / 2.0);
                        }
                    }

                    // this is where you decide how theta is counted relative to the x/y coordinate system
                    my_current_theta = (*curGenome)->myData[0] + (NOISE_STD * gaussrand());// currently theta = 0 points east, add "+ (M_PI / 2.0)" to make it point north

                    if ((my_current_theta >= DBL_MAX) || (my_current_theta <= -DBL_MAX)) {
                        my_current_theta = (2.0 * M_PI * randfloat()) - M_PI;
                        (*curGenome)->isNan = true;
                    }

                    // calculate agent's next move
                    my_current_speed = 1.0 + (NOISE_STD * gaussrand());
                    agent_move[0] = my_current_speed * cos(my_current_theta);
                    agent_move[1] = my_current_speed * sin(my_current_theta);


                    // update agent position
                    (*curGenome)->myPosition[0] += STEP_SIZE * agent_move[0];
                    (*curGenome)->myPosition[1] += STEP_SIZE * agent_move[1];

                    // check to see if the agent has hit a wall
                    // if so, agent bounces off and is set to die at next island reproduction event
                    hitWall = false;
                    if ((*curGenome)->myPosition[0] >= ((double)ISLAND_SIZE / 2.0)) {
                        (*curGenome)->myPosition[0] = ((double)ISLAND_SIZE - (*curGenome)->myPosition[0]) - 0.000001;
                        hitWall = true;
                    } else if ((*curGenome)->myPosition[0] <= ((double)ISLAND_SIZE / -2.0)) {
                        (*curGenome)->myPosition[0] = -1.0 * ((double)ISLAND_SIZE + (*curGenome)->myPosition[0] - 0.000001);
                        hitWall = true;
                    }
                    if ((*curGenome)->myPosition[1] >= ((double)ISLAND_SIZE / 2.0)) {
                        (*curGenome)->myPosition[1] = ((double)ISLAND_SIZE - (*curGenome)->myPosition[1]) - 0.000001;
                        hitWall = true;
                    } else if ((*curGenome)->myPosition[1] <= ((double)ISLAND_SIZE / -2.0)) {
                        (*curGenome)->myPosition[1] = -1.0 * ((double)ISLAND_SIZE + (*curGenome)->myPosition[1] - 0.000001);
                        hitWall = true;
                    }
                    if (hitWall) {
                        if (!(*curGenome)->isNan) {
                            hitWall_counter[k]++;
                        }
                        (*curGenome)->isNan = true;
                    }


                    // update position in grid
                    oldGridPosition[0] = (*curGenome)->mainGridLocation[0];
                    oldGridPosition[1] = (*curGenome)->mainGridLocation[1];
                    (*curGenome)->mainGridLocation[0] = (int)(((*curGenome)->myPosition[0] + ((double)ISLAND_SIZE / 2.0)) / NN_BLOCK_SIZE);
                    (*curGenome)->mainGridLocation[1] = (int)(((*curGenome)->myPosition[1] + ((double)ISLAND_SIZE / 2.0)) / NN_BLOCK_SIZE);
                    (*curGenome)->updateOnGrid(fullGrid[k], oldGridPosition);

                } else {
                    deleteList[counter] = true;
                }

                counter++;
            }

            // remove dead individuals from population
            i = 0;
            curGenome = myPop[k]->genomes.begin();
            while (curGenome != myPop[k]->genomes.end()) {
                if (deleteList[i]) {
                    (*curGenome)->removeFromGrid(fullGrid[k]);
                    delete (*curGenome);
                    curGenome = myPop[k]->genomes.erase(curGenome);
                    myPop[k]->pop_size--;
                } else {
                    ++curGenome;
                }
                i++;
            }
            delete [] deleteList;

            // calculate nearest neighbours, check for reproduction, do reproduction, re-calculate nearest neighbours (if reproduction occurs at least once)

            // check for incoming migrants (is migFlag set?), give random location, make sure they do not end up within reproductive distance
            // this allows for multiple migrants in the same timestep
            if ((myPop[0]->gen % MIGRATION_IN_INTERVAL) == 1) {

                stringstream sync_ss1;
                sync_ss1 << "/dev/shm/sync_dir_" << islIDs[k] << "/folderSync.txt";
                command = sync_ss1.str();
                // set file lock
                fl.l_type = F_WRLCK;
                fd = open(command.c_str(), O_RDWR);
                // waits until file lock is free (if it is already set), then sets file lock
                fcntl(fd, F_SETLKW, &fl);

                stringstream ss;
                ss << "/dev/shm/migFlag_" << islIDs[k] << ".txt";
                command = ss.str();
                intStat = stat(command.c_str(),&stFileInfo);;

                if (intStat == 0) {

                    stringstream ss2;
                    ss2 << "/dev/shm/migPop_" << islIDs[k] << ".pop";
                    command = ss2.str();
                    myTempPop = new Population();
                    myTempPop->loadPop(command.c_str());

                    stringstream ss;
                    ss << "/dev/shm/migFlag_" << islIDs[k] << ".txt";
                    command = ss.str();
                    remove(command.c_str());

                    // randomize the order of incoming migrants
                    while ((int)myTempPop->genomes.size() > 0) {
                        i = (int)(randfloat() * (int)myTempPop->genomes.size());
                        j = 0;
                        for (curGenome2 = myTempPop->genomes.begin(); curGenome2 != myTempPop->genomes.end(); ++curGenome2) {
                            if (i == j) {
                                break;
                            }
                            j++;
                        }
                        (*curGenome2)->isNan = false;
                        myMigPop[(k * (NUM_ISLANDS + 1)) + NUM_ISLANDS]->genomes.push_back(new Genome((*curGenome2)));
                        myMigPop[(k * (NUM_ISLANDS + 1)) + NUM_ISLANDS]->pop_size++;
                        delete (*curGenome2);
                        myTempPop->genomes.erase(curGenome2);
                    }
                    delete myTempPop;
                }
                // release lock file
                fl.l_type = F_UNLCK;
                fcntl(fd, F_SETLK, &fl);
                close(fd);

                migIncrement[k] = (double)myMigPop[(k * (NUM_ISLANDS + 1)) + NUM_ISLANDS]->pop_size / ((double)MIGRATION_IN_INTERVAL - 2.0);
                migIncrementCounter[k] = 0.0;
                migInTotal[k] = myMigPop[(k * (NUM_ISLANDS + 1)) + NUM_ISLANDS]->pop_size;
                migInCounter[k] = 0;
            }

            // migration code to take in all migrants, then release at regular intervals
            migIncrementCounter[k] += migIncrement[k];
            while ((migInCounter[k] < migInTotal[k]) && (migInCounter[k] < migIncrementCounter[k])) {
                curGenome2 = myMigPop[(k * (NUM_ISLANDS + 1)) + NUM_ISLANDS]->genomes.begin();
                (*curGenome2)->is_alive = true;
                newInds.push_back(new Genome(*curGenome2));
                migIn[k]++;
                migInCounter[k]++;
                delete (*curGenome2);
                myMigPop[(k * (NUM_ISLANDS + 1)) + NUM_ISLANDS]->genomes.pop_front();
                myMigPop[(k * (NUM_ISLANDS + 1)) + NUM_ISLANDS]->pop_size--;
            }

            // reproduction
            for (curGenome = myPop[k]->genomes.begin(); curGenome != myPop[k]->genomes.end(); ++curGenome) {
                // recalculate nearest neighbors (agents have moved)
                (*curGenome)->setNearestNeighbor(fullGrid[k], withinReproductionDistance, myPop[k]);
                // if the 1st potential parent hasn't yet reproduced and the two potential parents are within the reproduction distance, they produce one offspring
                if (((*curGenome)->numOffspringProduced == 0) && (withinReproductionDistance)) {
                    // call reproduction
                    offspring = (*curGenome)->reproduce((*curGenome)->closestInd, myPop[k]->my_isl_id, myPop[k]->max_tree_id, myPop[k]->beta);

                    offspring->initializeGenome();

                    // assign a temporary position
                    offspring->myPosition[0] = (*curGenome)->myPosition[0];
                    offspring->myPosition[1] = (*curGenome)->myPosition[1];

                    (*curGenome)->numOffspringProduced++;
                    (*curGenome)->closestInd->numOffspringProduced++;

                    (*curGenome)->fitness++;
                    (*curGenome)->closestInd->fitness++;

                    // parents will be moved to new random locations
                    moveList.push_back((*curGenome));
                    moveList.push_back((*curGenome)->closestInd);

                    // perform equation reduction
                    if (randfloat() < REDUCE_PROB) {
                        for (curTree = offspring->trees.begin(); curTree != offspring->trees.end(); ++curTree) {
                            (*curTree)->reduceEqn();
                            (*curTree)->recountNodes(0);
                        }
                    }

                    // check for migration, if so save offspring to disk and delete
                    // this allows for multiple migrants in the same timestep
                    if (!offspring->is_alive) { // if an offspring is born with more than MAX_TOTAL_NODES, it dies immediately
                        delete offspring;
                    } else if (randfloat() < (MIGRATION_PROB * myPop[k]->beta)) { // if the offspring is going to migrate, choose a random neighboring island to receive it
                        do {
                            i = (int)(3.0 * randfloat()) - 1;
                            j = (int)(3.0 * randfloat()) - 1;
                        } while (((i == 0) && (j == 0)) || ((abs(i) == 1) && (abs(j) == 1)));
                        wrapX = myX[k] + i;
                        if (wrapX < 0) {
                            wrapX += (int)sqrt(NUM_ISLANDS);
                        }
                        if (wrapX >= (int)sqrt(NUM_ISLANDS)) {
                            wrapX -= (int)sqrt(NUM_ISLANDS);
                        }
                        wrapY = myY[k] + j;
                        if (wrapY < 0) {
                            wrapY += (int)sqrt(NUM_ISLANDS);
                        }
                        if (wrapY >= (int)sqrt(NUM_ISLANDS)) {
                            wrapY -= (int)sqrt(NUM_ISLANDS);
                        }
                        migIsland = migInds[(wrapY * (int)sqrt(NUM_ISLANDS)) + wrapX];

                        myMigPop[(k * (NUM_ISLANDS + 1)) + migIsland]->genomes.push_back(offspring);
                        myMigPop[(k * (NUM_ISLANDS + 1)) + migIsland]->pop_size++;
                    } else { // offspring born on current island
                        newInds.push_back(offspring);
                        numBirths[k]++;
                    }
                }
            }

            myPop[k]->pop_size += (int)newInds.size();

            // save outgoing migrants to disk
            if ((myPop[0]->gen % MIGRATION_IN_INTERVAL) == (MIGRATION_IN_INTERVAL - 1)) {

                for (migIsland = 0; migIsland < NUM_ISLANDS; migIsland++) {
                    if ((migIsland != islIDs[k]) && (myMigPop[(k * (NUM_ISLANDS + 1)) + migIsland]->pop_size > 0)) {
                        stringstream ss2;
                        ss2 << "/dev/shm/migPop_" << migIsland << ".pop";
                        stringstream ss3;
                        ss3 << "/dev/shm/migFlag_" << migIsland << ".txt";

                        stringstream sync_ss1;
                        sync_ss1 << "/dev/shm/sync_dir_" << migIsland << "/folderSync.txt";
                        command = sync_ss1.str();
                        // set file lock
                        fl.l_type = F_WRLCK;
                        fd = open(command.c_str(), O_RDWR);
                        // waits until file lock is free (if it is already set), then sets file lock
                        fcntl(fd, F_SETLKW, &fl);
                        
                        myFilename = ss3.str();
                        intStat = stat(myFilename.c_str(),&stFileInfo);

                        if (intStat != 0) {
                            // set the flag to indicate that new migrants are ready
                            stringstream sync_ss3;
                            sync_ss3 << "touch /dev/shm/migFlag_" << migIsland << ".txt";
                            command = sync_ss3.str();
                            system(command.c_str());
                        } else {
                            // if the flag is already set, load the current new migrant population, add the migrants leaving from this island
                            myFilename = ss2.str();
                            myTempPop = new Population();
                            myTempPop->loadPop(myFilename.c_str());
                            for (curGenome2 = myTempPop->genomes.begin(); curGenome2 != myTempPop->genomes.end(); ++curGenome2) {
                                myMigPop[(k * (NUM_ISLANDS + 1)) + migIsland]->genomes.push_back(new Genome((*curGenome2)));
                                myMigPop[(k * (NUM_ISLANDS + 1)) + migIsland]->pop_size++;
                            }
                            delete myTempPop;
                        }

                        // save migrant population file
                        myFilename = ss2.str();
                        myMigPop[(k * (NUM_ISLANDS + 1)) + migIsland]->savePop(myFilename.c_str());

                        // release lock file
                        fl.l_type = F_UNLCK;
                        fcntl(fd, F_SETLK, &fl);
                        close(fd);

                        delete myMigPop[(k * (NUM_ISLANDS + 1)) + migIsland];
                        myMigPop[(k * (NUM_ISLANDS + 1)) + migIsland] = new Population();
                    }
                }
            }

            // move parents and new offspring/incoming migrants to random positions
            // make sure that these agents are at least one reproduction distance from their nearest neighbor
            if (((int)moveList.size() > 0) || ((int)newInds.size() > 0)) {
                for (curGenome = moveList.begin(); curGenome != moveList.end(); ++curGenome) {
                    (*curGenome)->initializeGenome();
                    oldGridPosition[0] = (*curGenome)->mainGridLocation[0];
                    oldGridPosition[1] = (*curGenome)->mainGridLocation[1];
                    do {
                        (*curGenome)->mainGridLocation[0] = oldGridPosition[0];
                        (*curGenome)->mainGridLocation[1] = oldGridPosition[1];
                        (*curGenome)->generateRandomPosition();
                        (*curGenome)->setNearestNeighbor(fullGrid[k], withinReproductionDistance, myPop[k]);

                    } while (withinReproductionDistance);
                    (*curGenome)->updateOnGrid(fullGrid[k], oldGridPosition);
                }
                for (curGenome = newInds.begin(); curGenome != newInds.end(); ++curGenome) {
                    (*curGenome)->initializeGenome();
                    oldGridPosition[0] = (*curGenome)->mainGridLocation[0];
                    oldGridPosition[1] = (*curGenome)->mainGridLocation[1];
                    do {
                        (*curGenome)->mainGridLocation[0] = oldGridPosition[0];
                        (*curGenome)->mainGridLocation[1] = oldGridPosition[1];
                        (*curGenome)->generateRandomPosition();
                        (*curGenome)->setNearestNeighbor(fullGrid[k], withinReproductionDistance, myPop[k]);
                    } while (withinReproductionDistance);
                    (*curGenome)->addToGrid(fullGrid[k]);
                    myPop[k]->genomes.push_back((*curGenome));
                }
                // recalculate all nearest neighbors
                for (curGenome = myPop[k]->genomes.begin(); curGenome != myPop[k]->genomes.end(); ++curGenome) {
                    (*curGenome)->setNearestNeighbor(fullGrid[k], myPop[k]);
                }
                moveList.clear();
                newInds.clear();
            }
        }
        myPop[0]->gen++;
    }

    for (k = 0; k < numIslsPerSlave; k++) {
        for (i = 0; i < ((int)(ISLAND_SIZE / NN_BLOCK_SIZE)*(int)(ISLAND_SIZE / NN_BLOCK_SIZE)); i++) {
            for (curGridEntry = fullGrid[k][i].begin(); curGridEntry != fullGrid[k][i].end(); ++curGridEntry) {
                delete (*curGridEntry);
            }
            fullGrid[k][i].clear();
        }
        delete [] fullGrid[k];
    }

    for (k = 0; k < numIslsPerSlave; k++) {
        delete myPop[k];
        
    }

    for (k = 0; k < numIslsPerSlave; k++) {
        for (i = 0; i < (NUM_ISLANDS + 1); i++) {
            delete myMigPop[(k * (NUM_ISLANDS + 1)) + i];
        }
    }
}
