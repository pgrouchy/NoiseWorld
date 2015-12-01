#include <sstream>
#include <sys/stat.h>
#include <sys/file.h>
#include <time.h>
#include <fcntl.h>

#include "mathNet.h"
#include "params.h"

#include <csignal>

void term_trap(int sig)
{
    stringstream final_clean_1;
    final_clean_1 << "killall EMM_island";
    system(final_clean_1.str().c_str());
    cout << "Islands submerged..." << endl;
    stringstream final_clean_2;
    final_clean_2 << "rm -rf /dev/shm/*";
    system(final_clean_2.str().c_str());
    cout << "Term was trapped. RAMdisk cleaned successfully (hopefully?).\n";
    exit(0);
}

int main(int argc, char* argv[]) {

    if (argc != 2) {
        cout << "Usage: EMM_SymComm [runID]" << endl;
        exit(-1);
    }

    signal(SIGTERM, term_trap);

    list<Genome*>::iterator curGenome, curGenome2;
    list<Population*>::iterator curPop;

    string myFilename;

    Population *initPop;

    int i, j, fd, intStat, bestIsl, numInSync, counter1 = 0, counter2 = 0, zipCounter = 0, bestFit = 0, era_counter = 0;
    int myInt[6], myBestInt[6];
    double myDouble[2], myBestDouble[2];
    double densityTotal = 0.0;

    int run_id = atoi(argv[1]);

    FILE *dataFile;
    struct stat stFileInfo;
    ofstream myFile;
    struct flock fl = {F_WRLCK, SEEK_SET,   0,      0,     0 };
    fl.l_pid = getpid();

    string tempString, tempString2;

    // Initialize randomizer so that different runs produce different results
    srand( (unsigned)(run_id*100000) );
    //srand ( time(NULL) );

    // make necessary directories
    for (i = 0; i < NUM_ISLANDS; i++) {
        stringstream isl_sync_dir;
        isl_sync_dir << "mkdir /dev/shm/sync_dir_" << i;
        system(isl_sync_dir.str().c_str());

        stringstream isl_sync_file;
        isl_sync_file << "touch /dev/shm/sync_dir_" << i << "/folderSync.txt";
        system(isl_sync_file.str().c_str());
    }
    stringstream data_dir_ss;
    data_dir_ss << "mkdir dataDir_" << run_id;
    system(data_dir_ss.str().c_str());

    // Start islands
    for (i = 0; i < NUM_ISLANDS; i++) {
        // Generate initial population
        initPop = new Population(i);
        stringstream isl_pop_ss;
        isl_pop_ss << "/dev/shm/initPop_" << i << ".pop";
        initPop->savePop(isl_pop_ss.str().c_str());
        delete initPop;
    }
    for (i = 0; i < NUM_SLAVES; i++) {
        // Run slaves
        stringstream isl_ss;
        isl_ss << "nice ./EMM_island " << i << " " << run_id << " &";
        system(isl_ss.str().c_str());
    }

    // we run the simulation until it is killed by the system or user
    while (true) {
        numInSync = 0;
        for (i = 0; i < NUM_SLAVES; i++) {
            stringstream ss;
            ss << "/dev/shm/fitFlag_" << i << ".txt";
            intStat = stat(ss.str().c_str(),&stFileInfo);
            if (intStat == 0) {
                numInSync++;
            } else {
                break;
            }
        }

        if (numInSync == NUM_SLAVES) {

            era_counter++;

            ofstream myFile;
            stringstream ss_gen_file;
            ss_gen_file << "dataDir_" << run_id << "/genOutput.csv";
            myFile.open(ss_gen_file.str().c_str(), ios::out | ios::app);

            myBestInt[3] = -1;

            for (i = 0; i < NUM_ISLANDS; i++) {

                stringstream sync_ss1;
                sync_ss1 << "/dev/shm/sync_dir_" << i << "/folderSync.txt";
                // set file lock
                fl.l_type = F_WRLCK;
                fd = open(sync_ss1.str().c_str(), O_RDWR);
                fcntl(fd, F_SETLKW, &fl);

                stringstream ss2;
                ss2 << "/dev/shm/currentStatus_" << i << ".txt";
                dataFile = fopen(ss2.str().c_str(), "r");
                fscanf(dataFile, "%d,%d,%d,%d,%d,%lf,%lf,%d\n", &myInt[0], &myInt[1], &myInt[2], &myInt[3], &myInt[4], &myDouble[0], &myDouble[1], &myInt[5]);
                fclose(dataFile);

                myFile << "Isl " << i << " t-step " << myInt[1] << " numBirths: " << myInt[3] << " NumMigIn: " << myInt[4] << " avgAge: " << myDouble[0] << " beta: " << myDouble[1] << " hitWall: " << myBestInt[5] << endl;

                cout << "Isl " << i << " t-step " << myInt[1] << " numBirths: " << myInt[3] << " NumMigIn: " << myInt[4] << " avgAge: " << myDouble[0] << " beta: " << myDouble[1] << " hitWall: " << myBestInt[5] << endl;

                if (myInt[3] > myBestInt[3]) {
                    bestIsl = i;
                    for (j = 0; j < 6; j++) {
                        if (j < 2) {
                            myBestDouble[j] = myDouble[j];
                        }
                        myBestInt[j] = myInt[j];
                    }
                }

                // release lock file
                fl.l_type = F_UNLCK;
                fcntl(fd, F_SETLK, &fl);
                close(fd);
                
                counter1 += myInt[3];
                counter2++;
                densityTotal += myDouble[1];
            }
            myFile.close();

            stringstream copy_best_pop;
            copy_best_pop << "cp /dev/shm/popSnapshot-" << (era_counter * ERA) << "_" << bestIsl << ".pop dataDir_" << run_id << "/curBestPop.pop";
            system(copy_best_pop.str().c_str());

            if (myBestInt[3] > bestFit) {
                bestFit = myBestInt[3];
                stringstream copy_best_pop_2;
                copy_best_pop_2 << "cp /dev/shm/popSnapshot-" << (era_counter * ERA) << "_" << bestIsl << ".pop dataDir_" << run_id << "/overallBestPop.pop";
                system(copy_best_pop_2.str().c_str());
            }

            stringstream ss_avg_file;
            ss_avg_file << "dataDir_" << run_id << "/avgPopSize.csv";
            myFile.open(ss_avg_file.str().c_str(), ios::out | ios::app);
            myFile << ((double)counter1 / (double)counter2) << "  " << (densityTotal / (double)counter2) << endl;
            myFile.close();
            counter1 = 0;
            counter2 = 0;
            densityTotal = 0.0;

            stringstream do_zip;
            do_zip << "/dev/shm/zip_time_flag.txt";
            intStat = stat(do_zip.str().c_str(),&stFileInfo);
            if (intStat == 0) {
                stringstream zip_ss;
                if (zipCounter == 0) {
                    zip_ss << "zip -q -j dataDir_" << run_id << "/popSnapshots_" << zipCounter << ".zip /dev/shm/initPop_* /dev/shm/popSnapshot* /dev/shm/mig*.pop";
                } else {
                    zip_ss << "zip -q -j dataDir_" << run_id << "/popSnapshots_" << zipCounter << ".zip /dev/shm/popSnapshot* /dev/shm/mig*.pop";
                }
                system(zip_ss.str().c_str());

                stringstream zip_rm;
                zip_rm << "rm /dev/shm/popSnapshot*";
                system(zip_rm.str().c_str());

                if (zipCounter == 0) {
                    stringstream del_ss;
                    del_ss << "rm /dev/shm/initPop_*";
                    system(del_ss.str().c_str());
                }

                zipCounter++;
                stringstream rm_zip;
                rm_zip << "rm /dev/shm/zip_time_flag.txt";
                system(rm_zip.str().c_str());
            }


            for (i = 0; i < NUM_SLAVES; i++) {
                stringstream ss;
                ss << "rm /dev/shm/fitFlag_" << i << ".txt";
                system(ss.str().c_str());
            }
        }
        numInSync = 0;
        for (i = 0; i < NUM_SLAVES; i++) {
            stringstream ss;
            ss << "/dev/shm/syncFlag_" << i << ".txt";
            intStat = stat(ss.str().c_str(),&stFileInfo);
            if (intStat == 0) {
                numInSync++;
            } else {
                break;
            }
        }
        if (numInSync == NUM_SLAVES) {
            for (i = 0; i < NUM_SLAVES; i++) {
                stringstream ss;
                ss << "rm /dev/shm/syncFlag_" << i << ".txt";
                system(ss.str().c_str());
            }
        }
        usleep(100000);
    }

    return(0);
}
