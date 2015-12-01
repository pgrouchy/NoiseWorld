#define LOAD_POP false

#define ISLAND_POP 50 // fixed island population size
#define ISLAND_SIZE 40 // each island will be ISLAND_SIZE x ISLAND_SIZE

#define NUM_ISLANDS 100 // ***IMPORTANT: sqrt(NUM_ISLANDS) must yield a whole number***
#define NUM_SLAVES 20 //   ***IMPORTANT: (NUM_ISLANDS / NUM_SLAVES) must yield a whole number***
#define NUM_INPUTS 3 // x, y, current comm in
#define NUM_OUTPUTS 2 // theta, current comm out
#define STEP_SIZE 0.0005
#define NOISE_STD 0.025 // percentage of total range input and output gaussian noise
#define REPRODUCTION_DIST 0.138888889
#define START_RADIUS 1.128379167
#define ERA 100000 // ***IMPORTANT: (ERA / SYNC_INTERVAL) must yield a whole number***
#define SYNC_INTERVAL 10000

#define NN_BLOCK_SIZE 0.25 // the world will get partitioned into blocks of NN_BLOCK_SIZE x NN_BLOCK_SIZE for nearest neighbour calculation purposes
                           // ***IMPORTANT: (ISLAND_SIZE / NN_BLOCK_SIZE) must yield a whole number***

#define MIGRATION_PROB 0.001 // probability that an offspring will be moved to a new random island
#define MIGRATION_IN_INTERVAL 10000 // number of timesteps between receiving migrants
#define REDUCE_PROB 0.1 // probability that an offsprings equations will undergo reduction

#define PERCENT_CONSTANT 0.5
#define PERCENT_TERMINAL 0.5
#define MAX_CONSTANT 5.0

#define MUT_RATE 0.025
#define TREE_SWAP_RATE 0.5
#define PROB_NON_TERM_XOVER 0.9 //terminal xover prob is 1 - PROB_NON_TERM_XOVER
#define PROB_POINT_MUT 0.5 //probability that a tree mutation will be a point mutation
#define PROB_MUT_CONST 0.5 //probability that a point mutation will mutate a constant
#define POINT_MUT_STD 0.5 //mutation standard dev
#define PROB_MUT_TERM_TYPE 0.5 //probability of switching from a constant to a variable or vise versa
#define MUT_INITIAL_COND_STD 0.25

#define MAX_TOTAL_NODES 200 // if offspring has more than TOTAL nodes during reproduction, it is discarded

#define MAX_BUFFER 1000
