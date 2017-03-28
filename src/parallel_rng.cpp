#include <boost/random/uniform_01.hpp>
#include <boost/random.hpp>
#include "parallel_rng.h"

typedef boost::mt19937 rng_t;
static rng_t master_rng;

threadsafe_rng create_uniform_01_threadsafe(){
    threadsafe_rng new_tsrng;

    rng_t* local_rng = new rng_t(master_rng()); //new random generator from seed
    boost::uniform_01<rng_t&>* local_uniform = new boost::uniform_01<rng_t&>(*local_rng);

    new_tsrng.local_rng = (void *) local_rng;
    new_tsrng.local_uniform = (void *) local_uniform;

    return new_tsrng;
}

void delete_uniform_01_threadsafe(threadsafe_rng del_tsrng){
    delete (rng_t *) del_tsrng.local_rng;
    delete (boost::uniform_01<rng_t&> *) del_tsrng.local_uniform;
}

double threadsafe_rng_next(threadsafe_rng tsrng){
    double du01;
    boost::uniform_01<rng_t&>* local_uniform = (boost::uniform_01<rng_t&> *) tsrng.local_uniform;
    du01 = (*local_uniform)();
    return du01;
}


void do_nothing(){
    int a = 1 + 2;
    a = a+1;
}
