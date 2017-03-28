#include <boost/random/uniform_01.hpp>
#include <boost/random.hpp>
#include "parallel_rng.h"

typedef boost::mt19937 rng_t;
typedef boost::uniform_01<rng_t&> rng_unif_t;
static rng_t master_rng;

threadsafe_rng create_threadsafe_rng(){
    threadsafe_rng new_tsrng;

    rng_t* local_rng = new rng_t(master_rng()); //new random generator from seed
    rng_unif_t* local_uniform = new rng_unif_t(*local_rng);

    new_tsrng.local_rng = (void *) local_rng;
    new_tsrng.local_uniform = (void *) local_uniform;

    return new_tsrng;
}

void delete_threadsafe_rng(threadsafe_rng del_tsrng){
    rng_t* to_del_rng = static_cast<rng_t*>(del_tsrng.local_rng);
    rng_unif_t* to_del_unif = static_cast<rng_unif_t*>(del_tsrng.local_uniform);
    delete to_del_unif;
    delete to_del_rng;
}

double threadsafe_rng_next(threadsafe_rng tsrng){
    double du01;
    rng_unif_t* local_uniform = static_cast<rng_unif_t*>(tsrng.local_uniform);
    du01 = (*local_uniform)();
    return du01;
}


void do_nothing(){
    int a = 1 + 2;
    a = a+1;
}
