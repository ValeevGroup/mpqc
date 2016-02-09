#include "../include/libint.h"
#include "../include/tiledarray.h"

int main(int argc, char** argv){
    auto &world = madness::initialize(argc, argv);
    world.gop.fence();
    return 0;
}
