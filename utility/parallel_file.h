//
// Created by Chong Peng on 3/24/16.
//

#ifndef TILECLUSTERCHEM_PARALLEL_FILE_H
#define TILECLUSTERCHEM_PARALLEL_FILE_H

namespace mpqc{
namespace utility{

inline void parallel_read_file(madness::World& world, const std::string& filename, char*& buffer){

    int size;
    std::string contents;
    if(world.rank() == 0){
        std::ifstream input_file(filename, std::ifstream::in);

        contents = std::string((std::istreambuf_iterator<char>(input_file)),
                                           std::istreambuf_iterator<char>());

        input_file.close();
        size = contents.size()+1;
    }

    world.gop.broadcast(size,0);
    buffer = new char[size];

    if(world.rank() == 0){
        strcpy(buffer, contents.c_str());
    }

    world.gop.broadcast(buffer,size,0);

}

} // end of namespace utility
} // end of namespace mpqc



#endif //TILECLUSTERCHEM_PARALLEL_FILE_H
