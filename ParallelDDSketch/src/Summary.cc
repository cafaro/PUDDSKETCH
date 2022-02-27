/**
 * @file Summary.cc
 * @author CM
 */

#include "Summary.h"



//****************************************************** SERIALIZATION

//*** Used to pack a positive and a negative sketch
char *pack_SketchesPair(double alpha, int collapses, int bound, std::map<int,long>& posiSketch, int posibins, long posipop, std::map<int,long>& negaSketch, int negabins, long negapop, int *final_size) {

    int buffer_size = sizeof(double) + sizeof(int) + sizeof(int);   // alpha + collapses + sketch bound
    buffer_size += sizeof(int) + sizeof(long);                      // posiBins + posiPopulation
    buffer_size += sizeof(int) + sizeof(long);                      // negaBins + negaPopulation

    //buffer_size += bound * (2 * sizeof(int));                     // max number of (key,count) tuples (PADDING the PACK)
    buffer_size += bound * (sizeof(int) + sizeof(long));            // max number of (key,count) tuples
        
    (*final_size) = buffer_size;

    char *localBuffer = (char *) malloc( sizeof(char) * buffer_size);
    if (!localBuffer) {
        return NULL;
    }
    memset(localBuffer, 0, sizeof(char) * buffer_size); 

    std::map<int,long>::iterator it;
    int buffer_position = 0;
    
    MPI_Pack(&alpha, 1, MPI_DOUBLE, localBuffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(&collapses, 1, MPI_INT, localBuffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(&bound, 1, MPI_INT, localBuffer, buffer_size, &buffer_position, MPI_COMM_WORLD);

    MPI_Pack(&posibins, 1, MPI_INT, localBuffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(&posipop, 1, MPI_LONG, localBuffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
    //MPI_Pack(&posipop, 1, MPI_INT, localBuffer, buffer_size, &buffer_position, MPI_COMM_WORLD);

    #ifndef NDEBUG
        long npop = 0, ppop = 0;
    #endif

    for(it=posiSketch.begin(); it != posiSketch.end(); ++it) {

        MPI_Pack(&(it->first), 1, MPI_INT, localBuffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
        MPI_Pack(&(it->second), 1, MPI_LONG, localBuffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
        //MPI_Pack(&(it->second), 1, MPI_INT, localBuffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
        
        #ifndef NDEBUG
            ppop += it->second;
        #endif
    }//posi
    
    MPI_Pack(&negabins, 1, MPI_INT, localBuffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(&negapop, 1, MPI_LONG, localBuffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
    //MPI_Pack(&negapop, 1, MPI_INT, localBuffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
    
    
    for(it=negaSketch.begin(); it != negaSketch.end(); ++it){
        MPI_Pack(&(it->first), 1, MPI_INT, localBuffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
        MPI_Pack(&(it->second), 1, MPI_LONG, localBuffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
        //MPI_Pack(&(it->second), 1, MPI_INT, localBuffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
        
        #ifndef NDEBUG
            npop += it->second;
        #endif
    }//nega

    //std::cout<<"ppop: "<<ppop<<", npop: "<<npop<<std::endl;
    #ifndef NDEBUG
        assert((npop + ppop) == (posipop+negapop));
    #endif

    return localBuffer;
}


//****************************************************** DESERIALIZATION

//*** Used to deserialize into a positive and a negative sketch representation
struct PackedPair *deserializeMsgPair(void *buffer, int size) {

    int pos = 0;
    
    struct PackedPair *res = (struct PackedPair *) malloc( sizeof(struct PackedPair) );
    if (!res) {
        return res;
    }
    memset(res, 0, sizeof(struct PackedPair));


    MPI_Unpack(buffer, size, &pos, &(res->alpha), 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(buffer, size, &pos, &(res->collapses), 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(buffer, size, &pos, &(res->bound), 1, MPI_INT, MPI_COMM_WORLD);
    
    MPI_Unpack(buffer, size, &pos, &(res->posibins), 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(buffer, size, &pos, &(res->posipop), 1, MPI_LONG, MPI_COMM_WORLD);
    //MPI_Unpack(buffer, size, &pos, &(res->posipop), 1, MPI_INT, MPI_COMM_WORLD);

    #ifndef NDEBUG
        long npop = 0, ppop = 0;
    #endif

    if (res->posibins > 0) {

        res->posi = (struct Bucket *)malloc(sizeof(struct Bucket)*(res->posibins));    

        //MPI_Unpack(buffer, size, &pos, res->posi, 2*(res->posibins), MPI_INT, MPI_COMM_WORLD);

        for(int b=0; b<res->posibins; ++b) {

            MPI_Unpack(buffer, size, &pos, &(res->posi[b].key), 1, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, size, &pos, &(res->posi[b].count), 1, MPI_LONG, MPI_COMM_WORLD);

            #ifndef NDEBUG
                ppop += res->posi[b].count;
            #endif

        }//for b
    }//fi posibins
    
    #ifndef NDEBUG
        assert(ppop == res->posipop);
    #endif
    
    MPI_Unpack(buffer, size, &pos, &(res->negabins), 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(buffer, size, &pos, &(res->negapop), 1, MPI_LONG, MPI_COMM_WORLD);
    //MPI_Unpack(buffer, size, &pos, &(res->negapop), 1, MPI_INT, MPI_COMM_WORLD);

    if (res->negabins > 0) {

        res->nega = (struct Bucket *)malloc(sizeof(struct Bucket)*res->negabins);
        
        //MPI_Unpack(buffer, size, &pos, res->nega, 2*(res->negabins), MPI_INT, MPI_COMM_WORLD);
        for(int b=0; b<res->negabins; ++b) {

            MPI_Unpack(buffer, size, &pos, &(res->nega[b].key), 1, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, size, &pos, &(res->nega[b].count), 1, MPI_LONG, MPI_COMM_WORLD);

            #ifndef NDEBUG
                npop += res->nega[b].count;
            #endif
        }//for b
    }//fi negabins

    #ifndef NDEBUG
        assert(npop == res->negapop);
    #endif

    return res;
}

