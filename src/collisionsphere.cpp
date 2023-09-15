#include "collisionsphere.h"


/*******************************************   Cylinder Sphere Collision Implementation ***********************************************/


CylinderCollisionSphere::CylinderCollisionSphere():small_sphere_list_end(0),big_sphere_list_end(0)
{
}

void CylinderCollisionSphere::pushToSmallSphere(unsigned i)
{
    //If i is already inside the "good" side we do nothing
    if(i < small_sphere_list_end || small_sphere_list_end == collision_list->size()){
        return;
    }

    unsigned jkr = collision_list->at(i);
    collision_list->at(i) = collision_list->at(small_sphere_list_end);
    collision_list->at(small_sphere_list_end) = jkr;
    small_sphere_list_end++;

    // WARNING small sphere size should never be greater than the big one.
    if(small_sphere_list_end > big_sphere_list_end)
        big_sphere_list_end = small_sphere_list_end;
}

void CylinderCollisionSphere::pushToBigSphere(unsigned i)
{
    //If i is already on the "other" side we do nothing
    if(i < big_sphere_list_end || big_sphere_list_end == collision_list->size()){
        return;
    }

    unsigned jkr = collision_list->at(i);
    collision_list->at(i)=collision_list->at(big_sphere_list_end);
    collision_list->at(big_sphere_list_end) = jkr;
    big_sphere_list_end++;

}


void CylinderCollisionSphere::popFromSmallSphere(unsigned i)
{
    //If i is already on the "other" side we do nothing
    if(i >= small_sphere_list_end || small_sphere_list_end == 0){
        return;
    }

    unsigned jkr = collision_list->at(i);
    collision_list->at(i)=collision_list->at(small_sphere_list_end-1);
    collision_list->at(small_sphere_list_end-1) = jkr;
    small_sphere_list_end--;
}

void CylinderCollisionSphere::popFromBigSphere(unsigned i)
{
    //If i is already on the "other" side we do nothing
    if(i >= big_sphere_list_end || big_sphere_list_end == 0){
        return;
    }

    unsigned jkr = collision_list->at(i);
    collision_list->at(i)=collision_list->at(big_sphere_list_end-1);
    collision_list->at(big_sphere_list_end-1) = jkr;
    big_sphere_list_end--;

    // WARNING small sphere size should never be greater than the big one.
    if(big_sphere_list_end < small_sphere_list_end)
        small_sphere_list_end = big_sphere_list_end;
}


void CylinderCollisionSphere::setBigSphereSize(float size){
    big_sphere_distance = size;
}

void CylinderCollisionSphere::setSmallSphereSize(float size){
    small_sphere_distance = size;
}

void CylinderCollisionSphere::push_index(unsigned int element)
{
        collision_list->push_back(element);
        list_size++;
}

/************************  Axon Collision Implementation  ******************/

AxonCollisionSphere::AxonCollisionSphere():small_sphere_list_end(0),big_sphere_list_end(0)
{
}

void AxonCollisionSphere::pushToSmallSphere(unsigned i)
{
    //If i is already inside the "good" side we do nothing
    if(i < small_sphere_list_end || small_sphere_list_end == collision_list->size()){
        return;
    }

    unsigned jkr = collision_list->at(i);
    collision_list->at(i) = collision_list->at(small_sphere_list_end);
    collision_list->at(small_sphere_list_end) = jkr;
    small_sphere_list_end++;

    // WARNING small sphere size should never be greater than the big one.
    if(small_sphere_list_end > big_sphere_list_end)
        big_sphere_list_end = small_sphere_list_end;
}

void AxonCollisionSphere::pushToBigSphere(unsigned i)
{
    //If i is already on the "other" side we do nothing
    if(i < big_sphere_list_end || big_sphere_list_end == collision_list->size()){
        return;
    }

    unsigned jkr = collision_list->at(i);
    collision_list->at(i)=collision_list->at(big_sphere_list_end);
    collision_list->at(big_sphere_list_end) = jkr;
    big_sphere_list_end++;

}


void AxonCollisionSphere::popFromSmallSphere(unsigned i)
{
    //If i is already on the "other" side we do nothing
    if(i >= small_sphere_list_end || small_sphere_list_end == 0){
        return;
    }

    unsigned jkr = collision_list->at(i);
    collision_list->at(i)=collision_list->at(small_sphere_list_end-1);
    collision_list->at(small_sphere_list_end-1) = jkr;
    small_sphere_list_end--;
}

void AxonCollisionSphere::popFromBigSphere(unsigned i)
{
    //If i is already on the "other" side we do nothing
    if(i >= big_sphere_list_end || big_sphere_list_end == 0){
        return;
    }

    unsigned jkr = collision_list->at(i);
    collision_list->at(i)=collision_list->at(big_sphere_list_end-1);
    collision_list->at(big_sphere_list_end-1) = jkr;
    big_sphere_list_end--;

    // WARNING small sphere size should never be greater than the big one.
    if(big_sphere_list_end < small_sphere_list_end)
        small_sphere_list_end = big_sphere_list_end;
}


void AxonCollisionSphere::setBigSphereSize(float size){
    big_sphere_distance = size;
}

void AxonCollisionSphere::setSmallSphereSize(float size){
    small_sphere_distance = size;
}

void AxonCollisionSphere::push_index(unsigned int element)
{
        collision_list->push_back(element);
        list_size++;
}

/************************  Neuron Collision Implementation  ******************/

NeuronCollisionSphere::NeuronCollisionSphere():small_sphere_list_end(0),big_sphere_list_end(0)
{
}

void NeuronCollisionSphere::pushToSmallSphere(unsigned i)
{
    //If i is already inside the "good" side we do nothing
    if(i < small_sphere_list_end || small_sphere_list_end == collision_list->size()){
        return;
    }

    unsigned jkr = collision_list->at(i);
    collision_list->at(i) = collision_list->at(small_sphere_list_end);
    collision_list->at(small_sphere_list_end) = jkr;
    small_sphere_list_end++;

    // WARNING small sphere size should never be greater than the big one.
    if(small_sphere_list_end > big_sphere_list_end)
        big_sphere_list_end = small_sphere_list_end;
}

void NeuronCollisionSphere::pushToBigSphere(unsigned i)
{
    //If i is already on the "other" side we do nothing
    if(i < big_sphere_list_end || big_sphere_list_end == collision_list->size()){
        return;
    }

    unsigned jkr = collision_list->at(i);
    collision_list->at(i)=collision_list->at(big_sphere_list_end);
    collision_list->at(big_sphere_list_end) = jkr;
    big_sphere_list_end++;

}


void NeuronCollisionSphere::popFromSmallSphere(unsigned i)
{
    //If i is already on the "other" side we do nothing
    if(i >= small_sphere_list_end || small_sphere_list_end == 0){
        return;
    }

    unsigned jkr = collision_list->at(i);
    collision_list->at(i)=collision_list->at(small_sphere_list_end-1);
    collision_list->at(small_sphere_list_end-1) = jkr;
    small_sphere_list_end--;
}

void NeuronCollisionSphere::popFromBigSphere(unsigned i)
{
    //If i is already on the "other" side we do nothing
    if(i >= big_sphere_list_end || big_sphere_list_end == 0){
        return;
    }

    unsigned jkr = collision_list->at(i);
    collision_list->at(i)=collision_list->at(big_sphere_list_end-1);
    collision_list->at(big_sphere_list_end-1) = jkr;
    big_sphere_list_end--;

    // WARNING small sphere size should never be greater than the big one.
    if(big_sphere_list_end < small_sphere_list_end)
        small_sphere_list_end = big_sphere_list_end;
}


void NeuronCollisionSphere::setBigSphereSize(float size){
    big_sphere_distance = size;
}

void NeuronCollisionSphere::setSmallSphereSize(float size){
    small_sphere_distance = size;
}

void NeuronCollisionSphere::push_index(unsigned int element)
{
        collision_list->push_back(element);
        list_size++;
}

/************************  PLY Sphere Collision Implementation  ******************/

PLYCollisionSphere::PLYCollisionSphere()
{

}

void PLYCollisionSphere::pushToSmallSphere(unsigned i, unsigned t)
{
    //Last position of the indexes list.
    unsigned list_end = small_sphere_list_end[i];
    //If i is already inside the "good" side we do nothing
    if(t < list_end || list_end == collision_list->at(i).size()){
        return;
    }

    unsigned jkr = collision_list->at(i)[t];
    collision_list->at(i)[t] = collision_list->at(i)[list_end];
    collision_list->at(i)[list_end] = jkr;
    small_sphere_list_end[i]++;

    // WARNING small sphere size should never be greater than the big one.
    if(list_end > big_sphere_list_end[i])
        big_sphere_list_end[i] = small_sphere_list_end[i];
}

void PLYCollisionSphere::pushToBigSphere(unsigned i, unsigned t)
{
    //Last position of the indexes list.
    unsigned list_end = big_sphere_list_end[i];

    //If t is already on the "other" side we do nothing
    if(t < list_end || list_end == collision_list->size()){
        return;
    }

    unsigned jkr = collision_list->at(i)[t];
    collision_list->at(i)[t]=collision_list->at(i)[list_end];
    collision_list->at(i)[list_end] = jkr;
    big_sphere_list_end[i]++;
}



void PLYCollisionSphere::popFromSmallSphere(unsigned i, unsigned t)
{
    //Last position of the indexes list.
    unsigned list_end = small_sphere_list_end[i];

    // if out index t is already on the "correct" side of the list
    if(t < list_end || list_end == collision_list->at(i).size()){
        return;
    }

    unsigned jkr = collision_list->at(i)[t];
    collision_list->at(i)[t] = collision_list->at(i)[list_end-1];
    collision_list->at(i)[list_end-1] = jkr;
    small_sphere_list_end[i]++;

    // WARNING small sphere size should never be greater than the big one.
    if(small_sphere_list_end[i] > big_sphere_list_end[i])
        big_sphere_list_end[i] = small_sphere_list_end[i];
}

void PLYCollisionSphere::popFromBigSphere(unsigned i, unsigned t)
{
    //Last position of the indexes list.
    unsigned list_end = big_sphere_list_end[i];
    //If i is already on the "other" side we do nothing
    if(t >= list_end || list_end == 0){
        return;
    }

    unsigned jkr = collision_list->at(i)[t];
    collision_list->at(i)[t] = collision_list->at(i)[list_end-1];
    collision_list->at(i)[list_end-1] = jkr;
    big_sphere_list_end[i]--;

    // WARNING small sphere size should never be greater than the big one.
    if(big_sphere_list_end < small_sphere_list_end)
        small_sphere_list_end = big_sphere_list_end;
}


void PLYCollisionSphere::setBigSphereSize(float size){
    big_sphere_distance = size;
}

void PLYCollisionSphere::setSmallSphereSize(float size){
    small_sphere_distance = size;
}

void PLYCollisionSphere::push_ply(std::vector<unsigned> list)
{
    collision_list->push_back(list);
    list_size++;
}



/*******************************************   Sphere Sphere Collision Implementation ***********************************************/


SphereCollisionSphere::SphereCollisionSphere():small_sphere_list_end(0),big_sphere_list_end(0)
{
}

void SphereCollisionSphere::pushToSmallSphere(unsigned i)
{
    //If i is already inside the "good" side we do nothing
    if(i < small_sphere_list_end || small_sphere_list_end == collision_list->size()){
        return;
    }

    unsigned jkr = collision_list->at(i);
    collision_list->at(i) = collision_list->at(small_sphere_list_end);
    collision_list->at(small_sphere_list_end) = jkr;
    small_sphere_list_end++;

    // WARNING small sphere size should never be greater than the big one.
    if(small_sphere_list_end > big_sphere_list_end)
        big_sphere_list_end = small_sphere_list_end;
}

void SphereCollisionSphere::pushToBigSphere(unsigned i)
{
    //If i is already on the "other" side we do nothing
    if(i < big_sphere_list_end || big_sphere_list_end == collision_list->size()){
        return;
    }

    unsigned jkr = collision_list->at(i);
    collision_list->at(i)=collision_list->at(big_sphere_list_end);
    collision_list->at(big_sphere_list_end) = jkr;
    big_sphere_list_end++;

}


void SphereCollisionSphere::popFromSmallSphere(unsigned i)
{
    //If i is already on the "other" side we do nothing
    if(i >= small_sphere_list_end || small_sphere_list_end == 0){
        return;
    }

    unsigned jkr = collision_list->at(i);
    collision_list->at(i)=collision_list->at(small_sphere_list_end-1);
    collision_list->at(small_sphere_list_end-1) = jkr;
    small_sphere_list_end--;
}

void SphereCollisionSphere::popFromBigSphere(unsigned i)
{
    //If i is already on the "other" side we do nothing
    if(i >= big_sphere_list_end || big_sphere_list_end == 0){
        return;
    }

    unsigned jkr = collision_list->at(i);
    collision_list->at(i)=collision_list->at(big_sphere_list_end-1);
    collision_list->at(big_sphere_list_end-1) = jkr;
    big_sphere_list_end--;

    // WARNING small sphere size should never be greater than the big one.
    if(big_sphere_list_end < small_sphere_list_end)
        small_sphere_list_end = big_sphere_list_end;
}


void SphereCollisionSphere::setBigSphereSize(float size){
    big_sphere_distance = size;
}

void SphereCollisionSphere::setSmallSphereSize(float size){
    small_sphere_distance = size;
}

void SphereCollisionSphere::push_index(unsigned int element)
{
        collision_list->push_back(element);
        list_size++;
}







