//
// Created by Ales Varabyou on 3/16/19.
//

#include "MinMap.h"

MinMap::MinMap() = default;
MinMap::~MinMap() = default;

MinMap::KmerMap::iterator MinMap::insert(const std::string &key, uint32_t chrID, uint32_t strand) {
    EVec *ev = new EVec(chrID,strand);
    return this->minmap.insert(std::make_pair(key,));
}
