#ifndef POLYMPMTEST_MATERIALPOINTS_H
#define POLYMPMTEST_MATERIALPOINTS_H


namespace polyMpmTest{

#include "utils.hpp"

class MaterialPoints {
  private:
    int count_;
    VectorView positioins_;
    IntView elmIDs_;
    BoolView isActive_;

  public:
    MaterialPoints();
    MaterialPoints(int count, VectorView position, IntView elmIDs, BoolView isActive);
 

    int getCount();
    VectorView getPositions();
    IntView getElmsIDs();
    BoolView isActive();

//TODO:MUTATOR
};

}
#endif

