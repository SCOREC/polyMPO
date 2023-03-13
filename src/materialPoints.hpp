#ifndef POLYMPMTEST_MATERIALPOINTS_H
#define POLYMPMTEST_MATERIALPOINTS_H

#include "utils.hpp"

namespace polyMpmTest{

class MaterialPoints {
  private:
    int count_;
    Vector2View positioins_;
    IntView elmIDs_;
    BoolView isActive_;

  public:
    MaterialPoints(){};
    MaterialPoints( int count,
                    Vector2View positions,
                    IntView elmIDs, 
                    BoolView isActive):
                    count_(count),
                    positions_(positions),
                    elmIDs_(elmIDs),
                    isActive_(isActive){};
 

    int getCount() { return count_; }
    Vector2View getPositions() {}
    IntView getElmsIDs();
    BoolView isActive();

//TODO:MUTATOR  
    void T2LTracking(Vector2View dx);
    
};

}
#endif

