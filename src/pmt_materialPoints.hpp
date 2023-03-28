#ifndef POLYMPMTEST_MATERIALPOINTS_H
#define POLYMPMTEST_MATERIALPOINTS_H

#include "pmt_utils.hpp"

namespace polyMpmTest{

class MaterialPoints {
  private:
    int count_;
    Vector2View positions_;
    BoolView isActive_;

  public:
    MaterialPoints(){};
    MaterialPoints( int count,
                    Vector2View positions,
                    BoolView isActive):
                    count_(count),
                    positions_(positions),
                    isActive_(isActive){};

    int getCount() { return count_; }
    Vector2View getPositions() { return positions_; }
    BoolView isActive() { return isActive_; }
};

}
#endif

