#ifndef SYNC_BACK_WATCHER_PY_WEDGE_HPP
#define SYNC_BACK_WATCHER_PY_WEDGE_HPP

#include <vec3.h>

template <typename Type>
class Wedge{
public:
    Wedge() = default;

    Wedge(Type angle_inp, Type x_0_inp, Type x_1_inp):x_0_(x_0_inp), x_1_(x_1_inp){
        full_openning_angle_=angle_inp;
        b_ = tan(angle_inp / 2.);
        hallf_max_width_ = Y_Plus(x_1_);
        side_length_ = sqrt(pow((x_1_-x_0_),2.) + pow(hallf_max_width_,2.));
    }

    Type full_openning_angle()const{return full_openning_angle_;}
    Type x_1()const{return x_1_;}
    Type x_0()const{return x_0_;}
    Type hallf_max_width()const{return hallf_max_width_;}
    Type b()const{return b_;}

    Type  Y_Plus(Type x){
        Type y_plus = b_ * (x - x_0_);
        return y_plus;
    }
    Type Y_Min(Type x){
        Type y_min = -b_ * (x - x_0_);
        return y_min;
    }

    void move_wedge(Type movement){
        x_1_ = x_1_ + movement;
        x_0_ = x_0_ + movement;
    }
    bool is_in_range(vec3<Type> const& pos){
        return (pos[0]>x_0_)&&(pos[0]<x_1_);
    }

    bool is_intersecting(vec3<Type> const& pos)
    {
        if (is_in_range(pos) && (((pos[1]>Y_Min(pos[0])) && (pos[1]<Y_Plus(pos[0]))))) {
            return true;
        }
        else {
            return false;
        }
    }

    Type shortest_distance(vec3<Type> const& pos){

        if(pos[0]<x_0_){
            return sqrt(pow(pos[0]-x_0_,2)+pow(pos[1],2.));
        }
        else if(is_in_range(pos)){
            Type sign = pos[1]>=0?1:-1;
            return (-(x_1_-x_0_)*pos[1]-(x_0_-pos[0])*sign*hallf_max_width_)/side_length_;
        }
        else{
            Type sign = pos[1]>=0?1:-1;
            return sqrt(pow(pos[0]-x_1_,2)+pow(pos[1]-sign*hallf_max_width_,2)); //Todo this is an approximation
        }
    }
private:
    Type full_openning_angle_;
    Type x_0_;
    Type x_1_;
    Type b_;
    Type side_length_, hallf_max_width_;
};

#endif //SYNC_BACK_WATCHER_PY_WEDGE_HPP
