#include <iostream>

class MtrCtrl{
    private:
        double target_speed_;

    public:
        MtrCtrl(): target_speed_(0){
            std::cout << "[C++] MotorController Created." << std::endl;
        }
        ~MtrCtrl() {
            std::cout << "[C++] MotorController Destroyed." << std::endl;
        }
        void set_speed(double speed) {
            target_speed_ = speed;
        }

        double get_speed() const {
            return target_speed_;
        }
};

extern "C"{
    void* MtrCtrl_new(){
        return new MtrCtrl();
    }

    void MtrCtrl_delete(void* obj){
        if(obj != nullptr)
            delete static_cast<MtrCtrl*>(obj);
    }

    void MtrCtrl_set_speed(void* obj, double speed){
        static_cast<MtrCtrl*>(obj)->set_speed(speed);
    }
    double MtrCtrl_get_speed(void* obj){
        return static_cast<MtrCtrl*>(obj)->get_speed();
    }
}
