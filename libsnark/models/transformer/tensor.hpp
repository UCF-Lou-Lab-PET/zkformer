#include <iostream>

#include <libff/common/utils.hpp>
#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>

namespace libsnark {

template <typename FieldT>
// in this impl, tensor is a vector of size 
// << tensor, is to print every dim
// tensor doesn't hold numeric values for a network, instead, the pb varaibles for a circuit
class pb_tensor : public pb_variable_array<FieldT> {
public:
    std::vector<size_t> shape;
    pb_tensor() :  pb_variable_array<FieldT>() {
    }

    pb_tensor(size_t size, const pb_variable<FieldT> &value) :  pb_variable_array<FieldT>(size, value) {
        shape.clear();
        shape.push_back(size);
    }

    friend std::ostream &operator<<(std::ostream &os, const pb_tensor<FieldT> &vec) {
        size_t cnt = 1;
        os << "[" << vec.size() << "]";
        os << "[";
        for (size_t i = 0; i < vec.shape.size(); ++i) {
            if (i)
                os << ", ";
            os << vec.shape[i];
            cnt *= vec.shape[i];
        }
        os << "]";
        // assert(cnt == vec.size());
        return os;
    }

};

}