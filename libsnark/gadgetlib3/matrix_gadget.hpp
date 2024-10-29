#include <iostream>

#include <libsnark/gadgetlib1/gadget.hpp>
#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>
#include <libff/algebra/curves/public_params.hpp>
#include "libff/algebra/fields/field_utils.hpp"
#include <libsnark/models/transformer/tensor.hpp>

namespace libsnark {

/**
 * construct a circuit for matrix multiplication
 * matrix_gadget(pb, mat_1, mat_2, mat_result, "")
*/
template <typename FieldT>
class matrix_gadget : public gadget<FieldT>
{
private:
    //intermediate product values of x00*a00
    pb_variable_array<FieldT> products;
    size_t output_h;
    size_t output_w;

public:
    pb_variable_array<FieldT> input;
    size_t input_h;
    size_t input_w;
    pb_variable_array<FieldT> kernel;
    size_t kernel_h;
    size_t kernel_w;
    pb_variable_array<FieldT> output;
    matrix_gadget(){};
    matrix_gadget(protoboard<FieldT> &pb,
                pb_variable_array<FieldT> &input,
                size_t i_h,
                size_t i_w,
                pb_variable_array<FieldT> &kernel,
                size_t k_h,
                size_t k_w,
                pb_variable_array<FieldT> &output,
                const std::string &annotation_prefix = "") : 
                gadget<FieldT>(pb, annotation_prefix), input(input), input_h(i_h), input_w(i_w), kernel(kernel), kernel_h(k_h), kernel_w(k_w), output(output)
    {
        output_h = input_h;
        output_w = kernel_w;
        // this impl uses prefix sum
        products.allocate(pb, output_h*output_w*(input_w-1), "intermeidiate products");
    };
    // construct from pb_tensor
    matrix_gadget(protoboard<FieldT> &pb,
                pb_tensor<FieldT> &input,
                pb_tensor<FieldT> &kernel,
                pb_tensor<FieldT> &output,
                const std::string &annotation_prefix = "") : 
                gadget<FieldT>(pb, annotation_prefix), input(input), kernel(kernel), output(output)
    {
        if(input.shape[0] == 1){
            input_h = input.shape[1];
            input_w = input.shape[2];
        }else{
            input_h = input.shape[0];
            input_w = input.shape[1] * input.shape[2];
        }
        if(kernel.shape[0] == 1){
            kernel_h = kernel.shape[1];
            kernel_w = kernel.shape[2];
        }else{
            kernel_h = kernel.shape[0];
            kernel_w = kernel.shape[1] * kernel.shape[2];
        }
        // output has c==1
        output_h = input_h;
        output_w = kernel_w;
        products.allocate(pb, output_h*output_w*(input_w-1), FMT(this->annotation_prefix, " intermeidiate products"));
    };

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};


template <typename FieldT>
void matrix_gadget<FieldT>::generate_r1cs_constraints()
{
    FieldT s = FieldT(2);
    FieldT t = qpow(s, kernel_w);
    for (size_t i = 0; i < input_w; ++i)
    {
        // FieldT t = qpow(s, kernel_w);
        linear_combination<FieldT> A, B, C;

        FieldT temp_s = FieldT::one();
        FieldT temp_t = FieldT::one();
        if(input_h * input_w > kernel_h * kernel_w){
            for(size_t x_i = 0; x_i < input_h; x_i ++){
                A.add_term(input[(x_i*input_w)+i], temp_t);
                temp_t *= t;
            }
            for(size_t a_j = 0; a_j < kernel_w; a_j ++){
                B.add_term(kernel[(i*kernel_w)+a_j], temp_s);
                temp_s *= s;
            }
        }else{
            for(size_t x_i = 0; x_i < input_h; x_i ++){
                B.add_term(input[(x_i*input_w)+i], temp_t);
                temp_t *= t;
            }
            for(size_t a_j = 0; a_j < kernel_w; a_j ++){
                A.add_term(kernel[(i*kernel_w)+a_j], temp_s);
                temp_s *= s;
            }
        }
        temp_s = FieldT::one();
        temp_t = FieldT::one();
        FieldT temp_minus_one = FieldT(-1);
        // std::cout<<"add for C"<<std::endl;
        for(size_t x_i=0; x_i<input_h; x_i++){
            for(size_t a_j=0; a_j<kernel_w; a_j++){
                // /sum_i x[x_i,i] * a[i,a_j] = y[x_i, a_j]
                if(i == input_w - 1){
                    C.add_term( output[((x_i*output_w)+a_j)],temp_s*temp_t);
                }else{
                    C.add_term( products[((x_i*output_w)+a_j)*(input_w-1)+i],temp_s*temp_t);
                }
                if(i!=0){
                    C.add_term( products[((x_i*output_w)+a_j)*(input_w-1)+i-1],temp_minus_one*temp_s*temp_t);
                }
                temp_s *= s;
            }
            temp_s = FieldT::one();
            temp_t *= t;
        }
        // std::cout<<std::endl;
        this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(A, B, C), "matrix multiplication");
        // std::cout<<std::endl;
    }
};

/*
compute the value of intermediate variables based on given in/out
specifically, products
*/
template <typename FieldT>
void matrix_gadget<FieldT>::generate_r1cs_witness()
{
    for(size_t y_i=0; y_i<output_h; y_i++){
        for(size_t y_j=0; y_j<output_w; y_j++){
            FieldT y = FieldT::zero();
            // pos: 1~n
            for(size_t pos=0; pos<input_w; pos++){
                FieldT temp_y = FieldT::zero();
                temp_y = this->pb.val(input[y_i*input_w + pos]) * this->pb.val(kernel[pos*kernel_w + y_j]);
                y += temp_y;
                if(pos == input_w - 1){
                    this->pb.val(output[y_i*output_w + y_j]) = y;
                }else{
                    this->pb.val(products[(y_i*output_w+y_j)*(input_w-1)+pos]) = y;
                }
            }
        }
        // std::cout<<std::endl;
    }
   
};

}