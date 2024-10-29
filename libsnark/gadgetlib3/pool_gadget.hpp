#include <iostream>

#include <libsnark/gadgetlib1/gadget.hpp>
#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>
#include <libff/algebra/curves/public_params.hpp>
#include "libff/algebra/fields/field_utils.hpp"



namespace libsnark {

template <typename FieldT>
class pool_gadget : public gadget<FieldT>
{
private:
    //intermediate product values of x00*a00
    size_t output_h;
    size_t output_w;
    size_t kernel_size;

public:
    pb_variable_array<FieldT> input;
    size_t input_h;
    size_t input_w;
    pb_variable_array<FieldT> output;
    size_t kernel_h;
    size_t kernel_w;
    pool_gadget(){};
    pool_gadget(protoboard<FieldT> &pb,
                pb_variable_array<FieldT> &input,
                size_t i_h,
                size_t i_w,
                pb_variable_array<FieldT> &output,
                size_t k_h,
                size_t k_w,
                const std::string &annotation_prefix = "") : 
                gadget<FieldT>(pb, annotation_prefix), input(input), input_h(i_h), input_w(i_w), output(output), kernel_h(k_h), kernel_w(k_w)
    {
        output_h = input_h / kernel_h;
        output_w = input_w / kernel_w;
        kernel_size = kernel_h * kernel_w;
    };

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

template <typename FieldT>
void pool_gadget<FieldT>::generate_r1cs_constraints()
{
    for(size_t y_i=0; y_i<output_h; y_i++){
        for(size_t y_j=0; y_j<output_w; y_j++){
            linear_combination<FieldT> A, B, C;
            for(size_t k_i=0; k_i<kernel_h; k_i++)
                for(size_t k_j=0; k_j<kernel_w; k_j++){
                    // for the top_left value in current window
                    // x_i = y_i*kernel_h, x_j= y_j*kernel_w
                    size_t x_index = (y_i*kernel_h + k_i)*input_w + (y_j*kernel_w+k_j);
                    A.add_term(input[x_index], 1);
                }
            B.add_term(0, 1);
            C.add_term(output[y_i*output_w + y_j], FieldT(kernel_size));
            this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(A, B, C),"pooling");
        }
    }
}


template <typename FieldT>
void pool_gadget<FieldT>::generate_r1cs_witness()
{
    FieldT divisor = FieldT(kernel_size).inverse();
    for(size_t y_i=0; y_i<output_h; y_i++){
        for(size_t y_j=0; y_j<output_w; y_j++){
            FieldT temp_y = FieldT::zero();
            for(size_t k_i=0; k_i<kernel_h; k_i++)
                for(size_t k_j=0; k_j<kernel_w; k_j++){
                    // for the top_left value in current window
                    // x_i = y_i*kernel_h, x_j= y_j*kernel_w
                    size_t x_index = (y_i*kernel_h + k_i)*input_w + (y_j*kernel_w+k_j);
                    temp_y += this->pb.val(input[x_index]);
                }
            this->pb.val(output[y_i * output_w + y_j]) = temp_y*divisor;
        }
    }
}

/**
 * a pool with stride and padding
 * for poolformer, where pooling is used as a token mixer
*/
template <typename FieldT>
class avg_pool_gadget : public gadget<FieldT>
{
private:
    pb_tensor<FieldT> padded_input;
    size_t kernel_size;

public:
    pb_tensor<FieldT> input;
    pb_tensor<FieldT> output;
    size_t kernel_h;
    size_t kernel_w;
    size_t stride;
    size_t padding;
    avg_pool_gadget(){};
    avg_pool_gadget(protoboard<FieldT> &pb,
                pb_tensor<FieldT> &input,
                pb_tensor<FieldT> &output,
                size_t kernel_h,
                size_t kernel_w,
                size_t stride,
                size_t padding,
                const std::string &annotation_prefix = "") : 
                gadget<FieldT>(pb, annotation_prefix), input(input), output(output), kernel_h(kernel_h), kernel_w(kernel_w), stride(stride), padding(padding)
    {
        //specify the private values
        // padded_input = new pb_tensor<FieldT> ((input.shape[1] + 2*padding)*(input.shape[2] + 2*padding), 0);
        pb_tensor<FieldT> padded((input.shape[1] + 2*padding)*(input.shape[2] + 2*padding), 0);
        padded_input = padded;
        padded_input.allocate(this->pb, padded_input.size(), "pooling padded input");
        // padded_input.allocate(this->pb, (input.shape[1] + 2*padding)*(input.shape[2] + 2*padding), "pooling padded input");
        pad(input, padded_input, padding);
        kernel_size = kernel_h * kernel_w;
    };

    void pad(pb_tensor<FieldT> &input, pb_tensor<FieldT> &padded_input, size_t padding);
    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

// input = C,H,W
template <typename FieldT>
void avg_pool_gadget<FieldT>::pad(pb_tensor<FieldT> &input, pb_tensor<FieldT> &padded_input, size_t padding)
{
    size_t height = input.shape[1];
    size_t width = input.shape[2];
    size_t new_h = height + 2*padding;
    size_t new_w = width + 2*padding;
    for (size_t i = 0; i < height; i++)
        for (size_t j = 0; j < width; j++){
            // std::cout <<"padded: [" << i + padding <<","<< j + padding <<"] == input[" <<i <<"," <<j<<"]"<<std::endl;
            padded_input[(i + padding)*new_w + (j + padding)] = input[i*width + j];
        }
    padded_input.shape = {input.shape[0], new_h, new_w};
    // std::cout <<"padded: " << padded_input <<std::endl;
};

// gen constraints with padded input
template <typename FieldT>
void avg_pool_gadget<FieldT>::generate_r1cs_constraints()
{
    size_t height = padded_input.shape[1];
    size_t width = padded_input.shape[2];
    size_t output_h = (height - kernel_h) / stride + 1;
    size_t output_w = (width - kernel_w) / stride + 1;
    for(size_t y_i=0; y_i<output_h; y_i++){
        for(size_t y_j=0; y_j<output_w; y_j++){
            linear_combination<FieldT> A, B, C;
            for (size_t k = 0; k < kernel_h; k++) {
                for (size_t l = 0; l < kernel_w; l++) {
                    // std::cout <<"+[" << y_i * stride + k << ", " <<y_j * stride + l <<"] ";
                    A.add_term(padded_input[(y_i * stride + k)*width + y_j * stride + l], 1);
                }
            }
            B.add_term(0, 1);
            C.add_term(output[y_i*output_w + y_j], FieldT(kernel_size));
            // C.add_term(output[y_i*output_w + y_j], 1);
            this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(A, B, C),"pooling");
        }
    }
};

template <typename FieldT>
void avg_pool_gadget<FieldT>::generate_r1cs_witness()
{
    
    size_t height = padded_input.shape[1];
    size_t width = padded_input.shape[2];

    // print padded
    // for(size_t i=0; i<height; i++){
    //     for(size_t j=0; j <width; j++){
    //         std::cout <<this->pb.val(padded_input[i*width + j]).as_ulong()<<" ";
    //     }
    //     std::cout<<std::endl;
    // }
    
    size_t pooledHeight = (height - kernel_h) / stride + 1;
    size_t pooledWidth = (width - kernel_w) / stride + 1;

    // allocate space for output, this is done outside here

    for (size_t i = 0; i < pooledHeight; i++) {
        for (size_t j = 0; j < pooledWidth; j++) {
            FieldT sum = FieldT::zero();;

            for (size_t k = 0; k < kernel_h; k++) {
                for (size_t l = 0; l < kernel_w; l++) {
                    // std::cout <<"+[" << i * stride + k << ", " <<j * stride + l <<"] " <<"("<< this->pb.val(padded_input[(i * stride + k)*width + j * stride + l]).as_ulong()<<")" <<" ";
                    sum += this->pb.val(padded_input[(i * stride + k)*width + j * stride + l]);
                }
            }
            // pooled[i][j] = sum / (kernel_h * kernel_w);
            this->pb.val(output[i*pooledWidth + j]) = sum * FieldT(kernel_size).inverse();
            // this->pb.val(output[i*pooledWidth + j]) = sum;
        }
    }

};


}