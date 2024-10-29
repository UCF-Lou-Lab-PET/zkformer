#include <iostream>

#include <libsnark/gadgetlib1/gadget.hpp>
#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>
#include <libff/algebra/curves/public_params.hpp>
#include "libff/algebra/fields/field_utils.hpp"

namespace libsnark {


template <typename FieldT>
class conv_gadget : public gadget<FieldT>
{
private:
    //intermediate product values of x00*a00
    // head + ys + tail = all_sum_products
    pb_variable_array<FieldT> sum_products;
    size_t output_h;
    size_t output_w;
    size_t input_h;
    size_t input_w;
    size_t kernel_h;
    size_t kernel_w;

public:
    pb_tensor<FieldT> input;
    pb_tensor<FieldT> kernel;
    pb_tensor<FieldT> output;
    conv_gadget(){};
    conv_gadget(protoboard<FieldT> &pb,
                pb_tensor<FieldT> &input,
                pb_tensor<FieldT> &kernel,
                pb_tensor<FieldT> &output,
                const std::string &annotation_prefix = "") : 
                gadget<FieldT>(pb, annotation_prefix), input(input), kernel(kernel), output(output)
    {
        input_h = input.shape[1];
        input_w = input.shape[2];
        kernel_h = kernel.shape[1];
        kernel_w = kernel.shape[2];
        output_h = input_h - kernel_h + 1;
        output_w = input_w - kernel_w + 1;
        // the sum contains y, it's not desired to re-define a lot of replicated values...
        // following == in + k - 1
        size_t intermediate_len = (input_h + kernel_h - 1) * (input_w + kernel_w -1)  - output_h * output_w;
        sum_products.allocate(pb, intermediate_len, "intermeidiate sum_products");
    };

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};


template <typename FieldT>
void conv_gadget<FieldT>::generate_r1cs_constraints()
{

    FieldT s = FieldT(2);
    FieldT t = FieldT(3);

    linear_combination<FieldT> A, B, C;

    FieldT temp_s = FieldT::one();
    FieldT temp_t = FieldT::one();

    for(size_t x_i=0; x_i<input_h;x_i++)
    {
        for(size_t x_j=0; x_j<(input_h);x_j++){
            A.add_term(input[(input_w*x_i)+x_j], temp_s*temp_t); 

            temp_t *= t;
        }
        temp_t = FieldT::one();
        temp_s *= s;
    }
    
    temp_s = FieldT::one();
    temp_t = FieldT::one();
    for(size_t a_i=kernel_h; a_i>0; --a_i)
    {
        for(size_t a_j=kernel_w; a_j>0; --a_j){
            B.add_term(kernel[(a_i-1)*kernel_w+a_j-1], temp_s*temp_t);
            temp_t *= t;
        }
        temp_t = FieldT::one();
        temp_s *= s;
    }

    
    temp_s = FieldT::one();
    temp_t = FieldT::one();
    size_t intermediate_index=0;
    size_t y_index=0;
    for(size_t y_i=0; y_i<(input_h+kernel_h-1);y_i++){
        for(size_t y_j=0; y_j<(input_w+kernel_w-1);y_j++){
            if(y_i < kernel_h - 1 || y_i > kernel_h + output_h - 2){
                C.add_term(sum_products[intermediate_index++], temp_s*temp_t);
            }else if (y_j < kernel_w - 1 || y_j > kernel_w + output_w - 2){
                C.add_term(sum_products[intermediate_index++], temp_s*temp_t);
            }else{
                C.add_term(output[y_index++], temp_s*temp_t);
            }
            temp_t *= t;
        }
        temp_t = FieldT::one();
        temp_s *= s;
    }
    this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(A, B, C), "a convolution");
};

template <typename FieldT>
void conv_gadget<FieldT>::generate_r1cs_witness()
{
    size_t intermediate_index=0;
    size_t y_index=0;
    for(size_t y_i=0; y_i<(input_h+kernel_h-1); y_i++){
        for(size_t y_j=0; y_j<(input_w+kernel_w-1); y_j++){
            FieldT y = FieldT::zero();
            for(size_t x_i=0; x_i<input_h; x_i++){
                for(size_t x_j=0; x_j<input_w; x_j++){
                    for(size_t a_i=0; a_i<kernel_h; a_i++){
                        for(size_t a_j=0; a_j<kernel_w; a_j++){
                            if((x_i+ (kernel_h-a_i-1)) == y_i && (x_j+(kernel_w-a_j-1)) == y_j)
                            {
                                y += this->pb.val(input[x_i*input_w+x_j]) * this->pb.val(kernel[a_i*kernel_w+a_j]);
                            }
                        }
                    }
                }
            }
            if(y_i < kernel_h - 1 || y_i > kernel_h + output_h - 2){
                this->pb.val(sum_products[intermediate_index++]) = y;
            }else if (y_j < kernel_w - 1 || y_j > kernel_w + output_w - 2){
                this->pb.val(sum_products[intermediate_index++]) = y;
            }else{
                this->pb.val(output[y_index++]) = y;
            }
        }
    }

};



template <typename FieldT>
class conv_3d_gadget : public gadget<FieldT>
{
private:
    std::vector<conv_gadget<FieldT>> convs;
    // num of 2D kernels == out_dim * in_dim
    std::vector<pb_tensor<FieldT>> kernels;
    std::vector<pb_tensor<FieldT>> intermediate_outputs;
    // output dim == kernel number
    size_t output_dim;
    size_t input_dim;
    // spatial information
    size_t kernel_size;
    // padding ensures the same spatial size of input and output
    // future implementation will extend to stride
    size_t padding_size;
    size_t input_h;
    size_t input_w;

public:
    pb_tensor<FieldT> input;
    pb_tensor<FieldT> output;
    conv_3d_gadget(){};
    conv_3d_gadget(protoboard<FieldT> &pb,
                pb_tensor<FieldT> &input,
                pb_tensor<FieldT> &output,
                size_t kernel_size,
                const std::string &annotation_prefix = "") : 
                gadget<FieldT>(pb, annotation_prefix), input(input), output(output), kernel_size(kernel_size)
    {
        input_dim = input.shape[0];
        output_dim = output.shape[0];
        input_h = input.shape[1];
        input_w = input.shape[2];
        // with kernel_size==3, we have the padding==1, ensuring same padding
        padding_size = (kernel_size - 1) / 2;

        // std::cout<<" padding: " <<padding_size<<std::endl;
        // kernel i means which output dimension
        for(size_t kernel_i=0; kernel_i < output_dim; kernel_i++){
            for(size_t feature_i=0; feature_i < input_dim; feature_i++){
                pb_tensor<FieldT> input_t;
                input_t.shape = {1, input.shape[1] + 2*padding_size, input.shape[2] + 2*padding_size};
                // pad a 2d feature into input_t
                size_t channel_offset = feature_i*input_h*input_w;
                for(size_t i=0; i<input_h+padding_size*2; i++){
                    if(i < padding_size || i > padding_size + input_h - 1){
                        input_t.insert(input_t.end(), input_w + 2*padding_size, 0);
                    }else{
                        input_t.insert(input_t.end(), padding_size, 0);
                        input_t.insert(input_t.end(), input.begin() + channel_offset + (i-padding_size)*input_w, input.begin()+ channel_offset + (i-padding_size+1)*input_w);
                        input_t.insert(input_t.end(), padding_size, 0);
                    }
                }
                size_t kernel_index = kernel_i * input_dim + feature_i;
                kernels.push_back(pb_tensor<FieldT>(kernel_size*kernel_size, 0));
                kernels[kernel_index].shape = {1, kernel_size, kernel_size};
                kernels[kernel_index].allocate(pb, kernel_size*kernel_size, " convolution weights ");
                
                intermediate_outputs.push_back(pb_tensor<FieldT>(input_h*input_w, 0));
                intermediate_outputs[kernel_index].shape = {1, input_h, input_w};
                intermediate_outputs[kernel_index].allocate(pb, input_h*input_w, " intermediate outputs on channel ");

                convs.push_back(conv_gadget<FieldT>(pb, input_t, kernels[kernel_index], intermediate_outputs[kernel_index], " one 2d conv in conv3d "));
            }
        }
        
    };

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

template <typename FieldT>
void conv_3d_gadget<FieldT>::generate_r1cs_constraints()
{
    for(size_t kernel_i=0; kernel_i < output_dim; kernel_i++){
        for(size_t feature_i=0; feature_i < input_dim; feature_i++){
            size_t conv_index = kernel_i * input_dim + feature_i;
            convs[conv_index].generate_r1cs_constraints();
        }
    }
    // sum the intermediate outputs to final outputs
}


template <typename FieldT>
void conv_3d_gadget<FieldT>::generate_r1cs_witness()
{
    //gen output kernels
    for(size_t kernel_i=0; kernel_i < output_dim; kernel_i++){
        for(size_t feature_i=0; feature_i < input_dim; feature_i++){
            size_t kernel_index = kernel_i * input_dim + feature_i;
            for(size_t i=0; i<kernel_size; i++){
                for(size_t j=0; j<kernel_size; j++){
                    size_t pos_index = i*kernel_size+j;
                    // this->pb.val(kernels[kernel_index][pos_index]) = FieldT((j+1) % 64);
                    this->pb.val(kernels[kernel_index][pos_index]) = FieldT(kernel_i+1);
                }
            }
        }
    }

    for(size_t kernel_i=0; kernel_i < output_dim; kernel_i++){
        for(size_t feature_i=0; feature_i < input_dim; feature_i++){
            size_t conv_index = kernel_i * input_dim + feature_i;
            convs[conv_index].generate_r1cs_witness();
        }
    }

    // sum the intermediate outputs to final outputs
    for(size_t kernel_i=0; kernel_i < output_dim; kernel_i++){
        // across input dimension
        for(size_t feature_i=0; feature_i < input_dim; feature_i++){
            size_t kernel_index = kernel_i * input_dim + feature_i;
            for(size_t i=0; i<input_h; i++){
                for(size_t j=0; j<input_w; j++){
                    size_t pos_index = i*input_w+j;
                    this->pb.val(output[kernel_i*input_h*input_w +  pos_index]) += this->pb.val(intermediate_outputs[kernel_index][pos_index]);
                }
            }
        }
    }

}

}