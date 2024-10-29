#include <iostream>

#include <libsnark/gadgetlib1/gadget.hpp>
#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>
#include <libff/algebra/curves/public_params.hpp>
#include "libff/algebra/fields/field_utils.hpp"

#include <libsnark/gadgetlib3/model_gadget.hpp>

// /**
//  * given a network
//  * construct a circuit
//  * for a transformer, building blocks are
//  * -patch embedding: (convolution based)
//  * -transformer block
//  *      -normalization
//  *      -token mixer (skip connect)
//  *      -normalization
//  *      -MLP (skip connect)
// */

namespace libsnark {

/**
 * patch embedding take as input [c, h, w] --> [4*c, h/2, w/2] --> (a kernel of [c', 4*c]) --> ouputs [c', h/2, w/2]
 * patch_embedding(pb, input_tensor, output_tensor"")
 * this can be done with dense connect (matrix multiplication) or convolution
 * don't allocate pb varaibles for unfolded_input
*/
template <typename FieldT>
class patch_embedding_gadget : public gadget<FieldT>
{
private:
    size_t input_feature;
    size_t output_feature;
    pb_tensor<FieldT> unfolded_input;
    pb_tensor<FieldT> kernel;
    std::shared_ptr<matrix_gadget<FieldT>> fc;

public:
    pb_tensor<FieldT> input;
    pb_tensor<FieldT> output;
    size_t embedding_dimension;
    size_t downscaling_factor;
    patch_embedding_gadget(){};
    patch_embedding_gadget(protoboard<FieldT> &pb,
                pb_tensor<FieldT> &input_t,
                pb_tensor<FieldT> &output_t,
                size_t embedding_dimension,
                size_t downscaling_factor,
                const std::string &annotation_prefix = "") : 
                gadget<FieldT>(pb, annotation_prefix), input(input_t), output(output_t), embedding_dimension(embedding_dimension), downscaling_factor(downscaling_factor)
    {
        //specify the private values
        input_feature = input.shape[0]*downscaling_factor*downscaling_factor;
        output_feature = embedding_dimension;

        // output.shape = {output_feature, input.shape[1]/downscaling_factor, input.shape[2]/downscaling_factor};

        unfolded_input.allocate(pb, input.size(), "unfolded input");
        // unfolded_input.shape = {input.shape[0]*downscaling_factor*downscaling_factor, input.shape[1] / downscaling_factor, input.shape[2] / downscaling_factor};

        kernel.allocate(this->pb, input_feature*output_feature, "embedding kernel");
        kernel.shape = {1, output_feature, input_feature};

        // std::cout << "----------test----in-------- "<< std::endl;
        // std::cout << unfolded_input << std::endl;
        // unfold here, somehow makes unfolded_input invisiable in dump, but it holds the correct value.
        unfold(input, unfolded_input, downscaling_factor);
        unfolded_input.shape = {input.shape[0]*downscaling_factor*downscaling_factor, input.shape[1] / downscaling_factor, input.shape[2] / downscaling_factor};
        fc.reset(new matrix_gadget<FieldT>(this->pb, kernel, unfolded_input, output, FMT(this->annotation_prefix, " fc for patch embeddin")));

    };

    void unfold(pb_tensor<FieldT> &input,  pb_tensor<FieldT> &unfolded_input, size_t downscaling_factor);
    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

// unfold and patch-fy
template <typename FieldT>
void patch_embedding_gadget<FieldT>::unfold(pb_tensor<FieldT> &input,  pb_tensor<FieldT> &unfolded_input, size_t downscaling_factor){
    // unfold and reshape input (patchify)
    // std::cout<<"called unfold"<<std::endl;
    // [3, 4, 4]
    auto &vec = input.shape;
    size_t ch_sum = vec[0], h = vec[1], w = vec[2];
    size_t new_h = h / downscaling_factor, new_w = w / downscaling_factor;
    // std::cout <<"--linear input: "<< input << " in_feature: " <<input.shape[0]<< std::endl;
    // [12, 2, 2]
    unfolded_input.shape = {(vec[0] * downscaling_factor * downscaling_factor), new_h, new_w};
    // every channel of original input
    size_t patch_index = 0;
    size_t input_index = 0;
    // [ch_sum, h, w] --> []
    for (size_t ch = 0; ch < ch_sum; ++ch) {
        // (i, j) indicates the location of the folded input
        for (size_t i = 0; i < new_h; ++i) {
            for (size_t j = 0; j < new_w; ++j) {
                // (x, y) indicates position in a patch
                for (size_t x = 0; x < downscaling_factor; ++x) {
                    for (size_t y = 0; y < downscaling_factor; ++y) {
                        // position in the new tensor [12, 2, 2]
                        patch_index = (ch*downscaling_factor*downscaling_factor + x*downscaling_factor + y)*new_h*new_w + i*new_w + j;
                        input_index = ch * h * w + (i * downscaling_factor + x) * w + j * downscaling_factor + y;
                        // std::cout <<"patch[" <<(ch+1)*(x*downscaling_factor+y) <<","<<i <<"," << j<<"]"<< " ";
                        unfolded_input[patch_index] = input[input_index];
                        // unfolded_input.insert(unfolded_input.begin()+patch_index, input.begin()+input_index, input.begin()+input_index+1);
                        // this->pb.val(unfolded_input[patch_index])
                        // = this->pb.val(input[ch * h * w + (i * downscaling_factor + x) * w + j * downscaling_factor + y]);
                    }
                }
            }
        }
    }
};


template <typename FieldT>
void patch_embedding_gadget<FieldT>::generate_r1cs_constraints()
{
    fc->generate_r1cs_constraints();
};

template <typename FieldT>
void patch_embedding_gadget<FieldT>::generate_r1cs_witness()
{
    // unfold and reshape input
    // std::cout <<"--linear temp: "<< tmp << " tmp_feature: " <<tmp.shape.back()<< std::endl;

    // size_t dim_t = unfolded_input.shape[1] * unfolded_input.shape[2];
    // for (size_t i = 0; i<unfolded_input.size(); i += dim_t) {
    //     for (size_t j = 0; j<dim_t; j++) {
    //         std::cout<<this->pb.val(unfolded_input[i+j]).as_ulong() << " ";
    //     }
    //     std::cout<<std::endl;
    // }

    // embed to the target dimension
    // kernel.shape = [1, h, w] = [1, out_feature, in_feature]
    for(size_t i=0; i<kernel.shape[1]; i++){
        for(size_t j=0; j<kernel.shape[2]; j++){
            if((i % kernel.shape[2]) == j)
                this->pb.val(kernel[i*kernel.shape[2] + j]) = FieldT(1);
        }
    }
    fc->generate_r1cs_witness();

    // std::cout << "---------tensor in the patch embedding----------" << std::endl;
    // std::cout <<"input: " << input << std::endl;
    // std::cout <<"patchify: " << unfolded_input << std::endl;
    // std::cout <<"kernel: " << kernel << std::endl;
    // std::cout <<"output"<< output << std::endl;
    // size_t dim2 = output.shape[1] * output.shape[2];
    // for (size_t i = 0; i<output.size(); i += dim2) {
    //     for (size_t j = 0; j<dim2; j++) {
    //         std::cout<<this->pb.val(output[i+j]).as_ulong() << " ";
    //     }
    //     std::cout<<std::endl;
    // }

    // some test code
    // unfolded = [12, 2, 2]
    // size_t dim_t = unfolded_input.shape[1] * unfolded_input.shape[2];
    // for (size_t i = 0; i<unfolded_input.size(); i += dim_t) {
    //     for (size_t j = 0; j<dim_t; j++) {
    //         std::cout<<this->pb.val(unfolded_input[i+j]).as_ulong() << " ";
    //     }
    //     std::cout<<std::endl;
    // }
    // std::cout<<"-------------output-------------"<<std::endl;
    // size_t dim2 = output.shape[1] * output.shape[2];
    // for (size_t i = 0; i<output.size(); i += dim2) {
    //     for (size_t j = 0; j<dim2; j++) {
    //         std::cout<<this->pb.val(output[i+j]).as_ulong() << " ";
    //     }
    //     std::cout<<std::endl;
    // }

    // std::cout <<"--linear output: "<< output << " out_feature: " <<output.shape.back()<< std::endl;
};




/**
 * layer_norm(pb, input_tensor, output_tensor, "")
 * with gamma,beta are learnale (initially 0,1), eps is a constant
 * input [c, h, w] as c feature maps of (h,w), do norm for each map
*/
template <typename FieldT>
class layer_norm_gadget : public gadget<FieldT>
{
private:
    size_t dim;
    size_t feature_height;
    size_t feature_width;
    pb_tensor<FieldT> avg;
    pb_tensor<FieldT> var;
public:
    pb_tensor<FieldT> input;
    pb_tensor<FieldT> output;
    layer_norm_gadget(){};
    layer_norm_gadget(protoboard<FieldT> &pb,
                pb_tensor<FieldT> &input,
                pb_tensor<FieldT> &output,
                const std::string &annotation_prefix = "") : 
                gadget<FieldT>(pb, annotation_prefix), input(input), output(output)
    {
        //specify the private values
        dim = input.shape[0];
        feature_height = input.shape[1];
        feature_width = input.shape[2];
        avg.allocate(pb, dim, "sum over each feature map");
        var.allocate(pb, dim, "var over each feature map");
    };

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

/**
 * the constraints are
 * sum = /sum x_i
 * avg = sum * dim.inverse
 * y_i = x_i - avg / sqrt(var)
 * number of x, y == input.size()
*/
template <typename FieldT>
void layer_norm_gadget<FieldT>::generate_r1cs_constraints()
{
    for(size_t ch=0; ch<dim; ch++){
        size_t dim_offset = ch*feature_height*feature_width;
        for(size_t i=0; i<feature_height; i++){
            for(size_t j=0; j<feature_width; j++){
                linear_combination<FieldT> A, B, C;
                A.add_term(output[dim_offset + i * feature_width + j], 1);
                B.add_term(var[ch], 1);
                C.add_term(input[dim_offset + i * feature_width + j], 1);
                C.add_term(avg[ch], -1);
                this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(A, B, C),"layer norm");
            }
        }
    }
};

template <typename FieldT>
void layer_norm_gadget<FieldT>::generate_r1cs_witness()
{
    // Var(x) = E(x^2) - E(x)^2
    // every loop means doing norm on one feature map
    for(size_t ch=0; ch<dim; ch++){
        size_t dim_offset = ch*feature_height*feature_width;
        // gather statics
        FieldT sum_t = FieldT::zero();
        FieldT sum_2_t = FieldT::zero();
        for(size_t i=0; i<feature_height; i++){
            for(size_t j=0; j<feature_width; j++){
                sum_t += this->pb.val(input[dim_offset + i * feature_width + j]);
                sum_2_t += this->pb.val(input[dim_offset + i * feature_width + j]) * this->pb.val(input[dim_offset + i * feature_width + j]);
            }
        }
        // compute avg and var wrt input
        this->pb.val(avg[ch]) = sum_t * FieldT(dim).inverse();
        this->pb.val(var[ch]) = sum_2_t * FieldT(dim).inverse() - this->pb.val(avg[ch]) * this->pb.val(avg[ch]);
        
        // for test, use some fixed value, st the norm result is simpler
        // this->pb.val(avg[ch]) = FieldT::zero();
        // this->pb.val(var[ch]) = FieldT::one();


        // do the norm
        for(size_t i=0; i<feature_height; i++){
            for(size_t j=0; j<feature_width; j++){
                this->pb.val(output[dim_offset + i * feature_width + j]) 
                = (this->pb.val(input[dim_offset + i * feature_width + j]) - this->pb.val(avg[ch])) * (this->pb.val(var[ch]).inverse());
            }
        }
    }

    // std::cout << "---------tensor in the larer norm----------" << std::endl;
    // std::cout <<"input: "<< input << std::endl;
    // std::cout <<"output: "<< output << std::endl;
    //     size_t dim2 = output.shape[1] * output.shape[2];
    // for (size_t i = 0; i<output.size(); i += dim2) {
    //     for (size_t j = 0; j<dim2; j++) {
    //         std::cout<<this->pb.val(output[i+j]).as_ulong() << " ";
    //     }
    //     std::cout<<std::endl;
    // }
};




/**
 * token mixer gadget with pool
 * based on metaformer
 * is essentially a pooling operation of 3*3 with stride==1 and padding, which ensure the tensor shape unchanged
*/
template <typename FieldT>
class pool_token_mixer_gadget : public gadget<FieldT>
{
private:
    size_t dim;
    size_t feature_size;
    std::vector<avg_pool_gadget<FieldT>> avg_pools;

public:
    // [c, h, w], pooling operates on (h,w), namely [i,:,:]
    pb_tensor<FieldT> input;
    pb_tensor<FieldT> output;
    pool_token_mixer_gadget(){};
    pool_token_mixer_gadget(protoboard<FieldT> &pb,
                pb_tensor<FieldT> &input,
                pb_tensor<FieldT> &output,
                const std::string &annotation_prefix = "") : 
                gadget<FieldT>(pb, annotation_prefix), input(input), output(output)
    {
        //specify the private values
        dim = input.shape[0];
        feature_size = input.shape[1] * input.shape[2];
        // mix features across token
        for (size_t i=0; i < dim; i++)
        {
            pb_tensor<FieldT> input_t;
            pb_tensor<FieldT> output_t;
            input_t.insert(input_t.end(), input.begin() + i*feature_size, input.begin() + (i+1)*feature_size);
            input_t.shape = {1, input.shape[1], input.shape[2]};

            output_t.insert(output_t.end(), output.begin() + i*feature_size, output.begin() + (i+1)*feature_size);
            output_t.shape = {1, output.shape[1], output.shape[2]};
            
            // for a avg_pool: input(input), output(output), kernel_h(kernel_h), kernel_w(kernel_w), stride(stride), padding(padding)
            avg_pools.push_back(avg_pool_gadget<FieldT>(pb, input_t, output_t, 3, 3, 1, 3/2, " avg pooling token mixer"));
            // output.insert(output.begin() + i*feature_size, output_t.begin(), output_t.end());

            // output.insert(output.begin() + i*feature_size, output_t.begin(), output_t.end());
        }   
        // output.shape = {input.shape[0], input.shape[1], input.shape[2]};
    };

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

template <typename FieldT>
void pool_token_mixer_gadget<FieldT>::generate_r1cs_constraints()
{
    for(size_t i=0; i < dim; i++){
        avg_pools[i].generate_r1cs_constraints();
    }
};

template <typename FieldT>
void pool_token_mixer_gadget<FieldT>::generate_r1cs_witness()
{
    for(size_t i=0; i < dim; i++){
        avg_pools[i].generate_r1cs_witness();
    }
    // std::cout << "---------tensor in the pool mixer----------" << std::endl;
    // std::cout <<"input: "<< input << std::endl;
    // std::cout <<"output: "<< output << std::endl;
    //     size_t dim2 = output.shape[1] * output.shape[2];
    // for (size_t i = 0; i<output.size(); i += dim2) {
    //     for (size_t j = 0; j<dim2; j++) {
    //         std::cout<<this->pb.val(output[i+j]).as_ulong() << " ";
    //     }
    //     std::cout<<std::endl;
    // }
};


/**
 * MLP : y = W2 * relu(W1^T * X + bias1) + bias2
 * where X [c,h,w] == [embedding_dim, N], W1, W2 [hidden_dim, embedding_dim]
*/
template <typename FieldT>
class mlp_gadget : public gadget<FieldT>
{
private:
    size_t feature_dim;
    size_t hidden_dim;
    size_t weight_size;
    pb_tensor<FieldT> weight_1;
    pb_tensor<FieldT> weight_2;
    pb_tensor<FieldT> bias_1;
    pb_tensor<FieldT> bias_2;
    pb_tensor<FieldT> intermediate_output;
    pb_tensor<FieldT> intermediate_output_2;
    std::vector<matrix_gadget<FieldT>> linear;
    std::shared_ptr<relu_gadget<FieldT>> relu;

public:
    pb_tensor<FieldT> input;
    pb_tensor<FieldT> output;
    mlp_gadget(){};
    mlp_gadget(protoboard<FieldT> &pb,
                pb_tensor<FieldT> &input,
                pb_tensor<FieldT> &output,
                const std::string &annotation_prefix = "") : 
                gadget<FieldT>(pb, annotation_prefix), input(input), output(output)
    {
        //specify the private values
        
        // size_t t_num_varaibles = this->pb.num_variables();

        feature_dim = input.shape[0];
        //hidden_dim should be a parameter
        // in vision transformer, it's typically set to be mlp_ratio * hidden dim where mlp_ratio is 4
        // hidden_dim = feature_dim / 2;
        hidden_dim =  feature_dim*4;
        weight_size = feature_dim * hidden_dim;
        weight_1.allocate(pb, weight_size, "FFN weights[1]");
        weight_1.shape = {1, hidden_dim, feature_dim};
        // [hidden_dim, embedding_dim]* [embedding_dim, N] = [hidden_dim, N]
        intermediate_output.allocate(pb, hidden_dim*input.shape[1]*input.shape[2], "FFN intermediate output of first mat");
        intermediate_output.shape = {1, hidden_dim, input.shape[1]*input.shape[2]};
        linear.push_back(matrix_gadget<FieldT>(this->pb, weight_1, input, intermediate_output, FMT(this->annotation_prefix, " FFN linear[1]")));
        intermediate_output_2.allocate(pb, hidden_dim*input.shape[1]*input.shape[2], "FFN intermediate output of relu");
        intermediate_output_2.shape = {1, hidden_dim, input.shape[1]*input.shape[2]};
        relu.reset(new relu_gadget<FieldT>(this->pb, intermediate_output, intermediate_output_2, 0, hidden_dim* input.shape[1]*input.shape[2], "FFN ReLU"));
        weight_2.allocate(pb, weight_size, "FFN weights[2]");
        weight_2.shape = {1, feature_dim, hidden_dim};
        linear.push_back(matrix_gadget<FieldT>(this->pb, weight_2, intermediate_output_2, output, FMT(this->annotation_prefix, " FFN linear[2]")));

        std::cout << "==>mlp kernel_1 " << weight_1 << std::endl;
        // std::cout << "==>mlp intermediate " << intermediate_output << std::endl;
        std::cout << "==>mlp kernel_2 " << weight_2 << std::endl;
    };

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

template <typename FieldT>
void mlp_gadget<FieldT>::generate_r1cs_constraints()
{
    linear[0].generate_r1cs_constraints();
    relu->generate_r1cs_constraints();
    linear[1].generate_r1cs_constraints();
};

template <typename FieldT>
void mlp_gadget<FieldT>::generate_r1cs_witness()
{
    // assign weights and bias matrix
    // kernel.shape = [1, h, w] = [1, out_feature, in_feature]
    for(size_t i=0; i<weight_1.shape[1]; i++){
        for(size_t j=0; j<weight_1.shape[2]; j++){
            if((i%weight_1.shape[2]) == j)
                this->pb.val(weight_1[i*weight_1.shape[2] + j]) = FieldT(1);
        }
    }
    // kernel.shape = [1, h, w] = [1, out_feature, in_feature]
    for(size_t i=0; i<weight_2.shape[1]; i++){
        for(size_t j=0; j<weight_2.shape[2]; j++){
            if((i%weight_2.shape[2]) == j)
                this->pb.val(weight_2[i*weight_2.shape[2] + j]) = FieldT(1);
        }
    }

    linear[0].generate_r1cs_witness();
    relu->generate_r1cs_witness();
    linear[1].generate_r1cs_witness();

    // std::cout << "---------tensor in the MLP----------" << std::endl;
    // std::cout <<"input: "<< input << std::endl;
    // std::cout <<"output: "<< output << std::endl;
    //     size_t dim2 = output.shape[1] * output.shape[2];
    // for (size_t i = 0; i<output.size(); i += dim2) {
    //     for (size_t j = 0; j<dim2; j++) {
    //         std::cout<<this->pb.val(output[i+j]).as_ulong() << " ";
    //     }
    //     std::cout<<std::endl;
    // }
};

// given [128, 28, 28], 224 tokens with dim==128
// three matrix Q,K,V for each head, and a FC to connect all z
// [128, d_k], where d_k==32
template <typename FieldT>
class attention_gadget : public gadget<FieldT>
{
private:
    // d_q==d_k
    size_t d_k;
    size_t d_v;
    size_t in_dim;
    size_t in_num;
    // x * weights = queries
    // softmax((q_s * k_s) / d^(1/2)) = similarity
    // similarity * v_s = outs
    // concate(outs) * weights = output
    std::vector<pb_tensor<FieldT>> q_weights;
    std::vector<pb_tensor<FieldT>> k_weights;
    std::vector<pb_tensor<FieldT>> v_weights;
    std::vector<pb_tensor<FieldT>> q_s;
    std::vector<pb_tensor<FieldT>> k_s;
    std::vector<pb_tensor<FieldT>> v_s;
    std::vector<matrix_gadget<FieldT>> q_linear;
    std::vector<matrix_gadget<FieldT>> k_linear;
    std::vector<matrix_gadget<FieldT>> v_linear;
    std::vector<matrix_gadget<FieldT>> att_linear_1;
    std::vector<matrix_gadget<FieldT>> att_linear_2;

    std::vector<pb_tensor<FieldT>> att_map;
    std::vector<relu_gadget<FieldT>> norm;
    std::vector<pb_tensor<FieldT>> norm_att_map;
    // output of every head
    pb_tensor<FieldT> z_s;
    // zs = [z0#z1#z2...], output = z_s*z_weights
    pb_tensor<FieldT> z_weights;
    std::shared_ptr<matrix_gadget<FieldT>> z_fc;

public:
    pb_tensor<FieldT> input;
    pb_tensor<FieldT> output;
    size_t head_num;
    attention_gadget(){};
    attention_gadget(protoboard<FieldT> &pb,
                pb_tensor<FieldT> &input,
                pb_tensor<FieldT> &output,
                size_t head_num,
                const std::string &annotation_prefix = "") : 
                gadget<FieldT>(pb, annotation_prefix), input(input), output(output), head_num(head_num)
    {
        //specify the private values
        d_k = 32;
        d_v = 32;
        in_dim = input.shape[0];
        in_num = input.shape[1] * input.shape[2];
        z_s.shape = {1, d_v*head_num, in_num};
        z_s.allocate(pb, d_v*head_num*in_num, "output from all heads");
        z_weights.shape = {1, output.shape[0], d_v*head_num};
        z_weights.allocate(pb, d_v*head_num*output.shape[0], "last FC");
        for(size_t i=0; i<head_num; i++){
            q_weights.push_back(pb_tensor<FieldT>(in_dim*d_k, 0));
            q_weights[i].shape={1, d_k, in_dim};
            q_weights[i].allocate(pb, in_dim*d_k, "Q weights[1]");
            k_weights.push_back(pb_tensor<FieldT>(in_dim*d_k, 0));
            k_weights[i].shape={1, d_k, in_dim};
            k_weights[i].allocate(pb, in_dim*d_k, "K weights[1]");
            v_weights.push_back(pb_tensor<FieldT>(in_dim*d_v, 0));
            v_weights[i].shape={1, d_v, in_dim};
            v_weights[i].allocate(pb, in_dim*d_v, "V weights[1]");

            q_s.push_back(pb_tensor<FieldT>(in_num*d_k, 0));
            q_s[i].shape={1, d_k, in_num};
            q_s[i].allocate(pb, in_num*d_k, "Query i [1]");
            k_s.push_back(pb_tensor<FieldT>(in_num*d_k, 0));
            k_s[i].shape={1, d_k, in_num};
            k_s[i].allocate(pb, in_num*d_k, "Key i [1]");
            v_s.push_back(pb_tensor<FieldT>(in_num*d_k, 0));
            v_s[i].shape={1, d_v, in_num};
            v_s[i].allocate(pb, in_num*d_v, "Value i [1]");

            q_linear.push_back(matrix_gadget<FieldT>(pb, q_weights[i], input, q_s[i], "compute Query"));
            k_linear.push_back(matrix_gadget<FieldT>(pb, k_weights[i], input, k_s[i], "compute Key"));
            v_linear.push_back(matrix_gadget<FieldT>(pb, v_weights[i], input, v_s[i], "compute Value"));
        

            att_map.push_back(pb_tensor<FieldT>(in_num*in_num, 0));
            att_map[i].shape = {1, in_num, in_num};
            att_map[i].allocate(pb, in_num*in_num, "attention values");
            
            norm_att_map.push_back(pb_tensor<FieldT>(in_num*in_num, 0));
            norm_att_map[i].shape = {1, in_num, in_num};
            norm_att_map[i].allocate(pb, in_num*in_num, "attention values");
            // z_i = softmax(qk/d^0.5) * v
            // [in_num, d_k] * [d_k, in_num]
            q_s[i].shape={1, in_num, d_k};
            att_linear_1.push_back(matrix_gadget<FieldT>(pb, q_s[i], k_s[i], att_map[i], "attntion value"));
            norm.push_back(relu_gadget<FieldT>(pb, att_map[i], norm_att_map[i], 0, in_num*in_num, "softmaxed attention value"));
            pb_tensor<FieldT> head_output_t;
            // [d_v, n] * [n ,n]
            head_output_t.shape = {1, d_v, in_num};
            head_output_t.insert(head_output_t.end(), z_s.begin() + i * d_v*in_num, z_s.begin() + (i+1) *d_v*in_num);
            att_linear_2.push_back(matrix_gadget<FieldT>(pb, v_s[i], norm_att_map[i], head_output_t, "compute zs"));
        }
        // z = linear (z1#z2#...)
        // std::cout << "==>z weightst: " << z_weights << std::endl;
        // std::cout << "==>z s: " << z_s << std::endl;
        z_fc.reset(new matrix_gadget<FieldT>(this->pb, z_weights, z_s, output, "fc in attention"));
        std::cout << "==>multihead attention: " << head_num << std::endl;
    };

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

template <typename FieldT>
void attention_gadget<FieldT>::generate_r1cs_constraints()
{
    for(size_t i=0; i<head_num; i++){
        q_linear[i].generate_r1cs_constraints();
        k_linear[i].generate_r1cs_constraints();
        v_linear[i].generate_r1cs_constraints();
        att_linear_1[i].generate_r1cs_constraints();
        norm[i].generate_r1cs_constraints();
        att_linear_2[i].generate_r1cs_constraints();
    }
    z_fc->generate_r1cs_constraints();
};

template <typename FieldT>
void attention_gadget<FieldT>::generate_r1cs_witness()
{
    // assign weights; [dk, in_dim]
    for(size_t i=0; i<head_num; i++){
        for(size_t weight_i=0; weight_i<d_k; weight_i++){
            for(size_t weight_j=0; weight_j<in_dim; weight_j++){
                if((weight_i%in_dim) == weight_j){
                    this->pb.val(q_weights[i][weight_i*in_dim + weight_j]) = FieldT(1);
                    this->pb.val(k_weights[i][weight_i*in_dim + weight_j]) = FieldT(1);
                    this->pb.val(v_weights[i][weight_i*in_dim + weight_j]) = FieldT(1);
                }
            }
        }
    }
    // std::cout << output.shape[0] << std::endl;
    for(size_t i=0; i<output.shape[0]; i++){
        for(size_t j=0; j<d_v*head_num; j++){
            if((i%(d_v*head_num)) == j)
                this->pb.val(z_weights[i*d_v*head_num + j]) = FieldT(1);
        }
    }
    // gen wit for gadgets
    for(size_t i=0; i<head_num; i++){
        q_linear[i].generate_r1cs_witness();
        k_linear[i].generate_r1cs_witness();
        v_linear[i].generate_r1cs_witness();
        att_linear_1[i].generate_r1cs_witness();
        norm[i].generate_r1cs_witness();
        att_linear_2[i].generate_r1cs_witness();
    }
    // std::cout << "test" << std::endl;
    z_fc->generate_r1cs_witness();
};



template <typename FieldT>
class linear_attention_gadget : public gadget<FieldT>
{
private:
    // d_q==d_k
    size_t d_k;
    size_t d_v;
    size_t in_dim;
    size_t in_num;
    // x * weights = queries
    // softmax((q_s * k_s) / d^(1/2)) = similarity
    // similarity * v_s = outs
    // concate(outs) * weights = output
    std::vector<pb_tensor<FieldT>> q_weights;
    std::vector<pb_tensor<FieldT>> k_weights;
    std::vector<pb_tensor<FieldT>> v_weights;
    std::vector<pb_tensor<FieldT>> q_s;
    std::vector<pb_tensor<FieldT>> k_s;
    std::vector<pb_tensor<FieldT>> v_s;
    std::vector<matrix_gadget<FieldT>> q_linear;
    std::vector<matrix_gadget<FieldT>> k_linear;
    std::vector<matrix_gadget<FieldT>> v_linear;
    std::vector<matrix_gadget<FieldT>> att_linear_1;
    std::vector<matrix_gadget<FieldT>> att_linear_2;
    // in efficient attention att_map is global context vector
    std::vector<pb_tensor<FieldT>> att_map;
    // output of every head
    pb_tensor<FieldT> z_s;
    // zs = [z0#z1#z2...], output = z_s*z_weights
    pb_tensor<FieldT> z_weights;
    std::shared_ptr<matrix_gadget<FieldT>> z_fc;

public:
    pb_tensor<FieldT> input;
    pb_tensor<FieldT> output;
    size_t head_num;
    linear_attention_gadget(){};
    linear_attention_gadget(protoboard<FieldT> &pb,
                pb_tensor<FieldT> &input,
                pb_tensor<FieldT> &output,
                size_t head_num,
                const std::string &annotation_prefix = "") : 
                gadget<FieldT>(pb, annotation_prefix), input(input), output(output), head_num(head_num)
    {
        //specify the private values
        d_k = 32;
        d_v = 32;
        in_dim = input.shape[0];
        in_num = input.shape[1] * input.shape[2];
        z_s.shape = {1, d_v*head_num, in_num};
        z_s.allocate(pb, d_v*head_num*in_num, "output from all heads");
        z_weights.shape = {1, output.shape[0], d_v*head_num};
        z_weights.allocate(pb, d_v*head_num*output.shape[0], "last FC");
        for(size_t i=0; i<head_num; i++){
            q_weights.push_back(pb_tensor<FieldT>(in_dim*d_k, 0));
            q_weights[i].shape={1, d_k, in_dim};
            q_weights[i].allocate(pb, in_dim*d_k, "Q weights[1]");
            k_weights.push_back(pb_tensor<FieldT>(in_dim*d_k, 0));
            k_weights[i].shape={1, d_k, in_dim};
            k_weights[i].allocate(pb, in_dim*d_k, "K weights[1]");
            v_weights.push_back(pb_tensor<FieldT>(in_dim*d_v, 0));
            v_weights[i].shape={1, d_v, in_dim};
            v_weights[i].allocate(pb, in_dim*d_v, "V weights[1]");

            q_s.push_back(pb_tensor<FieldT>(in_num*d_k, 0));
            q_s[i].shape={1, d_k, in_num};
            q_s[i].allocate(pb, in_num*d_k, "Query i [1]");
            k_s.push_back(pb_tensor<FieldT>(in_num*d_k, 0));
            k_s[i].shape={1, d_k, in_num};
            k_s[i].allocate(pb, in_num*d_k, "Key i [1]");
            v_s.push_back(pb_tensor<FieldT>(in_num*d_k, 0));
            v_s[i].shape={1, d_v, in_num};
            v_s[i].allocate(pb, in_num*d_v, "Value i [1]");

            // [d_k, n]
            q_linear.push_back(matrix_gadget<FieldT>(pb, q_weights[i], input, q_s[i], "compute Query"));
            k_linear.push_back(matrix_gadget<FieldT>(pb, k_weights[i], input, k_s[i], "compute Key"));
            v_linear.push_back(matrix_gadget<FieldT>(pb, v_weights[i], input, v_s[i], "compute Value"));
        
            // in linear attention, the computation order is different
            // q * (k*v)
            if(d_k*d_v < in_num *in_num){
                att_map.push_back(pb_tensor<FieldT>(d_k*d_v, 0));
                att_map[i].shape = {1, d_v, d_k};
                att_map[i].allocate(pb, d_v*d_k, "attention values");
                

                // z_i = softmax(qk/d^0.5) * v
                // [in_num, d_k] * [d_k, in_num]
                // v * k^T = [dv, innum] * [inmun, dk]
                k_s[i].shape={1, in_num, d_k};
                att_linear_1.push_back(matrix_gadget<FieldT>(pb, v_s[i], k_s[i], att_map[i], "attntion value"));

                pb_tensor<FieldT> head_output_t;
                // [d_v, n] * [n ,n]
                head_output_t.shape = {1, d_v, in_num};
                head_output_t.insert(head_output_t.end(), z_s.begin() + i * d_v*in_num, z_s.begin() + (i+1) *d_v*in_num);
                att_linear_2.push_back(matrix_gadget<FieldT>(pb, att_map[i], q_s[i], head_output_t, "compute zs"));
            } else{
                att_map.push_back(pb_tensor<FieldT>(in_num*in_num, 0));
                att_map[i].shape = {1, in_num, in_num};
                att_map[i].allocate(pb, in_num*in_num, "attention values");
                

                // z_i = softmax(qk/d^0.5) * v
                // [in_num, d_k] * [d_k, in_num]
                // v * k^T = [dv, innum] * [inmun, dk]
                q_s[i].shape={1, in_num, d_k};
                att_linear_1.push_back(matrix_gadget<FieldT>(pb, q_s[i], k_s[i], att_map[i], "attntion value"));

                pb_tensor<FieldT> head_output_t;
                // [d_v, n] * [n ,n]
                head_output_t.shape = {1, d_v, in_num};
                head_output_t.insert(head_output_t.end(), z_s.begin() + i * d_v*in_num, z_s.begin() + (i+1) *d_v*in_num);
                att_linear_2.push_back(matrix_gadget<FieldT>(pb, v_s[i], att_map[i], head_output_t, "compute zs"));
            }
        }
        // z = linear (z1#z2#...)
        // std::cout << "==>z weightst: " << z_weights << std::endl;
        // std::cout << "==>z s: " << z_s << std::endl;
        z_fc.reset(new matrix_gadget<FieldT>(this->pb, z_weights, z_s, output, "fc in attention"));
        std::cout << "==>multihead linear attention: " << head_num << std::endl;
    };

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

template <typename FieldT>
void linear_attention_gadget<FieldT>::generate_r1cs_constraints()
{
    for(size_t i=0; i<head_num; i++){
        q_linear[i].generate_r1cs_constraints();
        k_linear[i].generate_r1cs_constraints();
        v_linear[i].generate_r1cs_constraints();
        att_linear_1[i].generate_r1cs_constraints();
        att_linear_2[i].generate_r1cs_constraints();
    }
    z_fc->generate_r1cs_constraints();
};

template <typename FieldT>
void linear_attention_gadget<FieldT>::generate_r1cs_witness()
{
    // assign weights; [dk, in_dim]
    for(size_t i=0; i<head_num; i++){
        for(size_t weight_i=0; weight_i<d_k; weight_i++){
            for(size_t weight_j=0; weight_j<in_dim; weight_j++){
                if((weight_i%in_dim) == weight_j){
                    this->pb.val(q_weights[i][weight_i*in_dim + weight_j]) = FieldT(1);
                    this->pb.val(k_weights[i][weight_i*in_dim + weight_j]) = FieldT(1);
                    this->pb.val(v_weights[i][weight_i*in_dim + weight_j]) = FieldT(1);
                }
            }
        }
    }
    // std::cout << output.shape[0] << std::endl;
    for(size_t i=0; i<output.shape[0]; i++){
        for(size_t j=0; j<d_v*head_num; j++){
            if((i%(d_v*head_num)) == j)
                this->pb.val(z_weights[i*d_v*head_num + j]) = FieldT(1);
        }
    }
    // gen wit for gadgets
    for(size_t i=0; i<head_num; i++){
        q_linear[i].generate_r1cs_witness();
        k_linear[i].generate_r1cs_witness();
        v_linear[i].generate_r1cs_witness();
        att_linear_1[i].generate_r1cs_witness();
        att_linear_2[i].generate_r1cs_witness();
    }
    // std::cout << "test" << std::endl;
    z_fc->generate_r1cs_witness();
};



//depth conv for given [c,h,w], kernel is [c, k, k] doing 2D convfor each channel
//same padding is used, c convolutions 
// template <typename FieldT>
// class depth_conv_gadget : public gadget<FieldT>
// {
// private:
//     std::vector<conv_gadget<FieldT>> convs;
//     std::vector<pb_tensor<FieldT>> kernels;
//     size_t kernel_size;
//     size_t padding_size;
//     size_t input_dim;
//     size_t output_dim;

// public:
//     pb_tensor<FieldT> input;
//     depth_conv_gadget(){};
//     depth_conv_gadget(protoboard<FieldT> &pb,
//                 pb_tensor<FieldT> &input,
//                 const std::string &annotation_prefix = "") : 
//                 gadget<FieldT>(pb, annotation_prefix), input(input)
//     {
//         //specify the private values
//     };

//     void generate_r1cs_constraints();
//     void generate_r1cs_witness();
// };

// template <typename FieldT>
// void depth_conv_gadget<FieldT>::generate_r1cs_constraints()
// {
// };

// template <typename FieldT>
// void depth_conv_gadget<FieldT>::generate_r1cs_witness()
// {
// };


// input  [Cin, H, W]  (mat_mult with [cout, cin])
// output [Cout, H, W]
// template <typename FieldT>
// class spatial_mlp_gadget : public gadget<FieldT>
// {
// private:
//     size_t output_w;

// public:
//     pb_tensor<FieldT> input;
//     spatial_mlp_gadget(){};
//     spatial_mlp_gadget(protoboard<FieldT> &pb,
//                 pb_tensor<FieldT> &input,
//                 const std::string &annotation_prefix = "") : 
//                 gadget<FieldT>(pb, annotation_prefix), input(input)
//     {
//         //specify the private values
//     };

//     void generate_r1cs_constraints();
//     void generate_r1cs_witness();
// };

// template <typename FieldT>
// void spatial_mlp_gadget<FieldT>::generate_r1cs_constraints()
// {
// };

// template <typename FieldT>
// void spatial_mlp_gadget<FieldT>::generate_r1cs_witness()
// {
// };



template <typename FieldT>
class residual_gadget : public gadget<FieldT>
{
private:
    size_t output_w;

public:
    pb_tensor<FieldT> input;
    residual_gadget(){};
    residual_gadget(protoboard<FieldT> &pb,
                pb_tensor<FieldT> &input,
                const std::string &annotation_prefix = "") : 
                gadget<FieldT>(pb, annotation_prefix), input(input)
    {
        //specify the private values
    };

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

template <typename FieldT>
void residual_gadget<FieldT>::generate_r1cs_constraints()
{
};

template <typename FieldT>
void residual_gadget<FieldT>::generate_r1cs_witness()
{
};



// template <typename FieldT>
// class conv_gadget : public gadget<FieldT>
// {
// private:
//     size_t output_w;

// public:
//     pb_tensor<FieldT> input;
//     templ_gadget(){};
//     templ_gadget(protoboard<FieldT> &pb,
//                 pb_tensor<FieldT> &input,
//                 const std::string &annotation_prefix = "") : 
//                 gadget<FieldT>(pb, annotation_prefix), input(input), input_h(i_h), input_w(i_w), kernel(kernel), kernel_h(k_h), kernel_w(k_w), output(output)
//     {
//         //specify the private values
//     };

//     void generate_r1cs_constraints();
//     void generate_r1cs_witness();
// };

// template <typename FieldT>
// void templ_gadget<FieldT>::generate_r1cs_constraints()
// {
// };

// template <typename FieldT>
// void templ_gadget<FieldT>::generate_r1cs_witness()
// {
// };



}