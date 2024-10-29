/** @file
 *****************************************************************************

 Implementation of functions to sample R1CS examples with prescribed parameters
 (according to some distribution).

 See r1cs_examples.hpp .

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <cassert>

#include <libff/common/utils.hpp>
#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>
#include <libsnark/relations/constraint_satisfaction_problems/r1cs/r1cs.hpp>
#include <libsnark/relations/constraint_satisfaction_problems/r1cs/examples/r1cs_examples.hpp>

#include <libsnark/models/transformer/components.hpp>

namespace libsnark {

//constructing circuits from gadgets
// gadgets(pb, variables, anno)
//note that in networks, the input of current layer is output of last layer
// the assigment holds [network_output, network_input, weights, intermediate outputs, ...]
/** one forward in a transformer block
 * for a transformer, building blocks are
 * -patch embedding
 * -transformer block
 *      -normalization
 *      -token mixer (skip connect)
 *      -normalization
 *      -MLP (skip connect)
*/
template<typename FieldT>
r1cs_example<FieldT> generate_r1cs_example_from_model()
{   
    libff::enter_block("Call to generate_r1cs_example_from_model");
    // working with pb rather than cs itself
    protoboard<FieldT> pb;
    // size_t num_constraint = pb.num_constraints();
    // size_t num_varaibles = pb.num_variables();
    // std::vector<gadget<FieldT>> gadget_list;

    // primary input consists of network_input and out_put
    // pb_tensor<FieldT> model_output(1, 0);
    // model_output.allocate(pb, model_output.size(), "model output");
    // std::cout << "==>model output " << model_output << std::endl;
    
    // input [3, 224, 224]
    // a low resolution input
    size_t input_c=512, input_h=3, input_w=3;
    // size_t stage_n = 3;
    std::vector<size_t> downscaling = {4, 2, 2, 2};
    std::vector<size_t> embedding_dimension = {64, 128, 320, 512};

    // pb_tensor<FieldT> block_output(input_c*input_h*input_w, 0);
    // block_output.shape = {input_c, input_h, input_w};
    // std::cout << "==>block input " << block_output << std::endl;
    // block_output.allocate(pb, block_output.size(), " block output ");

    pb_tensor<FieldT> model_input(input_c*input_h*input_w, 0);
    model_input.shape = {input_c, input_h, input_w};
    std::cout << "==>model input " << model_input << std::endl;
    model_input.allocate(pb, model_input.size(), " model input ");

    
    libff::enter_block("NET: construct networl layers");
    // following is a forward pass in the network
    
    // --patch embedding, downscaling of [4, 2, 2, 2], embedding dim of [64, 128, 320, 512] ([96, 192, 384, 768] in medium model)
    // --patch embedding 1: [3, 224, 224] (reshape and unfold)--> [48, 56, 56] (dense or conv)--> [64, 56, 56]
    
    // libff::enter_block("NET: patch embedding gadget");
    // pb_tensor<FieldT> embedding_1(embedding_dimension[stage_n] * (input_h/downscaling[stage_n]) * (input_w/downscaling[stage_n]), 0);
    // embedding_1.allocate(pb, embedding_1.size(), "output of patch embedding");
    // embedding_1.shape = {embedding_dimension[stage_n], (input_h/downscaling[stage_n]), (input_w/downscaling[stage_n])};
    // patch_embedding_gadget<FieldT> patch_embedding_1(pb, model_input, embedding_1, embedding_dimension[stage_n], downscaling[stage_n], "patch embedding");
    // // gadget_list.push_back(patch_embedding_1);

    // patch_embedding_1.generate_r1cs_constraints();
    // // std::cout << "#constraints in layer: "<< pb.num_constraints() - num_constraint <<std::endl;
    // // num_constraint = pb.num_constraints();
    // // std::cout << "#variables in layer: "<< pb.num_variables() -  num_varaibles <<std::endl;
    // // num_varaibles = pb.num_variables();
    // libff::leave_block("NET: patch embedding gadget");
    // std::cout << "==>embedding output " << embedding_1 << std::endl;



    libff::enter_block("NET: layer norm gadget");
    //larnorm 1: [64, 56, 56] --> [64, 56, 56]
    pb_tensor<FieldT> block_1_norm_output_1(model_input.size(), 0);
    block_1_norm_output_1.shape = model_input.shape;
    block_1_norm_output_1.allocate(pb, block_1_norm_output_1.size(), "output of layer normalization");
    layer_norm_gadget<FieldT> layer_norm_1(pb, model_input, block_1_norm_output_1, "layer normalization");

    layer_norm_1.generate_r1cs_constraints();
    // std::cout << "#constraints in layer: "<< pb.num_constraints() - num_constraint <<std::endl;
    // num_constraint = pb.num_constraints();
    // std::cout << "#variables in layer: "<< pb.num_variables() -  num_varaibles <<std::endl;
    // num_varaibles = pb.num_variables();
    libff::leave_block("NET: layer norm gadget");
    std::cout << "==>norm output " << block_1_norm_output_1 << std::endl;



    libff::enter_block("NET: token mixer gadget");
    //token mixer: : [64, 56, 56] --> [64, 56, 56]
    pb_tensor<FieldT> block_1_feature_1(block_1_norm_output_1.size(), 0);
    block_1_feature_1.shape = model_input.shape;
    block_1_feature_1.allocate(pb, block_1_feature_1.size(), "output of token mixer");
    pool_token_mixer_gadget<FieldT> pool_mixer(pb, block_1_norm_output_1, block_1_feature_1, "token mixer");

    pool_mixer.generate_r1cs_constraints();
    // std::cout << "#constraints in layer: "<< pb.num_constraints() - num_constraint <<std::endl;
    // num_constraint = pb.num_constraints();
    // std::cout << "#variables in layer: "<< pb.num_variables() -  num_varaibles <<std::endl;
    // num_varaibles = pb.num_variables();
    libff::leave_block("NET: token mixer gadget");
    std::cout << "==>token mixer output " << block_1_feature_1 << std::endl;

   
   
    libff::enter_block("NET: layer norm gadget");
    //larnorm 1: : [64, 56, 56] --> [64, 56, 56]
    pb_tensor<FieldT> block_1_norm_output_2(block_1_feature_1.size(), 0);
    block_1_norm_output_2.shape = model_input.shape;
    block_1_norm_output_2.allocate(pb, block_1_norm_output_2.size(), "output of layer normalization");
    layer_norm_gadget<FieldT> layer_norm_2(pb,  block_1_feature_1, block_1_norm_output_2, "layer normalization");
    
    layer_norm_2.generate_r1cs_constraints();
    // std::cout << "#constraints in layer: "<< pb.num_constraints() - num_constraint <<std::endl;
    // num_constraint = pb.num_constraints();
    // std::cout << "#variables in layer: "<< pb.num_variables() -  num_varaibles <<std::endl;
    // num_varaibles = pb.num_variables();
    libff::leave_block("NET: layer norm gadget");
    std::cout << "==>norm output " << block_1_norm_output_2 << std::endl;

    libff::enter_block("NET: MLP gadget");
    //mlp: [64, 56, 56] --> [64, 56, 56]
    pb_tensor<FieldT> block_1_feature_2(block_1_norm_output_2.size(), 0);
    block_1_feature_2.allocate(pb, block_1_feature_2.size(), "output of MLP");
    block_1_feature_2.shape = model_input.shape;
    mlp_gadget<FieldT> mlp_1(pb, block_1_norm_output_2, block_1_feature_2, "MLP");

    mlp_1.generate_r1cs_constraints();
    // std::cout << "#constraints in layer: "<< pb.num_constraints() - num_constraint <<std::endl;
    // num_constraint = pb.num_constraints();
    // std::cout << "#variables in layer: "<< pb.num_variables() -  num_varaibles <<std::endl;
    // num_varaibles = pb.num_variables();
    libff::leave_block("NET: MLP gadget");
    std::cout << "==>mlp output " << block_1_feature_2 << std::endl;

    libff::leave_block("NET: construct networl layers");


    //gen constraints
    libff::enter_block("Generate constraints");

    libff::leave_block("Generate constraints");




    // gen witness
    // assign value for [network_output, network_input, weights, intermediate outputs, ...]
    libff::enter_block("Generate witness");
    // the prover needs to give an assignment to form a witness
    // primary part is model's output and input
    // private part is the model's weight
    // for(size_t i=0; i<model_output.size(); i++){
    //     pb.val(model_output[i]) = FieldT(1);
    // }
    // input h * input w
    size_t input_dim = model_input.shape[1] * model_input.shape[2];
    for(size_t i=0; i<model_input.size(); i += input_dim){
        for( size_t j=0; j<input_dim; j++){
            pb.val(model_input[i+j]) = FieldT((j+1)%64);
        }
    }
    // for(size_t i=0; i<16; i++){
    //     pb.val(kernel_1[i]) = FieldT(i+1);
    // }
    // witness for every gadget
    // for(size_t i=0; i<gadget_list.size(); i++){
    //     std::cout<< " in gen witness"<< std::endl;
        
    //     gadget_list[i].generate_r1cs_witness();
    //     //  std::cout<< gadget_list[i]<< std::endl;
    // }
    // patch_embedding_1.generate_r1cs_witness();
    layer_norm_1.generate_r1cs_witness();
    pool_mixer.generate_r1cs_witness();
    layer_norm_2.generate_r1cs_witness();
    mlp_1.generate_r1cs_witness();
    libff::leave_block("Generate witness");

    // for(auto g: gadget_list){
    //     delete g;
    // }






    /* sanity checks */
    // assert(cs.num_variables() == full_variable_assignment.size());
    // assert(cs.num_variables() >= cs.primary_input_size);
    // assert(cs.num_inputs() == cs.primary_input_size);
    // assert(cs.num_constraints() == num_constraints);
    // assert(cs.is_satisfied(primary_input, auxiliary_input));
    // pb.dump_variables();
    libff::leave_block("Call to generate_r1cs_example_from_model");
    // input_sizes mean the public input to this circuit
    pb.set_input_sizes(0);
    r1cs_constraint_system<FieldT> cs(pb.get_constraint_system());
    return r1cs_example<FieldT>(std::move(cs), std::move(pb.primary_input()), std::move(pb.auxiliary_input()));

}


// a template for test
template<typename FieldT>
r1cs_example<FieldT> generate_r1cs_example_from_fc(){

    libff::enter_block("Call to generate_r1cs_example_from_fc");
    // working with pb rather than cs itself
    protoboard<FieldT> pb;

    // size_t square_size = 10;
    size_t matrix_1_h=128, matrix_1_w=64;
    pb_tensor<FieldT> matrix_1(matrix_1_h*matrix_1_w, 0);
    matrix_1.shape = {1, matrix_1_h, matrix_1_w};
    matrix_1.allocate(pb, matrix_1.size(), "left matrix");
    std::cout << "matrix_1 " << matrix_1 << std::endl;

    size_t matrix_2_h=64, matrix_2_w=49;
    pb_tensor<FieldT> matrix_2(matrix_2_h*matrix_2_w, 0);
    matrix_2.shape = {1, matrix_2_h, matrix_2_w};
    matrix_2.allocate(pb, matrix_2.size(), "right matrix");
    std::cout << "matrix_2 " << matrix_2 << std::endl;

    libff::enter_block("Set up the matrix");
    pb_tensor<FieldT> matrix_output(matrix_1_h*matrix_2_w, 0);
    matrix_output.shape = {1, matrix_1_h, matrix_2_w};
    matrix_output.allocate(pb, matrix_output.size(), "matrix multiplication output");
    matrix_gadget<FieldT> matrix_mul(pb, matrix_1, matrix_2, matrix_output, "matrix multiplication");
    libff::leave_block("Set up the matrix");


    libff::enter_block("Generate constraints");
    matrix_mul.generate_r1cs_constraints();

    libff::leave_block("Generate constraints");




    libff::enter_block("Generate witness");
    size_t tensor_dim = matrix_1.shape[1] * matrix_1.shape[2];
    for(size_t i=0; i<matrix_1.size(); i += tensor_dim){
        for( size_t j=0; j<tensor_dim; j++){
            pb.val(matrix_1[i+j]) = FieldT(j+1);
        }
    }
    tensor_dim = matrix_2.shape[1] * matrix_2.shape[2];
    for(size_t i=0; i<matrix_2.size(); i += tensor_dim){
        for( size_t j=0; j<tensor_dim; j++){
            pb.val(matrix_2[i+j]) = FieldT(j+1);
        }
    }
    matrix_mul.generate_r1cs_witness();

    libff::leave_block("Generate witness");



    libff::leave_block("Call to generate_r1cs_example_from_fc");
    pb.set_input_sizes(0);
    r1cs_constraint_system<FieldT> cs(pb.get_constraint_system());
    return r1cs_example<FieldT>(std::move(cs), std::move(pb.primary_input()), std::move(pb.auxiliary_input()));

}


// a template for test
template<typename FieldT>
r1cs_example<FieldT> generate_r1cs_example_from_relu(){

    libff::enter_block("Call to generate_r1cs_example_from_relu");
    // working with pb rather than cs itself
    protoboard<FieldT> pb;
    size_t matrix_size = 4;

    pb_tensor<FieldT> input_matrix(matrix_size*matrix_size, 0);
    input_matrix.shape={1, matrix_size, matrix_size};
    input_matrix.allocate(pb, input_matrix.size(), "relu input");
    std::cout << "relu input " << input_matrix << std::endl;

    pb_tensor<FieldT> output_matrix(matrix_size*matrix_size, 0);
    output_matrix.shape={1, matrix_size, matrix_size};
    output_matrix.allocate(pb, output_matrix.size(), "relu output");
    std::cout << "relu output " << output_matrix << std::endl;


    libff::enter_block("NET: relu gadget");
    relu_gadget<FieldT> relu(pb, input_matrix, output_matrix, 0, input_matrix.size(), "relu gate");
    libff::leave_block("NET: relu gadget");


    libff::enter_block("Generate constraints");
    relu.generate_r1cs_constraints();


    libff::leave_block("Generate constraints");




    libff::enter_block("Generate witness");
    size_t tensor_dim = input_matrix.shape[2];
    for(size_t i=0; i<input_matrix.size(); i += tensor_dim){
        for( size_t j=0; j<tensor_dim; j++){
            pb.val(input_matrix[i+j]) = FieldT(j+1) * FieldT(4).inverse() ;
        }
    }


    // pb.val(input_matrix[i+j]).num_bits == 254
    for(size_t i=0; i<input_matrix.size(); i += tensor_dim){
        for( size_t j=0; j<tensor_dim; j++){
            
            std::cout << pb.val(input_matrix[i+j]).as_ulong() << " " << pb.val(input_matrix[i+j]).num_bits<<" ";
        }

        std::cout << std::endl;
    }

    relu.generate_r1cs_witness();


    libff::leave_block("Generate witness");



    libff::leave_block("Call to generate_r1cs_example_from_relu");

    r1cs_constraint_system<FieldT> cs(pb.get_constraint_system());
    return r1cs_example<FieldT>(std::move(cs), std::move(pb.primary_input()), std::move(pb.auxiliary_input()));

}



// a template for test
template<typename FieldT>
r1cs_example<FieldT> generate_r1cs_example_from_poolformer(){

    libff::enter_block("Call to generate_r1cs_example_from_poolformer");
    // working with pb rather than cs itself
    // put statement on pb
    protoboard<FieldT> pb;
    pb_tensor<FieldT> model_output(1, 0);
    model_output.allocate(pb, model_output.size(), "model output");
    std::cout << "==> model output " << model_output << std::endl;

    // input [3, 224, 224]
    // a low resolution input
    size_t input_c=1, input_h=16, input_w=16;
    std::vector<size_t> num_layers = {2, 2, 6, 2};
    // std::vector<size_t> downscaling = {4, 2, 2, 2};
    std::vector<size_t> downscaling = {2, 1, 1, 1};
    std::vector<size_t> embedding_dimension = {2, 2, 2, 2};
    // std::vector<size_t> embedding_dimension = {64, 128, 320, 512};

    pb_tensor<FieldT> model_input(input_c*input_h*input_w, 0);
    model_input.shape = {input_c, input_h, input_w};
    std::cout << "==> model input " << model_input << std::endl;
    model_input.allocate(pb, model_input.size(), "model input");


    
    libff::enter_block("NET: construct networl layers");
    // 4 stages, with #layers[2, 2, 6, 2]
    std::vector<pb_tensor<FieldT>> patch_embeddings;
    std::vector<patch_embedding_gadget<FieldT>> patch_embedding_gadgets;
    // [input, stage1_out, ...]
    std::vector<pb_tensor<FieldT>> stage_out_feature;
    stage_out_feature.push_back(model_input);
    for(size_t stage_i=0; stage_i<4; stage_i++){
        libff::enter_block( "NET: transformer stage: " + std::to_string(stage_i+1));
        // every stage is composed of [patch embedding] + num_layers_t*[transformer block]
        
        libff::enter_block( "NET: transformer patch embeding: " + std::to_string(stage_i+1));
        patch_embeddings.push_back(pb_tensor<FieldT>(embedding_dimension[stage_i] * (stage_out_feature[stage_i].shape[1]/downscaling[stage_i]) * (stage_out_feature[stage_i].shape[2]/downscaling[stage_i]), 0));
        patch_embeddings[stage_i].allocate(pb, patch_embeddings[stage_i].size(), " patch embedding output ");
        patch_embeddings[stage_i].shape = {embedding_dimension[stage_i], stage_out_feature[stage_i].shape[1]/downscaling[stage_i], stage_out_feature[stage_i].shape[2]/downscaling[stage_i]};
        patch_embedding_gadgets.push_back(patch_embedding_gadget<FieldT>(pb, stage_out_feature[stage_i], patch_embeddings[stage_i], embedding_dimension[stage_i], downscaling[stage_i], "patch embedding: "+  std::to_string(stage_i+1)));
        // patch_embedding_gadgets[stage_i].generate_r1cs_constraints();
        // patch_embedding_gadgets[stage_i].generate_r1cs_witness();
        libff::leave_block( "NET: transformer patch embeding: " + std::to_string(stage_i+1));


        std::vector<pb_tensor<FieldT>> norm_out_1;
        std::vector<layer_norm_gadget<FieldT>> norm_gadgets_1;
        std::vector<pb_tensor<FieldT>> token_mixer_out;
        std::vector<pool_token_mixer_gadget<FieldT>> pool_token_mixer_gadgets;
        std::vector<pb_tensor<FieldT>> norm_out_2;
        std::vector<layer_norm_gadget<FieldT>> norm_gadgets_2;
        // layer output is essentially the output of mlp
        std::vector<pb_tensor<FieldT>> layer_out_feature;
        std::vector<mlp_gadget<FieldT>> mlp_gadgets;
        layer_out_feature.push_back(patch_embeddings[stage_i]);
        size_t num_layers_t = num_layers[stage_i];
        for(size_t layer_i=0; layer_i<num_layers_t; layer_i++){
            libff::enter_block( "NET: transformer stage: " + std::to_string(stage_i+1) + " block: " + std::to_string(layer_i+1));
            // every layer means a transfomer block of [norm, token_mixer, norm, mlp]
            
            // norm_1
            norm_out_1.push_back(pb_tensor<FieldT>(patch_embeddings[stage_i].size(), 0));
            norm_out_1[layer_i].allocate(pb, patch_embeddings[stage_i].size(), " layer norm 1 output ");
            norm_out_1[layer_i].shape = patch_embeddings[stage_i].shape;
            norm_gadgets_1.push_back(layer_norm_gadget<FieldT>(pb, layer_out_feature[layer_i], norm_out_1[layer_i], " layer norm 1 "));
            // norm_gadgets_1[layer_i].generate_r1cs_constraints();
            // norm_gadgets_1[layer_i].generate_r1cs_witness();

            // token_mixer
            token_mixer_out.push_back(pb_tensor<FieldT>(patch_embeddings[stage_i].size(), 0));
            token_mixer_out[layer_i].allocate(pb, patch_embeddings[stage_i].size(), " token mixer output ");
            token_mixer_out[layer_i].shape = patch_embeddings[stage_i].shape;
            pool_token_mixer_gadgets.push_back(pool_token_mixer_gadget<FieldT>(pb, norm_out_1[layer_i], token_mixer_out[layer_i], " token mixer "));
            // pool_token_mixer_gadgets[layer_i].generate_r1cs_constraints();
            // pool_token_mixer_gadgets[layer_i].generate_r1cs_witness();

            // norm_2
            norm_out_2.push_back(pb_tensor<FieldT>(patch_embeddings[stage_i].size(), 0));
            norm_out_2[layer_i].allocate(pb, patch_embeddings[stage_i].size(), " layer norm 1 output ");
            norm_out_2[layer_i].shape = patch_embeddings[stage_i].shape;
            norm_gadgets_2.push_back(layer_norm_gadget<FieldT>(pb, token_mixer_out[layer_i], norm_out_2[layer_i], " layer norm 2 "));
            // norm_gadgets_2[layer_i].generate_r1cs_constraints();
            // norm_gadgets_2[layer_i].generate_r1cs_witness();

            // mlp
            if (layer_i == num_layers_t-1){
                stage_out_feature.push_back(pb_tensor<FieldT>(patch_embeddings[stage_i].size(), 0));
                stage_out_feature[stage_i+1].allocate(pb, patch_embeddings[stage_i].size(), " mlp out ");
                stage_out_feature[stage_i+1].shape = patch_embeddings[stage_i].shape;
                mlp_gadgets.push_back(mlp_gadget<FieldT>(pb, norm_out_2[layer_i], stage_out_feature[stage_i+1], " mlp "));
                // mlp_gadgets[layer_i].generate_r1cs_constraints();
                // mlp_gadgets[layer_i].generate_r1cs_witness();
            }else{
                layer_out_feature.push_back(pb_tensor<FieldT>(patch_embeddings[stage_i].size(), 0));
                layer_out_feature[layer_i+1].allocate(pb, patch_embeddings[stage_i].size(), " mlp out ");
                layer_out_feature[layer_i+1].shape = patch_embeddings[stage_i].shape;
                mlp_gadgets.push_back(mlp_gadget<FieldT>(pb, norm_out_2[layer_i], layer_out_feature[layer_i+1], " mlp "));
                // mlp_gadgets[layer_i].generate_r1cs_constraints();
                // mlp_gadgets[layer_i].generate_r1cs_witness();
            }

            libff::leave_block( "NET: transformer stage: " + std::to_string(stage_i+1) + " block: " + std::to_string(layer_i+1));
        }
        libff::leave_block( "NET: transformer stage: " + std::to_string(stage_i+1));
        std::cout << "==> stage " << stage_i << " output " << stage_out_feature[stage_i+1] << std::endl;
    }

    libff::leave_block("NET: construct networl layers");





    libff::enter_block("Generate constraints");


    pb.add_r1cs_constraint(r1cs_constraint<FieldT>(FieldT::one(), FieldT::one(), model_output[0]),"model output check");
    libff::leave_block("Generate constraints");




    libff::enter_block("Generate witness");

    for(size_t i=0; i<model_output.size(); i++){
        pb.val(model_output[i]) = FieldT::one();
    }
    // input h * input w
    size_t input_dim = model_input.shape[1] * model_input.shape[2];
    for(size_t i=0; i<model_input.size(); i += input_dim){
        for( size_t j=0; j<input_dim; j++){
            pb.val(model_input[i+j]) = FieldT((j+1)%64);
        }
    }


    libff::leave_block("Generate witness");



    libff::leave_block("Call to generate_r1cs_example_from_poolformer");
    pb.set_input_sizes(model_output.size() + model_input.size());
    r1cs_constraint_system<FieldT> cs(pb.get_constraint_system());
    return r1cs_example<FieldT>(std::move(cs), std::move(pb.primary_input()), std::move(pb.auxiliary_input()));

}


// a template for test

template<typename FieldT>
r1cs_example<FieldT> generate_r1cs_example_from_component(){

    libff::enter_block("Call to generate_r1cs_example_from_component");
    protoboard<FieldT> pb;

    // computation with with weights
    // size_t input_h=256, input_w=48, kernel_h=48, kernel_w=64;
    // pb_tensor<FieldT> unit_output(input_h*kernel_w, 0);
    // unit_output.shape = {1, input_h, kernel_w};
    // unit_output.allocate(pb, unit_output.size(), " component output ");
    // std::cout << "component output: " << unit_output << std::endl;

    // pb_tensor<FieldT> unit_input(input_h*input_w, 0);
    // unit_input.allocate(pb, unit_input.size(), " component input ");
    // unit_input.shape = {1, input_h, input_w};
    // std::cout << "component input: " << unit_input  << std::endl;
    // for(size_t i=0; i<unit_input.shape[1]; i++){
    //     for( size_t j=0; j<unit_input.shape[2]; j++){
            

    //         pb.val(unit_input[i*unit_input.shape[2]+j]) = FieldT((j+1)%64);

    //         // if((i%unit_input.shape[2]) == j){
    //         //     pb.val(unit_input[i*unit_input.shape[2] + j]) = FieldT(1);
    //         // }
    //     }
    // }

    // pb_tensor<FieldT> unit_kernel(kernel_h*kernel_w, 0);
    // unit_kernel.shape = {1, kernel_h, kernel_w};
    // unit_kernel.allocate(pb, unit_kernel.size(), " component weight ");
    // std::cout << "component weight: " << unit_kernel  << std::endl;

    // for(size_t i=0; i<unit_kernel.shape[1]; i++){
    //     for(size_t j=0; j<unit_kernel.shape[2]; j++){
    //         if((i % unit_kernel.shape[2]) == j)
    //             pb.val(unit_kernel[i*unit_kernel.shape[2] + j]) = FieldT(1);
    //         // pb.val(unit_kernel[i*unit_kernel.shape[2] + j]) = FieldT((j+1)%64);
    //     }
    // }

    // matrix_gadget<FieldT> mat_mult(pb, unit_input, unit_kernel, unit_output, " test mat_mult ");
    // mat_mult.generate_r1cs_constraints();
    // mat_mult.generate_r1cs_witness();



    // computation without weights
    size_t input_h=64, input_w=64;
    // size_t pool_kernel_size = 2;
    pb_tensor<FieldT> unit_output(input_h*input_w, 0);
    unit_output.shape = {1, input_h, input_w};
    unit_output.allocate(pb, unit_output.size(), " component output ");
    std::cout << "==>component output: " << unit_output << std::endl;

    pb_tensor<FieldT> unit_input(input_h*input_w, 0);
    unit_input.allocate(pb, unit_input.size(), " component input ");
    unit_input.shape = {1, input_h, input_w};
    std::cout << "==>component input: " << unit_input  << std::endl;
    for(size_t i=0; i<unit_input.shape[1]; i++){
        for( size_t j=0; j<unit_input.shape[2]; j++){
            pb.val(unit_input[i*unit_input.shape[2]+j]) = FieldT((j+1)%64);
        }
    }

    // relu_gadget<FieldT> relu_test(pb, unit_input, unit_output, 0, unit_input.size(), " test relu ");
    // relu_test.generate_r1cs_constraints();
    // relu_test.generate_r1cs_witness();

    // relu_approx_gadget<FieldT> relu_approx_test(pb, unit_input, unit_output, " test relu approximation ");
    // relu_approx_test.generate_r1cs_constraints();
    // relu_approx_test.generate_r1cs_witness();


    // pool_gadget<FieldT> pool_test(pb, unit_input, input_h, input_w, unit_output, pool_kernel_size, pool_kernel_size, " test 2d pool ");
    // pool_test.generate_r1cs_constraints();
    // pool_test.generate_r1cs_witness();

    sigmoid_approx_gadget<FieldT> sigmoid_approx_test(pb, unit_input, unit_output, " test relu approximation ");
    sigmoid_approx_test.generate_r1cs_constraints();
    sigmoid_approx_test.generate_r1cs_witness();


    // test 3d convolution
    // size_t input_c=3, input_h=32, input_w=32;
    // size_t output_c=32;
    // pb_tensor<FieldT> unit_input(input_c*input_h*input_w, 0);
    // unit_input.allocate(pb, unit_input.size(), " component input ");
    // unit_input.shape = {input_c, input_h, input_w};
    // std::cout << "==>component input: " << unit_input  << std::endl;
    // for(size_t i=0; i<unit_input.shape[0]; i++){
    //     for( size_t j=0; j<unit_input.shape[1] * unit_input.shape[2]; j++){
    //         pb.val(unit_input[i * unit_input.shape[1] *unit_input.shape[2]+j]) = FieldT((j+1)%64);
    //     }
    // }

    // pb_tensor<FieldT> unit_output(output_c*input_h*input_w, 0);
    // unit_output.allocate(pb, unit_output.size(), " component output ");
    // unit_output.shape = {output_c, input_h, input_w};
    // std::cout << "==>component output: " << unit_output  << std::endl;

    // // in, out, kernel_size
    // conv_3d_gadget<FieldT> conv_3d_test(pb, unit_input, unit_output, 3, " test 3d conv ");
    // conv_3d_test.generate_r1cs_constraints();
    // conv_3d_test.generate_r1cs_witness();



    libff::enter_block("Generate constraints");
    libff::leave_block("Generate constraints");




    libff::enter_block("Generate witness");
    libff::leave_block("Generate witness");



    libff::leave_block("Call to generate_r1cs_example_from_component");
    pb.set_input_sizes(unit_output.size());
    // pb.set_input_sizes(0);
    // pb.dump_variables();
    r1cs_constraint_system<FieldT> cs(pb.get_constraint_system());
    return r1cs_example<FieldT>(std::move(cs), std::move(pb.primary_input()), std::move(pb.auxiliary_input()));

}



// a template for test
template<typename FieldT>
r1cs_example<FieldT> generate_r1cs_example_from_layer(){

    libff::enter_block("Call to generate_r1cs_example_from_layer");
    protoboard<FieldT> pb;
    // public: output of layer
    // size_t input_c=320, input_h=4, input_w=4;

    // primary input consists of network_input and out_put
    // pb_tensor<FieldT> model_output(1, 0);
    // model_output.allocate(pb, model_output.size(), "model output");
    // std::cout << "==>model output " << model_output << std::endl;

    size_t input_c=512, input_h=2, input_w=2;
    // size_t kernel_c=0, kernel_h=0, kernel_w=0;
    size_t stage_n = 0;
    std::vector<size_t> downscaling = {4, 2, 2, 2};
    std::vector<size_t> embedding_dimension = {64, 128, 320, 512};
    // std::vector<size_t> embedding_dimension = {16, 16, 16, 16};
    
    pb_tensor<FieldT> layer_input(input_c*input_h*input_w, 0);
    layer_input.shape = {input_c, input_h, input_w};
    layer_input.allocate(pb, layer_input.size(), " layer input ");
    std::cout << "==>layer input " << layer_input << std::endl;
    for(size_t i=0; i<layer_input.shape[0]; i++){
        for(size_t j=0; j<layer_input.shape[1] * layer_input.shape[2]; j++){
            pb.val(layer_input[i*layer_input.shape[1] * layer_input.shape[2] + j]) = FieldT((j+1) % 64);
        }
    }



    // test for embedding layer which downsamples
    // pb_tensor<FieldT> layer_output(embedding_dimension[stage_n] * (input_h/downscaling[stage_n]) * (input_w/downscaling[stage_n]), 0);
    // layer_output.shape = {embedding_dimension[stage_n], (input_h/downscaling[stage_n]), (input_w/downscaling[stage_n])};
    // layer_output.allocate(pb, layer_output.size(), "layer output");
    // std::cout << "==>layer output " << layer_output << std::endl;

    // test for other layers, the shape of output is consistent with input
    pb_tensor<FieldT> layer_output(input_c*input_h*input_w, 0);
    layer_output.shape = {input_c, input_h, input_w};
    layer_output.allocate(pb, layer_output.size(), "layer output");
    std::cout << "==>layer output " << layer_output << std::endl;

    // pb_tensor<FieldT> layer_input(input_c*input_h*input_w, 0);
    // layer_input.shape = {input_c, input_h, input_w};
    // layer_input.allocate(pb, layer_input.size(), " layer input ");
    // std::cout << "==>layer input " << layer_input << std::endl;
    // for(size_t i=0; i<layer_input.shape[0]; i++){
    //     for(size_t j=0; j<layer_input.shape[1] * layer_input.shape[2]; j++){
    //         pb.val(layer_input[i*layer_input.shape[1] * layer_input.shape[2] + j]) = FieldT((j+1) % 64);
    //     }
    // }

    // patch_embedding_gadget<FieldT> patch_embedding_test(pb, layer_input, layer_output, embedding_dimension[stage_n], downscaling[stage_n], " patch embedding ");
    // patch_embedding_test.generate_r1cs_constraints();
    // patch_embedding_test.generate_r1cs_witness();

    layer_norm_gadget<FieldT> layer_norm_test(pb, layer_input, layer_output, " test layer normalization ");
    layer_norm_test.generate_r1cs_constraints();
    layer_norm_test.generate_r1cs_witness();

    // pool_token_mixer_gadget<FieldT> pool_mixer_test(pb, layer_input, layer_output, "token mixer");
    // pool_mixer_test.generate_r1cs_constraints();
    // pool_mixer_test.generate_r1cs_witness();

    // mlp_gadget<FieldT> mlp_test(pb, layer_input, layer_output, " test MLP ");
    // mlp_test.generate_r1cs_constraints();
    // mlp_test.generate_r1cs_witness();
    

    libff::enter_block("Generate constraints");



    libff::leave_block("Generate constraints");




    libff::enter_block("Generate witness");



    libff::leave_block("Generate witness");



    libff::leave_block("Call to generate_r1cs_example_from_layer");

    // pb.set_input_sizes(model_output.size() + layer_input.size());
    pb.set_input_sizes(0);
    r1cs_constraint_system<FieldT> cs(pb.get_constraint_system());
    return r1cs_example<FieldT>(std::move(cs), std::move(pb.primary_input()), std::move(pb.auxiliary_input()));

}






// a template for test
template<typename FieldT>
r1cs_example<FieldT> generate_r1cs_example_from_mlp(){

    libff::enter_block("Call to generate_r1cs_example_from_mlp");
    // working with pb rather than cs itself
    protoboard<FieldT> pb;
    // say the label
    pb_tensor<FieldT> model_output(1, 0);
    model_output.allocate(pb, model_output.size(), "model output");
    std::cout << "==>output " << model_output << std::endl;
    for(size_t i=0; i<model_output.size(); i++){
        pb.val(model_output[i]) = FieldT(1);
    }

    // input for mnits, 784 = 28*28
    size_t input_c=1, input_h=1, input_w=784;
    size_t hidden_dim_1=512, hidden_dim_2=512, output_dim=10;

    pb_tensor<FieldT> model_input(input_c*input_h*input_w, 0);
    model_input.shape = {input_c, input_h, input_w};
    std::cout << "==>model input " << model_input << std::endl;
    model_input.allocate(pb, model_input.size(), " model input ");
    for(size_t i=0; i<model_input.shape[1]; i++){
        for( size_t j=0; j<model_input.shape[2]; j++){
            pb.val(model_input[i*model_input.shape[2]+j]) = FieldT((j+1)%64);
        }
    }

    // three matrix [784, 512] [512, 512] [512, 10]
    // 1
    pb_tensor<FieldT> kernel_1(784*hidden_dim_1, 0);
    kernel_1.shape = {1, 784, hidden_dim_1};
    kernel_1.allocate(pb, kernel_1.size(), " kernel 1 ");
    std::cout << "==>model kernel 1 " << kernel_1 << std::endl;
    for(size_t i=0; i<kernel_1.shape[1]; i++){
        for( size_t j=0; j<kernel_1.shape[2]; j++){
            if((i%kernel_1.shape[2]) == j)
                pb.val(kernel_1[i*kernel_1.shape[2] + j]) = FieldT(1);
        }
    }
    // 2
    pb_tensor<FieldT> kernel_2(hidden_dim_1*hidden_dim_2, 0);
    kernel_2.shape = {1, hidden_dim_1, hidden_dim_2};
    kernel_2.allocate(pb, kernel_2.size(), " kernel 2 ");
    std::cout << "==>model kernel 2 " << kernel_2 << std::endl;
    for(size_t i=0; i<kernel_2.shape[1]; i++){
        for( size_t j=0; j<kernel_2.shape[2]; j++){
            if((i%kernel_2.shape[2]) == j)
                pb.val(kernel_2[i*kernel_2.shape[2] + j]) = FieldT(1);
        }
    }
    // 3
    pb_tensor<FieldT> kernel_3(hidden_dim_2*output_dim, 0);
    kernel_3.shape = {1, hidden_dim_2, output_dim};
    kernel_3.allocate(pb, kernel_3.size(), " kernel 3 ");
    std::cout << "==>model kernel 3 " << kernel_3 << std::endl;
    for(size_t i=0; i<kernel_3.shape[1]; i++){
        for( size_t j=0; j<kernel_3.shape[2]; j++){
            if((i%kernel_3.shape[2]) == j)
                pb.val(kernel_3[i*kernel_3.shape[2] + j]) = FieldT(1);
        }
    }


    libff::enter_block("NET: MLP for MNIST");

    pb_tensor<FieldT> lay_1_out(hidden_dim_1, 0);
    lay_1_out.shape = {1, 1, hidden_dim_1};
    lay_1_out.allocate(pb, hidden_dim_1, " hidden1 ");
    matrix_gadget<FieldT> fc_1(pb, model_input, kernel_1, lay_1_out, " fc1 ");
    fc_1.generate_r1cs_constraints();
    fc_1.generate_r1cs_witness();

    pb_tensor<FieldT> lay_1_out_act(hidden_dim_1, 0);
    lay_1_out_act.shape = {1, 1, hidden_dim_1};
    lay_1_out_act.allocate(pb, hidden_dim_1, " hidden1 ");
    // relu_approx_gadget<FieldT> relu_1(pb, lay_1_out, lay_1_out_act, " relu1 ");
    relu_gadget<FieldT> relu_1(pb, lay_1_out, lay_1_out_act,0, lay_1_out.size(), " relu1 ");
    relu_1.generate_r1cs_constraints();
    relu_1.generate_r1cs_witness();

    pb_tensor<FieldT> lay_2_out(hidden_dim_2, 0);
    lay_2_out.shape = {1, 1, hidden_dim_2};
    lay_2_out.allocate(pb, hidden_dim_2, " hidden1 ");
    matrix_gadget<FieldT> fc_2(pb, lay_1_out_act, kernel_2, lay_2_out, " fc2 ");
    fc_2.generate_r1cs_constraints();
    fc_2.generate_r1cs_witness();

    pb_tensor<FieldT> lay_2_out_act(hidden_dim_2, 0);
    lay_2_out_act.shape = {1, 1, hidden_dim_2};
    lay_2_out_act.allocate(pb, hidden_dim_2, " hidden2 ");
    // relu_approx_gadget<FieldT> relu_2(pb, lay_2_out, lay_2_out_act, " relu2 ");
    relu_gadget<FieldT> relu_2(pb, lay_2_out, lay_2_out_act, 0, lay_2_out.size(), " relu2 ");
    relu_2.generate_r1cs_constraints();
    relu_2.generate_r1cs_witness();

    pb_tensor<FieldT> lay_3_out(output_dim, 0);
    lay_3_out.shape = {1, 1, output_dim};
    lay_3_out.allocate(pb, output_dim, " hidden3 ");
    matrix_gadget<FieldT> fc_3(pb, lay_2_out_act, kernel_3, lay_3_out, " fc3 ");
    fc_3.generate_r1cs_constraints();
    fc_3.generate_r1cs_witness();
    
   
    libff::leave_block("NET: MLP for MNIST");


    libff::enter_block("Generate constraints");



    libff::leave_block("Generate constraints");




    libff::enter_block("Generate witness");



    libff::leave_block("Generate witness");



    libff::leave_block("Call to generate_r1cs_example_from_mlp");
    pb.set_input_sizes(model_output.size() + model_input.size());
    r1cs_constraint_system<FieldT> cs(pb.get_constraint_system());
    return r1cs_example<FieldT>(std::move(cs), std::move(pb.primary_input()), std::move(pb.auxiliary_input()));

}




















// a template for test
template<typename FieldT>
r1cs_example<FieldT> generate_r1cs_example_from_cnn(){

    libff::enter_block("Call to generate_r1cs_example_from_cnn");
    // working with pb rather than cs itself
    protoboard<FieldT> pb;
    // say the label
    pb_tensor<FieldT> model_output(1, 0);
    model_output.allocate(pb, model_output.size(), "model output");
    std::cout << "==>output " << model_output << std::endl;
    for(size_t i=0; i<model_output.size(); i++){
        pb.val(model_output[i]) = FieldT(1);
    }

    // input for cifar, [3, 32, 32]
    size_t input_c=3, input_h=32, input_w=32;
    size_t hidden_dim_1 = 32, hidden_dim_2=64;
    size_t fc_dim_1=128;
    pb_tensor<FieldT> model_input(input_c*input_h*input_w, 0);
    model_input.shape = {input_c, input_h, input_w};
    std::cout << "==>model input " << model_input << std::endl;
    model_input.allocate(pb, model_input.size(), " model input ");
    for(size_t i=0; i<model_input.shape[0]; i++){
        for( size_t j=0; j<model_input.shape[1] * model_input.shape[2]; j++){
            pb.val(model_input[i*model_input.shape[1] * model_input.shape[2]+j]) = FieldT((j+1)%64);
        }
    }


    // define fc
    pb_tensor<FieldT> kernel_1(4096*fc_dim_1, 0);
    kernel_1.shape = {1, 4096, fc_dim_1};
    kernel_1.allocate(pb, kernel_1.size(), " kernel 1 ");
    std::cout << "==>model kernel 1 " << kernel_1 << std::endl;
    for(size_t i=0; i<kernel_1.shape[1]; i++){
        for( size_t j=0; j<kernel_1.shape[2]; j++){
            if((i%kernel_1.shape[2]) == j)
                pb.val(kernel_1[i*kernel_1.shape[2] + j]) = FieldT(1);
        }
    }
    // 2
    pb_tensor<FieldT> kernel_2(fc_dim_1*10, 0);
    kernel_2.shape = {1, fc_dim_1, 10};
    kernel_2.allocate(pb, kernel_2.size(), " kernel 2 ");
    std::cout << "==>model kernel 2 " << kernel_2 << std::endl;
    for(size_t i=0; i<kernel_2.shape[1]; i++){
        for( size_t j=0; j<kernel_2.shape[2]; j++){
            if((i%kernel_2.shape[2]) == j)
                pb.val(kernel_2[i*kernel_2.shape[2] + j]) = FieldT(1);
        }
    }



    /**
     * [3, 32, 32]
     * conv: [32, 32, 32], [32, 32, 32]
     * pool: [32, 16, 16]
     * conv: [64, 16, 16], [64, 16, 16]
     * pool: [64, 8, 8]
     * fc: [1024, 512]
     * fc: [512, 10]
    */
    libff::enter_block("NET: CNN for Cifar");


    // conv + relu
    pb_tensor<FieldT> lay_1_out(hidden_dim_1*input_h*input_w, 0);
    lay_1_out.shape = {hidden_dim_1, input_h, input_w};
    lay_1_out.allocate(pb, lay_1_out.size(), " out of conv1 ");
    conv_3d_gadget<FieldT> conv_1(pb, model_input, lay_1_out, 3, " conv-1 ");
    conv_1.generate_r1cs_constraints();
    conv_1.generate_r1cs_witness();
    std::cout << "==>conv 1 output " << lay_1_out << std::endl;

    pb_tensor<FieldT> lay_1_out_act(lay_1_out.size(), 0);
    lay_1_out_act.shape = lay_1_out.shape;
    lay_1_out_act.allocate(pb, lay_1_out.size(), " hidden1 act ");
    relu_approx_gadget<FieldT> relu_1(pb, lay_1_out, lay_1_out_act, " relu1 ");
    relu_1.generate_r1cs_constraints();
    relu_1.generate_r1cs_witness();
    std::cout << "==>relu 1 output " << lay_1_out_act << std::endl;

    // conv + relu
    pb_tensor<FieldT> lay_2_out(lay_1_out.size(), 0);
    lay_2_out.shape = lay_1_out.shape;
    lay_2_out.allocate(pb, lay_2_out.size(), " out of conv2 ");
    conv_3d_gadget<FieldT> conv_2(pb, lay_1_out_act, lay_2_out, 3, " conv-2 ");
    conv_2.generate_r1cs_constraints();
    conv_2.generate_r1cs_witness();
    std::cout << "==>conv 2 output " << lay_2_out << std::endl;

    pb_tensor<FieldT> lay_2_out_act(lay_2_out.size(), 0);
    lay_2_out_act.shape = lay_2_out.shape;
    lay_2_out_act.allocate(pb, lay_2_out_act.size(), " hidden2 act ");
    relu_approx_gadget<FieldT> relu_2(pb, lay_2_out, lay_2_out_act, " relu2 ");
    relu_2.generate_r1cs_constraints();
    relu_2.generate_r1cs_witness();
    std::cout << "==>relu 2 output " << lay_2_out_act << std::endl;

    // pool
    pb_tensor<FieldT> lay_2_out_pool(lay_2_out.size()/4, 0);
    lay_2_out_pool.shape = {lay_2_out.shape[0], lay_2_out.shape[1]/2, lay_2_out.shape[2]/2};
    lay_2_out_pool.allocate(pb, lay_2_out_pool.size(), " hidden2 pool ");
    avg_pool_gadget<FieldT> pool_1(pb, lay_2_out_act, lay_2_out_pool, 2, 2, 2, 0, " avgpool ");
    pool_1.generate_r1cs_constraints();
    pool_1.generate_r1cs_witness();
    std::cout << "==>pool output " << lay_2_out_pool << std::endl;


    //stage 2


    // conv + relu
    pb_tensor<FieldT> lay_3_out(hidden_dim_2*lay_2_out_pool.shape[1]*lay_2_out_pool.shape[2], 0);
    lay_3_out.shape = {hidden_dim_2, lay_2_out_pool.shape[1], lay_2_out_pool.shape[2]};
    lay_3_out.allocate(pb, lay_3_out.size(), " out of conv3 ");
    conv_3d_gadget<FieldT> conv_3(pb, lay_2_out_pool, lay_3_out, 3, " conv-3 ");
    conv_3.generate_r1cs_constraints();
    conv_3.generate_r1cs_witness();
    std::cout << "==>conv 3 output " << lay_3_out << std::endl;

    pb_tensor<FieldT> lay_3_out_act(lay_3_out.size(), 0);
    lay_3_out_act.shape = lay_3_out.shape;
    lay_3_out_act.allocate(pb, lay_3_out.size(), " hidden3 act ");
    relu_approx_gadget<FieldT> relu_3(pb, lay_3_out, lay_3_out_act, " relu3 ");
    relu_3.generate_r1cs_constraints();
    relu_3.generate_r1cs_witness();
    std::cout << "==>relu 3 output " << lay_3_out_act << std::endl;

    // conv + relu
    pb_tensor<FieldT> lay_4_out(lay_3_out.size(), 0);
    lay_4_out.shape = lay_3_out.shape;
    lay_4_out.allocate(pb, lay_4_out.size(), " out of conv4 ");
    conv_3d_gadget<FieldT> conv_4(pb, lay_3_out_act, lay_4_out, 3, " conv-4 ");
    conv_4.generate_r1cs_constraints();
    conv_4.generate_r1cs_witness();
    std::cout << "==>conv 4 output " << lay_4_out << std::endl;

    pb_tensor<FieldT> lay_4_out_act(lay_4_out.size(), 0);
    lay_4_out_act.shape = lay_4_out.shape;
    lay_4_out_act.allocate(pb, lay_4_out_act.size(), " hidden2 act ");
    relu_approx_gadget<FieldT> relu_4(pb, lay_4_out, lay_4_out_act, " relu4 ");
    relu_4.generate_r1cs_constraints();
    relu_4.generate_r1cs_witness();
    std::cout << "==>relu 4 output " << lay_4_out_act << std::endl;

    // pool
    pb_tensor<FieldT> lay_4_out_pool(lay_4_out.size()/4, 0);
    lay_4_out_pool.shape = {lay_4_out.shape[0], lay_4_out.shape[1]/2, lay_4_out.shape[2]/2};
    lay_4_out_pool.allocate(pb, lay_4_out_pool.size(), " hidden2 pool ");
    avg_pool_gadget<FieldT> pool_2(pb, lay_4_out_act, lay_4_out_pool, 2, 2, 2, 0, " avgpool ");
    pool_2.generate_r1cs_constraints();
    pool_2.generate_r1cs_witness();
    std::cout << "==>pool output " << lay_4_out_pool << std::endl;


    // fc1
    lay_4_out_pool.shape = {1, 1, lay_4_out_pool.size()};
    std::cout << "==>pool output " << lay_4_out_pool << std::endl;
    pb_tensor<FieldT> fc_1_out(fc_dim_1, 0);
    fc_1_out.shape = {1, 1, fc_dim_1};
    fc_1_out.allocate(pb, fc_1_out.size(), " hidden1 ");
    matrix_gadget<FieldT> fc_1(pb, lay_4_out_pool, kernel_1, fc_1_out, " fc1 ");
    fc_1.generate_r1cs_constraints();
    fc_1.generate_r1cs_witness();
    std::cout << "==>fc 1 output " << fc_1_out << std::endl;


    pb_tensor<FieldT> fc_1_out_act(fc_1_out.size(), 0);
    fc_1_out_act.shape = fc_1_out.shape;
    fc_1_out_act.allocate(pb, fc_1_out_act.size(), " hidden2 act ");
    relu_approx_gadget<FieldT> relu_5(pb, fc_1_out, fc_1_out_act, " relu5 ");
    relu_5.generate_r1cs_constraints();
    relu_5.generate_r1cs_witness();
    std::cout << "==>relu 5 output " << fc_1_out_act << std::endl;


    // fc2
    pb_tensor<FieldT> fc_2_out(10, 0);
    fc_2_out.shape = {1, 1, 10};
    fc_2_out.allocate(pb, fc_2_out.size(), " hidden1 ");
    matrix_gadget<FieldT> fc_2(pb, fc_1_out_act, kernel_2, fc_2_out, " fc2 ");
    fc_2.generate_r1cs_constraints();
    fc_2.generate_r1cs_witness();
    std::cout << "==>fc 2 output " << fc_2_out << std::endl;
    
   
    libff::leave_block("NET: CNN for Cifar");


    libff::enter_block("Generate constraints");



    libff::leave_block("Generate constraints");




    libff::enter_block("Generate witness");



    libff::leave_block("Generate witness");



    libff::leave_block("Call to generate_r1cs_example_from_cnn");
    pb.set_input_sizes(model_output.size() + model_input.size());
    r1cs_constraint_system<FieldT> cs(pb.get_constraint_system());
    return r1cs_example<FieldT>(std::move(cs), std::move(pb.primary_input()), std::move(pb.auxiliary_input()));

}







// // a template for test
// template<typename FieldT>
// r1cs_example<FieldT> generate_r1cs_example_from_fc(){

//     libff::enter_block("Call to generate_r1cs_example");
//     // working with pb rather than cs itself
//        protoboard<FieldT> pb;
    //     pb_tensor<FieldT> model_output(1, 0);
    // model_output.allocate(pb, model_output.size(), "model output");
    // std::cout << "output " << model_output << " of size " << model_output.size() << std::endl;

//     libff::enter_block("NET: patch embedding gadget");
   
//     libff::leave_block("NET: patch embedding gadget");


//     libff::enter_block("Generate constraints");



//     libff::leave_block("Generate constraints");




//     libff::enter_block("Generate witness");



//     libff::leave_block("Generate witness");



//     libff::leave_block("Call to generate_r1cs_example_from_model");
//     pb.set_input_sizes(model_output.size() + model_input.size());
//     r1cs_constraint_system<FieldT> cs(pb.get_constraint_system());
//     return r1cs_example<FieldT>(std::move(cs), std::move(pb.primary_input()), std::move(pb.auxiliary_input()));

// }


}//libsnark