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

#ifndef R1CS_EXAMPLES_TCC_
#define R1CS_EXAMPLES_TCC_

#include <cassert>

#include <libff/common/utils.hpp>

namespace libsnark {

template<typename FieldT>
r1cs_example<FieldT> generate_r1cs_example_with_field_input(const size_t num_constraints,
                                                            const size_t num_inputs)
{
    libff::enter_block("Call to generate_r1cs_example_with_field_input");

    assert(num_inputs <= num_constraints + 2);

    r1cs_constraint_system<FieldT> cs;
    cs.primary_input_size = num_inputs;
    /***
     * in this example, the computation is determined by the number of constraints
     * where as if certain computation wants verifying, this dependent relation should be reversed
     * the number of aux depends on how many constaints are there.
     * and in this example a and b are defined, if more assignments are deemed as input, witness (aux) number goes down
    */
    cs.auxiliary_input_size = 2 + num_constraints - num_inputs; // TODO: explain this

    r1cs_variable_assignment<FieldT> full_variable_assignment;
    // FieldT a = FieldT::random_element();
    // FieldT b = FieldT::random_element();
    FieldT a = FieldT(3);
    FieldT b = FieldT(2);
    full_variable_assignment.push_back(a);
    full_variable_assignment.push_back(b);

    for (size_t i = 0; i < num_constraints-1; ++i)
    {
        linear_combination<FieldT> A, B, C;

        if (i % 2)
        {
            // a * b = c
            A.add_term(i+1, 1);
            B.add_term(i+2, 1);
            C.add_term(i+3, 1);
            FieldT tmp = a*b;
            full_variable_assignment.push_back(tmp);
            a = b; b = tmp;
        }
        else
        {
            // a + b = c
            B.add_term(0, 1);
            A.add_term(i+1, 1);
            A.add_term(i+2, 1);
            C.add_term(i+3, 1);
            FieldT tmp = a+b;
            full_variable_assignment.push_back(tmp);
            a = b; b = tmp;
        }

        cs.add_constraint(r1cs_constraint<FieldT>(A, B, C));
    }

    linear_combination<FieldT> A, B, C;
    FieldT fin = FieldT::zero();
    for (size_t i = 1; i < cs.num_variables(); ++i)
    {
        A.add_term(i, 1);
        B.add_term(i, 1);
        fin = fin + full_variable_assignment[i-1];
    }
    C.add_term(cs.num_variables(), 1);
    cs.add_constraint(r1cs_constraint<FieldT>(A, B, C));
    full_variable_assignment.push_back(fin.squared());

    /* split variable assignment */
    r1cs_primary_input<FieldT> primary_input(full_variable_assignment.begin(), full_variable_assignment.begin() + num_inputs);
    r1cs_primary_input<FieldT> auxiliary_input(full_variable_assignment.begin() + num_inputs, full_variable_assignment.end());


    /* sanity checks */
    assert(cs.num_variables() == full_variable_assignment.size());
    assert(cs.num_variables() >= num_inputs);
    assert(cs.num_inputs() == num_inputs);
    assert(cs.num_constraints() == num_constraints);
    assert(cs.is_satisfied(primary_input, auxiliary_input));

    libff::leave_block("Call to generate_r1cs_example_with_field_input");

    return r1cs_example<FieldT>(std::move(cs), std::move(primary_input), std::move(auxiliary_input));
}

template<typename FieldT>
r1cs_example<FieldT> generate_r1cs_example_with_binary_input(const size_t num_constraints,
                                                             const size_t num_inputs)
{
    libff::enter_block("Call to generate_r1cs_example_with_binary_input");

    assert(num_inputs >= 1);

    r1cs_constraint_system<FieldT> cs;
    cs.primary_input_size = num_inputs;
    cs.auxiliary_input_size = num_constraints; /* we will add one auxiliary variable per constraint */

    r1cs_variable_assignment<FieldT> full_variable_assignment;
    for (size_t i = 0; i < num_inputs; ++i)
    {
        full_variable_assignment.push_back(FieldT(std::rand() % 2));
    }

    size_t lastvar = num_inputs-1;
    for (size_t i = 0; i < num_constraints; ++i)
    {
        ++lastvar;
        const size_t u = (i == 0 ? std::rand() % num_inputs : std::rand() % i);
        const size_t v = (i == 0 ? std::rand() % num_inputs : std::rand() % i);

        /* chose two random bits and XOR them together:
           res = u + v - 2 * u * v
           2 * u * v = u + v - res
        */
        linear_combination<FieldT> A, B, C;
        A.add_term(u+1, 2);
        B.add_term(v+1, 1);
        if (u == v)
        {
            C.add_term(u+1, 2);
        }
        else
        {
            C.add_term(u+1, 1);
            C.add_term(v+1, 1);
        }
        C.add_term(lastvar+1, -FieldT::one());

        cs.add_constraint(r1cs_constraint<FieldT>(A, B, C));
        full_variable_assignment.push_back(full_variable_assignment[u] + full_variable_assignment[v] - full_variable_assignment[u] * full_variable_assignment[v] - full_variable_assignment[u] * full_variable_assignment[v]);
    }

    /* split variable assignment */
    r1cs_primary_input<FieldT> primary_input(full_variable_assignment.begin(), full_variable_assignment.begin() + num_inputs);
    r1cs_primary_input<FieldT> auxiliary_input(full_variable_assignment.begin() + num_inputs, full_variable_assignment.end());

    /* sanity checks */
    assert(cs.num_variables() == full_variable_assignment.size());
    assert(cs.num_variables() >= num_inputs);
    assert(cs.num_inputs() == num_inputs);
    assert(cs.num_constraints() == num_constraints);
    assert(cs.is_satisfied(primary_input, auxiliary_input));

    libff::leave_block("Call to generate_r1cs_example_with_binary_input");

    return r1cs_example<FieldT>(std::move(cs), std::move(primary_input), std::move(auxiliary_input));
}

template<typename FieldT>
r1cs_example<FieldT> generate_r1cs_example_with_conv_1(const size_t num_inputs,
                                                        const std::vector<FieldT> inputs,
                                                        const size_t num_kernels,
                                                        const std::vector<FieldT> kernels)
{
    libff::enter_block("Call to generate_r1cs_conv_1");

    assert(num_inputs >= 1);

    r1cs_constraint_system<FieldT> cs;
    cs.primary_input_size = num_inputs + num_kernels;
    cs.auxiliary_input_size = num_inputs + num_kernels - 1;
    size_t num_constraints = 2*num_inputs + 2*num_kernels - 1;
    cs.num_convol = 0;


    r1cs_variable_assignment<FieldT> full_variable_assignment;

    
    std::cout<<"kernels :";
    for(size_t i=0; i< num_kernels;i++){
        full_variable_assignment.push_back(kernels[i]);
    }

    std::cout<<"\ninputs :";
    for (size_t i = 0; i < num_inputs; ++i)
    {
        full_variable_assignment.push_back(inputs[i]);
    }

    std::cout<<"\noutputs :";
    for(size_t i=0; i<num_inputs + num_kernels-1;i++){
        FieldT y= FieldT::zero();
        for(size_t k=0; k<num_inputs;k++){
            for(size_t l=0;l<num_kernels;l++){
                if((k+l) == i)
                {
                    std::cout<<"["<<k<<"]["<<l<<"]("<<(inputs[k]*kernels[l]).as_ulong()<<")";
                    y += inputs[k] * kernels[l];
                }
            }
        }
        full_variable_assignment.push_back(y);
    }
    //std::cout<<"\n";
    
    cs.add_convol_constraint(num_inputs, num_kernels); //num_inputs + num_kernels - 1);
    std::cout<<"convol size = "<<cs.num_convol<<"\nconvol output size = "<<cs.convol_outputs_size[0]<<std::endl;
    /* split variable assignment */
    r1cs_primary_input<FieldT> primary_input(full_variable_assignment.begin(), full_variable_assignment.begin() + num_inputs + num_kernels);
    r1cs_primary_input<FieldT> auxiliary_input(full_variable_assignment.begin() + num_inputs + num_kernels, full_variable_assignment.end());

    /* sanity checks */
    assert(cs.num_variables() == full_variable_assignment.size());
    assert(cs.num_variables() >= num_inputs);
    assert(cs.num_inputs() == num_inputs);
    assert(cs.num_constraints() == num_constraints);
    assert(cs.is_satisfied(primary_input, auxiliary_input));

    libff::leave_block("Call to generate_r1cs_conv_1");

    return r1cs_example<FieldT>(std::move(cs), std::move(primary_input), std::move(auxiliary_input));
}

template<typename FieldT>
r1cs_example<FieldT> generate_r1cs_example_with_conv_1_opt(const size_t num_inputs,
                                                        const std::vector<FieldT> inputs,
                                                        const size_t num_kernels,
                                                        const std::vector<FieldT> kernels)
{

}

template<typename FieldT>
r1cs_example<FieldT> generate_r1cs_example_with_conv_2(const size_t input_h,
                                                        const size_t input_w,
                                                        const size_t kernel_h,
                                                        const size_t kernel_w)
{
    libff::enter_block("Call to generate_r1cs_example_with_conv_2");
    r1cs_constraint_system<FieldT> cs;
    //number of variable in kernel and input
    cs.primary_input_size = (input_w)*(input_h);
    //not number of y_s in output + intermideate products between a and x
    // cs.auxiliary_input_size = (kernel_w+input_w-1)*(kernel_h+input_h-1); // TODO: explain this
    size_t output_h = input_h-kernel_h+1;
    size_t output_w = input_w-kernel_w+1;
    cs.auxiliary_input_size = output_h*output_w * (1 + kernel_h*kernel_w) + (kernel_h)*(kernel_w);
    // have to verify every multiplication in convolution
    size_t num_constraints = output_h*output_w * (1 + kernel_h*kernel_w);
    r1cs_variable_assignment<FieldT> full_variable_assignment;
    
    libff::enter_block("set variables");
    size_t sumk=0;
    for(size_t i=0; i<(kernel_h)*(kernel_w); i++){
        full_variable_assignment.push_back(i+1);
        sumk+=(i+1);
    }
    size_t sumx = 0;
    for(size_t i=0; i<(input_w)*(input_h); i++){
        full_variable_assignment.push_back(i+1);
        sumx+=(i+1);
    }
    libff::leave_block("set variables");

    
    libff::enter_block("set constraints");
    // A to hold kernel, B for input and C for all the products and ys
    size_t kernel_size = kernel_h * kernel_w;
    for (size_t i = 0; i < 1; ++i)
    {
        for(size_t slide_i=0; slide_i < output_h; slide_i++)
            for (size_t slide_j=0; slide_j < output_w; slide_j++){
                size_t slide_num = (slide_i*output_w) + slide_j;
                // every slide corresponds to one (kernel * input)
                for(size_t a_i=0; a_i < kernel_h; a_i++)
                    for(size_t a_j=0; a_j < kernel_w; a_j++){
                        linear_combination<FieldT> A, B, C;
                        A.add_term((a_i*kernel_w)+a_j+1, 1);
                        size_t x_i = a_i + slide_i;
                        size_t x_j = a_j + slide_j;
                        B.add_term((kernel_h*kernel_w)+(x_i*input_w)+x_j+1, 1);
                        C.add_term((kernel_h*kernel_w) + (input_h*input_w)+ (slide_num*(kernel_size+1)) +(a_i*kernel_w)+a_j+1, 1);
                        cs.add_constraint(r1cs_constraint<FieldT>(A, B, C));
                    }
                linear_combination<FieldT> A, B, C;
                B.add_term(0,1);
                for(size_t a_i=0; a_i < kernel_h; a_i++)
                    for(size_t a_j=0; a_j < kernel_w; a_j++){
                        A.add_term((kernel_h*kernel_w) + (input_h*input_w)+ (slide_num*(kernel_size+1)) +(a_i*kernel_w)+a_j+1, 1);
                    }
                C.add_term((kernel_h*kernel_w) + (input_h*input_w)+ (slide_num*(kernel_size+1)) + kernel_size+1, 1);
                cs.add_constraint(r1cs_constraint<FieldT>(A, B, C));
            }
    }
    libff::leave_block("set constraints");

    libff::enter_block("Compute y variables");
    for(size_t slide_i=0; slide_i < output_h; slide_i++){
        for (size_t slide_j=0; slide_j < output_w; slide_j++){
            FieldT y = FieldT::zero();
            // std::cout<<"y("<<slide_i<<", "<<slide_j<<")=";
            for(size_t a_i=0; a_i < kernel_h; a_i++)
                for(size_t a_j=0; a_j < kernel_w; a_j++){
                    size_t x_i = a_i + slide_i;
                    size_t x_j = a_j + slide_j;
                    FieldT temp_y = FieldT::zero();
                    temp_y = FieldT(a_i*kernel_w+a_j+1)*FieldT(x_i*input_w+x_j+1);
                    //push all products
                    full_variable_assignment.push_back(temp_y);
                    y += temp_y;
                }
            full_variable_assignment.push_back(y);
        }
    }
    libff::leave_block("Compute y variables");

    /* split variable assignment */
    r1cs_primary_input<FieldT> primary_input(full_variable_assignment.begin(), full_variable_assignment.begin() + cs.primary_input_size);
    r1cs_primary_input<FieldT> auxiliary_input(full_variable_assignment.begin() + cs.primary_input_size, full_variable_assignment.end());

    /* sanity checks */
    assert(cs.num_variables() == full_variable_assignment.size());
    assert(cs.num_variables() >= cs.primary_input_size);
    assert(cs.num_inputs() == cs.primary_input_size);
    assert(cs.num_constraints() == num_constraints);
    assert(cs.is_satisfied(primary_input, auxiliary_input));

    libff::leave_block("Call to generate_r1cs_example_with_conv_2");

    return r1cs_example<FieldT>(std::move(cs), std::move(primary_input), std::move(auxiliary_input));
}

template<typename FieldT>
r1cs_example<FieldT> generate_r1cs_example_with_conv_2_opt(const size_t input_h,
                                                        const size_t input_w,
                                                        const size_t kernel_h,
                                                        const size_t kernel_w)
{
    libff::enter_block("Call to generate_r1cs_example_with_conv_2_opt");

    r1cs_constraint_system<FieldT> cs;
    //number of variable in kernel and input
    cs.primary_input_size = (input_w)*(input_h);
    //not number of y_s in output + intermideate y_s
    cs.auxiliary_input_size = (kernel_w+input_w-1)*(kernel_h+input_h-1)+(kernel_h)*(kernel_w); // TODO: explain this
    size_t num_constraints = 1;

    libff::enter_block("set variables");
    ///TODO need to change
    //should assign real value in the input and kernel matrix (when push back)
    r1cs_variable_assignment<FieldT> full_variable_assignment;
    size_t sumk=0;
    for(size_t i=0; i<(kernel_h)*(kernel_w); i++){
        full_variable_assignment.push_back(i+1);
        sumk+=(i+1);
    }
    size_t sumx = 0;
    for(size_t i=0; i<(input_w)*(input_h); i++){
        full_variable_assignment.push_back(i+1);
        sumx+=(i+1);
    }
    libff::leave_block("set variables");

    libff::enter_block("set constraints");
    // A to hold kernel, B for input and C for all the y_s
    // or should be only one constraint as "pruduct of sum"
    for (size_t i = 0; i < num_constraints; ++i)
    // for (size_t i = 0; i < 1; ++i)
    {
        FieldT s = FieldT(i+2);
        FieldT t = FieldT(i+3);
        linear_combination<FieldT> A, B, C;

        FieldT temp_s = FieldT::one();
        FieldT temp_t = FieldT::one();
        for(size_t a_i=kernel_h; a_i>0; --a_i)
        {
            for(size_t a_j=kernel_w; a_j>0; --a_j){
                A.add_term(((a_i-1)*kernel_w)+a_j-1+1, temp_s*temp_t);
                temp_t *= t;
            }
            temp_t = FieldT::one();
            temp_s *= s;
        }

        temp_s = FieldT::one();
        temp_t = FieldT::one();
        for(size_t x_i=0; x_i<input_h;x_i++)
        {
            for(size_t x_j=0; x_j<(input_h);x_j++){
                B.add_term((kernel_w*kernel_h)+(input_w*x_i)+x_j+1, temp_s*temp_t); 
                temp_t *= t;
            }
            temp_t = FieldT::one();
            temp_s *= s;
        }
        temp_s = FieldT::one();
        temp_t = FieldT::one();
        for(size_t y_i=0; y_i<(input_h+kernel_h-1);y_i++){
            for(size_t y_j=0; y_j<(input_w+kernel_w-1);y_j++){
                C.add_term((kernel_h*kernel_w)+(input_w*input_h)+(input_w+kernel_w-1)*y_i+y_j+1, temp_s*temp_t);
                temp_t *= t;
            }
            temp_t = FieldT::one();
            temp_s *= s;
        }

        cs.add_constraint(r1cs_constraint<FieldT>(A, B, C));
    }
    libff::leave_block("set constraints");

    libff::enter_block("Compute y variables");
    for(size_t y_i=0; y_i<(input_h+kernel_h-1); y_i++){
        for(size_t y_j=0; y_j<(input_w+kernel_w-1); y_j++){
            FieldT y = FieldT::zero();

            for(size_t x_i=0; x_i<input_h; x_i++){
                for(size_t x_j=0; x_j<input_w; x_j++){
                    for(size_t a_i=0; a_i<kernel_h; a_i++){
                        for(size_t a_j=0; a_j<kernel_w; a_j++){
                            if((x_i+ (kernel_h-a_i-1)) == y_i && (x_j+(kernel_w-a_j-1)) == y_j)
                            {
                                y += FieldT(a_i*kernel_w+a_j+1)*FieldT(x_i*input_w+x_j+1);
                            }
                        }
                    }
                }
            }
            full_variable_assignment.push_back(y);
        }
    }
    libff::leave_block("Compute y variables");

    /* split variable assignment */
    r1cs_primary_input<FieldT> primary_input(full_variable_assignment.begin(), full_variable_assignment.begin() + cs.primary_input_size);
    r1cs_primary_input<FieldT> auxiliary_input(full_variable_assignment.begin() + cs.primary_input_size, full_variable_assignment.end());

    /* sanity checks */
    assert(cs.num_variables() == full_variable_assignment.size());
    assert(cs.num_variables() >= cs.primary_input_size);
    assert(cs.num_inputs() == cs.primary_input_size);
    assert(cs.num_constraints() == num_constraints);
    assert(cs.is_satisfied(primary_input, auxiliary_input));

    libff::leave_block("Call to generate_r1cs_example_with_conv_2_opt");

    return r1cs_example<FieldT>(std::move(cs), std::move(primary_input), std::move(auxiliary_input));
}

template<typename FieldT>
r1cs_example<FieldT> generate_r1cs_example_with_matrix(const size_t matrix_1_h,
                                                        const size_t matrix_1_w,
				                                        const size_t matrix_2_h,
                                                        const size_t matrix_2_w)
{
    libff::enter_block("Call to generate_r1cs_example_with_matrix");
    r1cs_constraint_system<FieldT> cs;
    cs.primary_input_size = matrix_1_h*matrix_2_w;
    cs.auxiliary_input_size =  (matrix_1_h)*(matrix_1_w) + matrix_1_h*matrix_2_w*(matrix_1_w+1) + (matrix_2_w)*(matrix_2_h) - matrix_1_h*matrix_2_w;
    
    
    size_t output_h = matrix_1_h;
    size_t output_w = matrix_2_w;
    size_t output_size = output_h * output_w;
    size_t num_constraints = matrix_1_h*matrix_2_w*(matrix_1_w+1);

    libff::enter_block("set variables");
    r1cs_variable_assignment<FieldT> full_variable_assignment;
    for(size_t i=0; i<(matrix_1_h)*(matrix_1_w); i++){
        full_variable_assignment.push_back(i+1);
    }
    for(size_t i=0; i<(matrix_2_w)*(matrix_2_h); i++){
        full_variable_assignment.push_back(i+1);
    }
    libff::leave_block("set variables");

    libff::enter_block("Compute y variables");
    for(size_t y_i=0; y_i<output_h; y_i++){
        for(size_t y_j=0; y_j<output_w; y_j++){
            FieldT y = FieldT::zero();
            for(size_t pos=0; pos<matrix_1_w; pos++){
                FieldT temp_y = FieldT::zero();
                temp_y = FieldT(y_i*matrix_1_w + pos + 1) * FieldT(pos*matrix_2_w + y_j + 1);
                full_variable_assignment.push_back(temp_y);
                y += temp_y;
            }
            full_variable_assignment.push_back(y);
        }
    }
    libff::leave_block("Compute y variables");

    libff::enter_block("set constraints");

    for(size_t y_i=0; y_i<output_h; y_i++){
        for(size_t y_j=0; y_j<output_w; y_j++){
            for(size_t pos=0; pos<matrix_1_w; pos++){
                linear_combination<FieldT> A, B, C;
                A.add_term((y_i*matrix_1_w)+pos+1, 1);
                B.add_term((matrix_1_h*matrix_1_w)+(pos*matrix_2_w)+y_j+1,1);
                C.add_term((matrix_1_h*matrix_1_w)+(matrix_2_h*matrix_2_w)+((y_i*output_w)+y_j)*(matrix_1_w+1)+pos+1, 1);
                cs.add_constraint(r1cs_constraint<FieldT>(A, B, C));
            }
            // constraint for one sum
            linear_combination<FieldT> A, B, C;
            B.add_term(0,1);
            // A terms of sum
            for(size_t pos=0; pos < matrix_1_w; pos++){
                A.add_term((matrix_1_h*matrix_1_w)+(matrix_2_h*matrix_2_w)+((y_i*output_w)+y_j)*(matrix_1_w+1)+pos+1, 1);
            }
            C.add_term((matrix_1_h*matrix_1_w)+(matrix_2_h*matrix_2_w)+((y_i*output_w)+y_j)*(matrix_1_w+1)+matrix_1_w+1, 1);
            cs.add_constraint(r1cs_constraint<FieldT>(A, B, C));
        }
    }

    libff::leave_block("set constraints");

    /* split variable assignment */
    r1cs_primary_input<FieldT> primary_input(full_variable_assignment.begin(), full_variable_assignment.begin() + cs.primary_input_size);
    r1cs_primary_input<FieldT> auxiliary_input(full_variable_assignment.begin() + cs.primary_input_size, full_variable_assignment.end());

    /* sanity checks */
    std::cout<<"Circuit parameters::  #constraints == "<<cs.num_constraints() << ";    #variables ==  "<<cs.num_variables()  <<std::endl;
    assert(cs.num_variables() == full_variable_assignment.size());
    assert(cs.num_variables() >= cs.primary_input_size);
    assert(cs.num_inputs() == cs.primary_input_size);
    assert(cs.num_constraints() == num_constraints);
    assert(cs.is_satisfied(primary_input, auxiliary_input));

    libff::leave_block("Call to generate_r1cs_example_with_matrix");

    return r1cs_example<FieldT>(std::move(cs), std::move(primary_input), std::move(auxiliary_input));
}



// quick power to express ib 
template<typename FieldT>
FieldT qpow(FieldT a, size_t n)
{
    if (n == 0)
        return FieldT::one();
    else if (n % 2 == 1)
        return qpow(a, n - 1) * a;
    else
    {
        FieldT temp = qpow(a, n / 2);
        return temp * temp;
    }
}

template<typename FieldT>
r1cs_example<FieldT> generate_r1cs_example_with_matrix_opt(const size_t matrix_1_h,
                                                        const size_t matrix_1_w,
				                                        const size_t matrix_2_h,
                                                        const size_t matrix_2_w)
{
    libff::enter_block("Call to generate_r1cs_example_with_matrix_opt");
    r1cs_constraint_system<FieldT> cs;
    cs.primary_input_size = 0;
    cs.auxiliary_input_size = (matrix_1_h)*(matrix_1_w) + matrix_1_h*matrix_2_w*(matrix_1_w+1) + (matrix_2_w)*(matrix_2_h);
    
    size_t output_h = matrix_1_h;
    size_t output_w = matrix_2_w;

    size_t num_constraints = matrix_1_w + 1;
    // size_t num_constraints = matrix_1_w ;

    libff::enter_block("set variables");

    r1cs_variable_assignment<FieldT> full_variable_assignment;
    for(size_t i=0; i<(matrix_1_h)*(matrix_1_w); i++){
        full_variable_assignment.push_back(i+1);
    }
    for(size_t i=0; i<(matrix_2_w)*(matrix_2_h); i++){
        full_variable_assignment.push_back(i+1);
    }

    libff::enter_block("Compute y variables");
    // std::cout<< "Y-matrix" <<std::endl;
    for(size_t y_i=0; y_i<output_h; y_i++){
        for(size_t y_j=0; y_j<output_w; y_j++){
            FieldT y = FieldT::zero();
            // pos: 1~n
            for(size_t pos=0; pos<matrix_1_w; pos++){
                FieldT temp_y = FieldT::zero();
                // x[y_i, pos] * a[pos, y_j]
                temp_y = FieldT(y_i*matrix_1_w + pos + 1) * FieldT(pos*matrix_2_w + y_j + 1);
                full_variable_assignment.push_back(temp_y);
                y += temp_y;
            }
            full_variable_assignment.push_back(y);
        }

    }
    libff::leave_block("Compute y variables");

    libff::leave_block("set variables");
    
    libff::enter_block("set constraints");
    //introduce Z
    // s for z; t for z^b
    FieldT s = FieldT(2);
    FieldT t = qpow(s, matrix_2_w);

    linear_combination<FieldT> A_fin, B_fin, C_fin;

    for (size_t i = 0; i < matrix_1_w; ++i)
    {

        linear_combination<FieldT> A, B, C;

        FieldT temp_t = FieldT::one();
        for(size_t x_i = 0; x_i < matrix_1_h; x_i ++){
            A.add_term((x_i*matrix_1_w)+i+1, temp_t);
            temp_t *= t;
        }
        FieldT temp_s = FieldT::one();
        for(size_t a_j = 0; a_j < matrix_2_w; a_j ++){
            B.add_term((matrix_1_h*matrix_1_w)+(i*matrix_2_w)+a_j+1, temp_s);

            temp_s *= s;
        }

        temp_s = FieldT::one();
        temp_t = FieldT::one();
        for(size_t x_i=0; x_i<matrix_1_h; x_i++){
            for(size_t a_j=0; a_j<matrix_2_w; a_j++){
                C.add_term((matrix_1_h*matrix_1_w)+(matrix_2_h*matrix_2_w)+((x_i*output_w)+a_j)*(matrix_1_w+1)+i+1,temp_s*temp_t);
                A_fin.add_term((matrix_1_h*matrix_1_w)+(matrix_2_h*matrix_2_w)+((x_i*output_w)+a_j)*(matrix_1_w+1)+i+1,temp_s*temp_t);
                temp_s *= s;
            }
            temp_s = FieldT::one();
            temp_t *= t;
        }
        cs.add_constraint(r1cs_constraint<FieldT>(A, B, C));
        // std::cout<<std::endl;
    }

    B_fin.add_term(0,1);
    FieldT temp_s = FieldT::one();
    for(size_t y_i=0; y_i<output_h; y_i++){
        for(size_t y_j=0; y_j<output_w; y_j++){
            C_fin.add_term((matrix_1_h*matrix_1_w)+(matrix_2_h*matrix_2_w)+((y_i*output_w)+y_j)*(matrix_1_w+1)+matrix_1_w+1,temp_s);
            temp_s *= s;
        }
    }
    cs.add_constraint(r1cs_constraint<FieldT>(A_fin, B_fin, C_fin));



    libff::leave_block("set constraints");

    /* split variable assignment */
    r1cs_primary_input<FieldT> primary_input(full_variable_assignment.begin(), full_variable_assignment.begin() + cs.primary_input_size);
    r1cs_primary_input<FieldT> auxiliary_input(full_variable_assignment.begin() + cs.primary_input_size, full_variable_assignment.end());

    /* sanity checks */
    std::cout<<"Circuit parameters::  #constraints == "<<cs.num_constraints() << ";    #variables ==  "<<cs.num_variables()  <<std::endl;
    assert(cs.num_variables() == full_variable_assignment.size());
    assert(cs.num_variables() >= cs.primary_input_size);
    assert(cs.num_inputs() == cs.primary_input_size);
    assert(cs.num_constraints() == num_constraints);
    assert(cs.is_satisfied(primary_input, auxiliary_input));

    libff::leave_block("Call to generate_r1cs_example_with_matrix_opt");

    return r1cs_example<FieldT>(std::move(cs), std::move(primary_input), std::move(auxiliary_input));
}


template<typename FieldT>
r1cs_example<FieldT> generate_r1cs_example_with_matrix_opt_2(const size_t matrix_1_h,
                                                        const size_t matrix_1_w,
				                                        const size_t matrix_2_h,
                                                        const size_t matrix_2_w)
{
    libff::enter_block("Call to generate_r1cs_example_with_matrix_opt_2");
    // matrix_1_w == matrix_2_h    i.e. (a,n)*(n,b)
    r1cs_constraint_system<FieldT> cs;

    size_t output_h = matrix_1_h;
    size_t output_w = matrix_2_w;
    size_t x_block = matrix_1_w + matrix_2_h * matrix_2_w - 1;
    size_t y_num = matrix_1_h * x_block;

    cs.primary_input_size = (matrix_1_h)*(matrix_1_w);
    cs.auxiliary_input_size = y_num + (matrix_2_w)*(matrix_2_h);
    size_t num_constraints = 1;
    
    libff::enter_block("set variables");
    //TODO: need to change
    //should assign real value in the input and kernel matrix (when push back)
    r1cs_variable_assignment<FieldT> full_variable_assignment;

    for(size_t i=0; i<(matrix_1_h)*(matrix_1_w); i++){
        full_variable_assignment.push_back(i+1);
    }
    for(size_t i=0; i<(matrix_2_w)*(matrix_2_h); i++){
        full_variable_assignment.push_back(i+1);
    }
    libff::leave_block("set variables");

    libff::enter_block("Compute y variables");
    std::cout<< y_num <<std::endl;
    for(size_t x_i=0; x_i<matrix_1_h; x_i++){
        for(size_t y_i=(x_i)*x_block; y_i<(x_i+1)*x_block; y_i++){
            FieldT y = FieldT::zero();
            for(size_t x_j=0; x_j<matrix_1_w; x_j++)
                for(size_t a_i=0; a_i<matrix_2_h; a_i++)
                    for(size_t a_j=0; a_j<matrix_2_w; a_j++){
                        size_t x_index = x_i*x_block + (matrix_1_w - x_j - 1);
                        size_t a_index = a_j*matrix_2_h + a_i;
                        if(x_index + a_index == y_i){
                            FieldT temp_y = FieldT(x_i*matrix_1_w + x_j + 1) * FieldT(a_i*matrix_2_w + a_j + 1);
                            y += temp_y;
                        }
                    }
            full_variable_assignment.push_back(y);
        }
    }

    libff::leave_block("Compute y variables");

    libff::enter_block("set constraints");
    linear_combination<FieldT> A, B, C;
    FieldT s = FieldT(2);
    FieldT t = qpow(s, x_block);
    FieldT temp_t = FieldT::one();
    FieldT temp_s = FieldT::one();
    for(size_t x_i=0; x_i<matrix_1_h; x_i++){
        for(size_t x_j=matrix_1_w; x_j>0; --x_j){
            A.add_term((x_i*matrix_1_w)+x_j-1+1, temp_s*temp_t);
            temp_s *= s;
        }
        temp_s = FieldT::one();
        temp_t *= t;
    }

    temp_s = FieldT::one();
    for(size_t a_j=0; a_j<matrix_2_w; a_j++){
        for(size_t a_i=0; a_i<matrix_2_h; a_i++){
            B.add_term((matrix_1_h*matrix_1_w)+(a_i*matrix_2_w)+a_j+1, temp_s);
            temp_s *= s;
        }
    }

    temp_s = FieldT::one();
    for(size_t y_i=0; y_i<y_num; y_i++){
        C.add_term((matrix_1_h*matrix_1_w)+(matrix_2_h*matrix_2_w)+y_i+1, temp_s);
        temp_s *= s;
    }
    cs.add_constraint(r1cs_constraint<FieldT>(A, B, C));
    libff::leave_block("set constraints");

    /* split variable assignment */
    r1cs_primary_input<FieldT> primary_input(full_variable_assignment.begin(), full_variable_assignment.begin() + cs.primary_input_size);
    r1cs_primary_input<FieldT> auxiliary_input(full_variable_assignment.begin() + cs.primary_input_size, full_variable_assignment.end());

    /* sanity checks */
    std::cout<<"cs.num_constraint = "<<cs.num_constraints() <<"      defined: "<< num_constraints <<std::endl;
    std::cout<<"cs.num_variables = "<<cs.num_variables() <<"      defined: "<< full_variable_assignment.size() <<std::endl;
    std::cout<<"cs.num_inputs = "<<cs.num_inputs() <<"      defined: "<< cs.primary_input_size <<std::endl;
    assert(cs.num_variables() == full_variable_assignment.size());
    assert(cs.num_variables() >= cs.primary_input_size);
    assert(cs.num_inputs() == cs.primary_input_size);
    assert(cs.num_constraints() == num_constraints);
    assert(cs.is_satisfied(primary_input, auxiliary_input));

    libff::leave_block("Call to generate_r1cs_example_with_matrix_opt_2");

    return r1cs_example<FieldT>(std::move(cs), std::move(primary_input), std::move(auxiliary_input));

}


template<typename FieldT>
r1cs_example<FieldT> generate_r1cs_example_with_matrix_opt_3(const size_t matrix_1_h,
                                                        const size_t matrix_1_w,
				                                        const size_t matrix_2_h,
                                                        const size_t matrix_2_w)
{
    libff::enter_block("Call to generate_r1cs_example_with_matrix_opt_3");
    r1cs_constraint_system<FieldT> cs;
    cs.primary_input_size = (matrix_1_h)*(matrix_1_w)+(matrix_2_w)*(matrix_2_h);
    size_t output_h = matrix_1_h;
    size_t output_w = matrix_2_w;
    cs.auxiliary_input_size = matrix_1_h*matrix_2_w*(matrix_1_w+1);
    size_t num_constraints = 1;


    libff::enter_block("set variables");
    r1cs_variable_assignment<FieldT> full_variable_assignment;
    for(size_t i=0; i<(matrix_1_h)*(matrix_1_w); i++){
        full_variable_assignment.push_back(i+1);
    }
    for(size_t i=0; i<(matrix_2_w)*(matrix_2_h); i++){
        full_variable_assignment.push_back(i+1);
    }

    libff::enter_block("Compute y variables");
    for(size_t y_i=0; y_i<output_h; y_i++){
        for(size_t y_j=0; y_j<output_w; y_j++){
            FieldT y = FieldT::zero();
            for(size_t pos=0; pos<matrix_1_w; pos++){
                FieldT temp_y = FieldT::zero();
                temp_y = FieldT(y_i*matrix_1_w + pos + 1) * FieldT(pos*matrix_2_w + y_j + 1);
                full_variable_assignment.push_back(temp_y);
                y += temp_y;
            }
            full_variable_assignment.push_back(y);
        }
    }
    libff::leave_block("Compute y variables");

    libff::leave_block("set variables");


    libff::enter_block("set constraints");
    FieldT s = FieldT(2);
    FieldT t = qpow(s, matrix_2_w);

    linear_combination<FieldT> A_fin, B_fin, C_fin;

    FieldT temp_s = FieldT::one();
    FieldT temp_t = FieldT::one();
    for (size_t i = 0; i < matrix_1_w; ++i){
        temp_s = FieldT::one();
        temp_t = FieldT::one();
        for(size_t x_i=0; x_i<matrix_1_h; x_i++){
            for(size_t a_j=0; a_j<matrix_2_w; a_j++){
                    A_fin.add_term((matrix_1_h*matrix_1_w)+(matrix_2_h*matrix_2_w)+((x_i*output_w)+a_j)*(matrix_1_w+1)+i+1,temp_s*temp_t);
                    temp_s *= s;
            }
            temp_s = FieldT::one();
            temp_t *= t;
        }
    }

    B_fin.add_term(0,1);
    temp_s = FieldT::one();
    for(size_t y_i=0; y_i<output_h; y_i++){
        for(size_t y_j=0; y_j<output_w; y_j++){
            C_fin.add_term((matrix_1_h*matrix_1_w)+(matrix_2_h*matrix_2_w)+((y_i*output_w)+y_j)*(matrix_1_w+1)+matrix_1_w+1,temp_s);
            temp_s *= s;
        }
    }
    cs.add_constraint(r1cs_constraint<FieldT>(A_fin, B_fin, C_fin));



    libff::leave_block("set constraints");

    /* split variable assignment */
    r1cs_primary_input<FieldT> primary_input(full_variable_assignment.begin(), full_variable_assignment.begin() + cs.primary_input_size);
    r1cs_primary_input<FieldT> auxiliary_input(full_variable_assignment.begin() + cs.primary_input_size, full_variable_assignment.end());

    std::cout<<"cs.num_cs = "<<cs.num_constraints() <<"      defined"<< num_constraints <<std::endl;
    std::cout<<"cs.num_variables = "<<cs.num_variables() <<"      defined: "<< full_variable_assignment.size() <<std::endl;
    std::cout<<"cs.num_inputs = "<<cs.num_inputs() <<"      defined"<< cs.primary_input_size <<std::endl;
    assert(cs.num_variables() == full_variable_assignment.size());
    assert(cs.num_variables() >= cs.primary_input_size);
    assert(cs.num_inputs() == cs.primary_input_size);
    assert(cs.num_constraints() == num_constraints);
    assert(cs.is_satisfied(primary_input, auxiliary_input));

    libff::leave_block("Call to generate_r1cs_example_with_matrix_opt_3");

    return r1cs_example<FieldT>(std::move(cs), std::move(primary_input), std::move(auxiliary_input));
}


/**
 * optimize with prefix sum for A query
*/
template<typename FieldT>
r1cs_example<FieldT> generate_r1cs_example_with_matrix_opt_4(const size_t matrix_1_h,
                                                        const size_t matrix_1_w,
				                                        const size_t matrix_2_h,
                                                        const size_t matrix_2_w)
{
    libff::enter_block("Call to generate_r1cs_example_with_matrix_opt_4");
    r1cs_constraint_system<FieldT> cs;
    cs.primary_input_size = 0;
    cs.auxiliary_input_size = (matrix_1_h)*(matrix_1_w) + matrix_1_h*matrix_2_w*(matrix_1_w) + (matrix_2_w)*(matrix_2_h);
    
    size_t output_h = matrix_1_h;
    size_t output_w = matrix_2_w;

    size_t num_constraints = matrix_1_w;

    libff::enter_block("set variables");
    r1cs_variable_assignment<FieldT> full_variable_assignment;
    // push back input_x
    for(size_t i=0; i<(matrix_1_h)*(matrix_1_w); i++){
        full_variable_assignment.push_back(i+1);
    }
    // push back weight matrix
    for(size_t i=0; i<(matrix_2_w)*(matrix_2_h); i++){
        full_variable_assignment.push_back(i+1);
    }

    libff::enter_block("Compute y variables");
    // std::cout<< "Y-matrix" <<std::endl;
    for(size_t y_i=0; y_i<output_h; y_i++){
        for(size_t y_j=0; y_j<output_w; y_j++){
            FieldT y = FieldT::zero();
            for(size_t pos=0; pos<matrix_1_w; pos++){
                FieldT temp_y = FieldT::zero();
                temp_y = FieldT(y_i*matrix_1_w + pos + 1) * FieldT(pos*matrix_2_w + y_j + 1);
                y += temp_y;
                full_variable_assignment.push_back(y);
            }
        }
    }
    libff::leave_block("Compute y variables");

    libff::leave_block("set variables");

    libff::enter_block("set constraints");
    FieldT s = FieldT(2);
    FieldT t = qpow(s, matrix_2_w);

    for (size_t i = 0; i < matrix_1_w; ++i)
    {
        // FieldT t = qpow(s, matrix_2_w);
        linear_combination<FieldT> A, B, C;

        FieldT temp_t = FieldT::one();
        for(size_t x_i = 0; x_i < matrix_1_h; x_i ++){
            A.add_term((x_i*matrix_1_w)+i+1, temp_t);
            temp_t *= t;
        }
        FieldT temp_s = FieldT::one();
        for(size_t a_j = 0; a_j < matrix_2_w; a_j ++){
            B.add_term((matrix_1_h*matrix_1_w)+(i*matrix_2_w)+a_j+1, temp_s);
            // std::cout<<"a["<<i<<","<<a_j<<"]"<<"*"<<temp_s.as_ulong()<<" + ";
            temp_s *= s;
        }
        temp_s = FieldT::one();
        temp_t = FieldT::one();
        FieldT temp_minus_one = FieldT(-1);
        // std::cout<<"add for C"<<std::endl;
        for(size_t x_i=0; x_i<matrix_1_h; x_i++){
            for(size_t a_j=0; a_j<matrix_2_w; a_j++){
                C.add_term((matrix_1_h*matrix_1_w)+(matrix_2_h*matrix_2_w)+((x_i*output_w)+a_j)*(matrix_1_w)+i+1,temp_s*temp_t);
                if(i!=0){
                    C.add_term((matrix_1_h*matrix_1_w)+(matrix_2_h*matrix_2_w)+((x_i*output_w)+a_j)*(matrix_1_w)+i,temp_minus_one*temp_s*temp_t);
                }
                temp_s *= s;
            }
            temp_s = FieldT::one();
            temp_t *= t;
        }
        // std::cout<<std::endl;
        cs.add_constraint(r1cs_constraint<FieldT>(A, B, C));
        // std::cout<<std::endl;
    }

    libff::leave_block("set constraints");

    /* split variable assignment */
    r1cs_primary_input<FieldT> primary_input(full_variable_assignment.begin(), full_variable_assignment.begin() + cs.primary_input_size);
    r1cs_primary_input<FieldT> auxiliary_input(full_variable_assignment.begin() + cs.primary_input_size, full_variable_assignment.end());

    /* sanity checks */
    std::cout<<"Circuit parameters::  #constraints == "<<cs.num_constraints() << ";    #variables ==  "<<cs.num_variables()  <<std::endl;
    assert(cs.num_variables() == full_variable_assignment.size());
    assert(cs.num_variables() >= cs.primary_input_size);
    assert(cs.num_inputs() == cs.primary_input_size);
    assert(cs.num_constraints() == num_constraints);
    assert(cs.is_satisfied(primary_input, auxiliary_input));

    libff::leave_block("Call to generate_r1cs_example_with_matrix_opt_4");

    return r1cs_example<FieldT>(std::move(cs), std::move(primary_input), std::move(auxiliary_input));
}




/**
 * generate r1cs for a pooling operation
*/
template<typename FieldT>
r1cs_example<FieldT> generate_r1cs_example_from_pooling_matrix(const size_t input_h,
                                                        const size_t input_w,
				                                        const size_t kernel_h,
                                                        const size_t kernel_w)
{
    libff::enter_block("Call to generate_r1cs_example_from_pooling_matrix");
    
    r1cs_constraint_system<FieldT> cs;

    size_t divisor = kernel_h * kernel_w;
    cs.primary_input_size = 0;
    size_t output_h = (input_h / kernel_h);
    size_t output_w = (input_w / kernel_w);
    cs.auxiliary_input_size = output_h*output_w  + input_h*input_w;

    size_t num_constraints = output_h*output_w;

    libff::enter_block("set variables");
    
    r1cs_variable_assignment<FieldT> full_variable_assignment;

    std::vector<int> temp_vector;
    for(size_t x_i=0; x_i<input_h; x_i++){
        for(size_t x_j=0; x_j<input_w; x_j++){
            size_t x_index = x_i*input_w + x_j;
            temp_vector.push_back(x_index);
            std::cout<< temp_vector[x_index] <<" ";
        }
        std::cout<<std::endl;
    }

    // push back output matrix
    libff::enter_block("Compute y variables");
    for(size_t y_i=0; y_i<output_h; y_i++){
        for(size_t y_j=0; y_j<output_w; y_j++){
            size_t temp_y=0;
            for(size_t k_i=0; k_i<kernel_h; k_i++)
                for(size_t k_j=0; k_j<kernel_w; k_j++){
                    size_t x_index = (y_i*kernel_h + k_i)*input_w + (y_j*kernel_w+k_j);
                    temp_y += temp_vector[x_index];
                }
            full_variable_assignment.push_back(FieldT(temp_y));
        }
    }
    libff::leave_block("Compute y variables");
    // push back input_x
    for(size_t x_i=0; x_i<input_h; x_i++){
        for(size_t x_j=0; x_j<input_w; x_j++){
            size_t x_index = x_i*input_w + x_j;
            full_variable_assignment.push_back(FieldT(temp_vector[x_index]));
        }
    }

    libff::leave_block("set variables");
    
    libff::enter_block("set constraints");
    for(size_t y_i=0; y_i<output_h; y_i++){
        for(size_t y_j=0; y_j<output_w; y_j++){
            linear_combination<FieldT> A, B, C;
            for(size_t k_i=0; k_i<kernel_h; k_i++)
                for(size_t k_j=0; k_j<kernel_w; k_j++){
                    size_t x_index = (y_i*kernel_h + k_i)*input_w + (y_j*kernel_w+k_j);
                    A.add_term(output_h*output_w + x_index+1, 1);
                }
            B.add_term(0, 1);
            C.add_term(y_i*output_w + y_j + 1, 1);
            cs.add_constraint(r1cs_constraint<FieldT>(A, B, C));
        }
    }


    libff::leave_block("set constraints");

    /* split variable assignment */
    r1cs_primary_input<FieldT> primary_input(full_variable_assignment.begin(), full_variable_assignment.begin() + cs.primary_input_size);
    r1cs_primary_input<FieldT> auxiliary_input(full_variable_assignment.begin() + cs.primary_input_size, full_variable_assignment.end());

    /* sanity checks */
    std::cout<<"cs.num_cs = "<<cs.num_constraints() <<"      defined"<< num_constraints <<std::endl;
    std::cout<<"cs.num_variables = "<<cs.num_variables() <<"      defined: "<< full_variable_assignment.size() <<std::endl;
    std::cout<<"cs.num_inputs = "<<cs.num_inputs() <<"      defined"<< cs.primary_input_size <<std::endl;
    assert(cs.num_variables() == full_variable_assignment.size());
    assert(cs.num_variables() >= cs.primary_input_size);
    assert(cs.num_inputs() == cs.primary_input_size);
    assert(cs.num_constraints() == num_constraints);
    assert(cs.is_satisfied(primary_input, auxiliary_input));

    libff::leave_block("Call to generate_r1cs_example_from_pooling_matrix");

    return r1cs_example<FieldT>(std::move(cs), std::move(primary_input), std::move(auxiliary_input));
}


} // libsnark

#endif // R1CS_EXAMPLES_TCC
