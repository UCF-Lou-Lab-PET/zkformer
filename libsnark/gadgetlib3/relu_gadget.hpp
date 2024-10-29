#include <iostream>

#include <libsnark/gadgetlib1/gadget.hpp>
#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>
#include <libff/algebra/curves/public_params.hpp>
#include "libff/algebra/fields/field_utils.hpp"


namespace libsnark {

template <typename FieldT>
class ReLU_gadget : public gadget<FieldT>
{
private:
    pb_variable<FieldT> less_or_eq;
    std::shared_ptr<comparison_gadget<FieldT>> comparison;

public:
    pb_variable<FieldT> out;
    pb_variable<FieldT> in;
    pb_variable<FieldT> threshold;
    pb_variable<FieldT> less;
    ReLU_gadget(){};
    ReLU_gadget(protoboard<FieldT> &pb,
                const pb_variable<FieldT> &out,
                const pb_variable<FieldT> &in,
                const pb_variable<FieldT> &threshold,
                const std::string &annotation_prefix = "") : gadget<FieldT>(pb, annotation_prefix), out(out), in(in), threshold(threshold)
    {
        less.allocate(pb, FMT(this->annotation_prefix, "less"));
        less_or_eq.allocate(pb, FMT(this->annotation_prefix, "less_or_eq"));
        comparison.reset(new comparison_gadget<FieldT>(pb, 10, in, threshold, less, less_or_eq,
                                                       FMT(this->annotation_prefix, "comparison")));
    };

    // varaibles in the pb: [one, out, in, threshold, sym_1...] 
    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

template <typename FieldT>
void ReLU_gadget<FieldT>::generate_r1cs_constraints()
{
    comparison->generate_r1cs_constraints();

    this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(in, FieldT::one() - less, out),
                                 FMT(this->annotation_prefix, "y = x * (! (x<0) )"));
};
/*
compute the value of intermediate variables based on given in/out
*/
template <typename FieldT>
void ReLU_gadget<FieldT>::generate_r1cs_witness()
{
    comparison->generate_r1cs_witness();
};


/**
 * relu gadget for arrays/vectors
 * relu(pb, OUT, IN, threshold, "out ? max{in, 0}");
 * for the relu computation, we only care about whether x<0 or not, less_or_eq is not needed
*/
template <typename FieldT>
class ReLU_vector_gadget : public gadget<FieldT>
{
private:
    pb_variable_array<FieldT> less_or_eq;
    std::vector<comparison_gadget<FieldT>> comparison;

public:
    size_t n;
    pb_variable_array<FieldT> out;
    pb_variable_array<FieldT> in;
    pb_variable<FieldT> threshold;
    pb_variable_array<FieldT> less;

    ReLU_vector_gadget(protoboard<FieldT> &pb,
                        pb_variable_array<FieldT> &out,
                        pb_variable_array<FieldT> &in,
                        const pb_variable<FieldT> &threshold,
                        size_t n,
                        const std::string &annotation_prefix = "") : gadget<FieldT>(pb, annotation_prefix), out(out), in(in), threshold(threshold), n(n)
    {
        less.allocate(pb, n, FMT(this->annotation_prefix, "less"));
        less_or_eq.allocate(pb, n, FMT(this->annotation_prefix, "less_or_eq"));
        for (size_t i=0; i < n; i++)
        {
            // the second parameter means how many bits are A and B
            comparison.push_back(comparison_gadget<FieldT>(pb,  "avg pooling token mixer"));
        }
    };

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

template <typename FieldT>
void ReLU_vector_gadget<FieldT>::generate_r1cs_constraints()
{

    for(size_t i = 0; i < n; i++){
        comparison[i].generate_r1cs_constraints();
        this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(in[i], FieldT::one() - less[i], out[i]),
                                    FMT(this->annotation_prefix, "y[i] = x[i] * (! (x[i]<0) )"));
    }
};

template <typename FieldT>
void ReLU_vector_gadget<FieldT>::generate_r1cs_witness()
{
    for(size_t i = 0; i < n; i++){
        comparison[i].generate_r1cs_witness();
    }
};





// writing a new relu which functions the same but work better with other components
template <typename FieldT>
class relu_gadget : public gadget<FieldT>
{
private:
    pb_variable_array<FieldT> less_or_eq;
    std::vector<comparison_gadget<FieldT>> comparison;

public:
    pb_variable_array<FieldT> in;
    pb_variable_array<FieldT> out;
    pb_variable<FieldT> threshold;
    size_t n;
    pb_variable_array<FieldT> less;

    relu_gadget(protoboard<FieldT> &pb,
                        pb_variable_array<FieldT> &in,
                        pb_variable_array<FieldT> &out,
                        const pb_variable<FieldT> &threshold,
                        size_t n,
                        const std::string &annotation_prefix = "") : gadget<FieldT>(pb, annotation_prefix), in(in), out(out), threshold(threshold), n(n)
    {
        less.allocate(pb, n, FMT(this->annotation_prefix, "'less"));
        less_or_eq.allocate(pb, n, FMT(this->annotation_prefix, "'less_or_eq"));
        for (size_t i=0; i < n; i++)
        {
            // the second parameter (10) means how many bits are A and B
            comparison.push_back(comparison_gadget<FieldT>(pb, 32, in[i], threshold, less[i], less_or_eq[i], "'comparison"));
        }
    };

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

template <typename FieldT>
void relu_gadget<FieldT>::generate_r1cs_constraints()
{

    for(size_t i = 0; i < n; i++){
        comparison[i].generate_r1cs_constraints();
        this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(in[i], FieldT::one() - less[i], out[i]),
                                    FMT(this->annotation_prefix, "y[i] = x[i] * (! (x[i]<0) )"));
    }
};

template <typename FieldT>
void relu_gadget<FieldT>::generate_r1cs_witness()
{
    for(size_t i = 0; i < n; i++){
        comparison[i].generate_r1cs_witness();
    }
    for(size_t i = 0; i < n; i++){
        this->pb.val(out[i]) = (this->pb.val(in[i]).as_ulong() <= 0) ? FieldT::zero() : this->pb.val(in[i]);
    }
};


template <typename FieldT>
class relu_approx_gadget : public gadget<FieldT>
{
private:
//no private
    size_t length;
public:
    pb_tensor<FieldT> in;
    pb_tensor<FieldT> out;

    relu_approx_gadget(protoboard<FieldT> &pb,
                        pb_tensor<FieldT> &in,
                        pb_tensor<FieldT> &out,
                        const std::string &annotation_prefix = "") : gadget<FieldT>(pb, annotation_prefix), in(in), out(out)
    {
        length = in.size();
    };

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

template <typename FieldT>
void relu_approx_gadget<FieldT>::generate_r1cs_constraints()
{

    for(size_t i = 0; i < length; i++){
        this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(in[i], FieldT::one() + in[i], out[i]),
                                    FMT(this->annotation_prefix, "y[i] = x[i] * (x[i] + 1) )"));
    }
};

template <typename FieldT>
void relu_approx_gadget<FieldT>::generate_r1cs_witness()
{
    for(size_t i = 0; i < length; i++){
        this->pb.val(out[i]) = this->pb.val(in[i]) * (this->pb.val(in[i]) + FieldT::one());
    }
};













/*an proximation of sigmoid function
    0.5 + 0.2159198015 - 0.0082176259 + 0.0001825597 - 0.0000018848 + 0.0000000072
*/
template <typename FieldT>
class sigmoid_approx_gadget : public gadget<FieldT>
{
private:
//no private
    size_t length;
    pb_tensor<FieldT> in_2;
    pb_tensor<FieldT> in_3;
    pb_tensor<FieldT> in_5;
    pb_tensor<FieldT> in_7;
    pb_tensor<FieldT> in_9;

public:
    pb_tensor<FieldT> in;
    pb_tensor<FieldT> out;

    sigmoid_approx_gadget(protoboard<FieldT> &pb,
                        pb_tensor<FieldT> &in,
                        pb_tensor<FieldT> &out,
                        const std::string &annotation_prefix = "") : gadget<FieldT>(pb, annotation_prefix), in(in), out(out)
    {
        length = in.size();
        in_2.allocate(pb, length, " input[i] ^ 2 ");
        in_3.allocate(pb, length, " input[i] ^ 3 ");
        in_5.allocate(pb, length, " input[i] ^ 5 ");
        in_7.allocate(pb, length, " input[i] ^ 7 ");
        in_9.allocate(pb, length, " input[i] ^ 9 ");
    };

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

template <typename FieldT>
void sigmoid_approx_gadget<FieldT>::generate_r1cs_constraints()
{
    for(size_t i = 0; i < length; i++){
        
        linear_combination<FieldT> A0, B0, C0;
        A0.add_term(in[i], 1);
        B0.add_term(in[i], 1);
        C0.add_term(in_2[i], 1);
        this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(A0, B0, C0), " orders ");

        linear_combination<FieldT> A1, B1, C1;
        A1.add_term(in_2[i], 1);
        B1.add_term(in[i], 1);
        C1.add_term(in_3[i], 1);
        this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(A1, B1, C1), " orders ");

        linear_combination<FieldT> A2, B2, C2;
        A2.add_term(in_2[i], 1);
        B2.add_term(in_3[i], 1);
        C2.add_term(in_5[i], 1);
        this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(A2, B2, C2), " orders ");

        linear_combination<FieldT> A3, B3, C3;
        A3.add_term(in_2[i], 1);
        B3.add_term(in_5[i], 1);
        C3.add_term(in_7[i], 1);
        this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(A3, B3, C3), " orders ");


        linear_combination<FieldT> A4, B4, C4;
        A4.add_term(in_2[i], 1);
        B4.add_term(in_7[i], 1);
        C4.add_term(in_9[i], 1);
        this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(A4, B4, C4), " orders ");


        linear_combination<FieldT> A, B, C;
        A.add_term(0, FieldT(500000000));
        A.add_term(in[i], FieldT(2159198015));
        A.add_term(in_3[i], FieldT(821762590) * FieldT(-1));
        A.add_term(in_5[i], FieldT(1825597));
        A.add_term(in_7[i], FieldT(18848) * FieldT(-1));
        A.add_term(in_9[i], FieldT(72));
        B.add_term(0, 1);
        C.add_term(out[i], FieldT(10000000000));
        this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(A, B, C), " sigmoid ");
    }
};

template <typename FieldT>
void sigmoid_approx_gadget<FieldT>::generate_r1cs_witness()
{
    FieldT divisor = FieldT(10000000000).inverse();
    for(size_t i = 0; i < length; i++){

        this->pb.val(in_2[i]) =  this->pb.val(in[i]) * this->pb.val(in[i]);
        this->pb.val(in_3[i]) =  this->pb.val(in[i]) * this->pb.val(in_2[i]);
        this->pb.val(in_5[i]) =  this->pb.val(in_3[i]) * this->pb.val(in_2[i]);
        this->pb.val(in_7[i]) =  this->pb.val(in_5[i]) * this->pb.val(in_2[i]);
        this->pb.val(in_9[i]) =  this->pb.val(in_7[i]) * this->pb.val(in_2[i]);

        this->pb.val(out[i]) = 
            FieldT(500000000) * divisor
            + FieldT(2159198015) * divisor *  this->pb.val(in[i])
            - FieldT(821762590) * divisor * this->pb.val(in_3[i])
            + FieldT(1825597) * divisor * this->pb.val(in_5[i])
            - FieldT(18848) * divisor * this->pb.val(in_7[i])
            + FieldT(72) * divisor * this->pb.val(in_9[i]);
    }
};





/**
 * create a look up table in relu gadget
 * st. all comparison gadget share the same space for alpha
 * this doesn't work out
*/
// template <typename FieldT>
// class ReLU_vector_gadget_opt : public gadget<FieldT>
// {
// private:
//     pb_variable_array<FieldT> less_or_eq;
//     pb_variable_array<FieldT> alpha_table;
//     // std::vector<comparison_gadget<FieldT>> comparison;
//     std::vector<comparison_gadget_opt<FieldT>> comparison;

// public:
//     size_t n;
//     pb_variable_array<FieldT> out;
//     pb_variable_array<FieldT> in;
//     pb_variable<FieldT> threshold;
//     pb_variable_array<FieldT> less;

//     ReLU_vector_gadget_opt(protoboard<FieldT> &pb,
//                         pb_variable_array<FieldT> &out,
//                         pb_variable_array<FieldT> &in,
//                         const pb_variable<FieldT> &threshold,
//                         size_t n,
//                         const std::string &annotation_prefix = "") : gadget<FieldT>(pb, annotation_prefix), out(out), in(in), threshold(threshold), n(n)
//     {
//         less.allocate(pb, n, FMT(this->annotation_prefix, "less"));
//         less_or_eq.allocate(pb, n, FMT(this->annotation_prefix, "less_or_eq"));
//         alpha_table.allocate(pb, 2*n, FMT(this->annotation_prefix, "alpha look up table"));
//         for (size_t i=0; i < n; i++)
//         {
//             // method1: let comparison gadget creat alpha
//             // comparison.push_back(comparison_gadget<FieldT>(pb, 10, in[i], threshold, less[i], less_or_eq[i], "comparison"));
//             // method2: use pre-defined alpha
//             comparison.push_back(comparison_gadget<FieldT>(pb, 10, in[i], threshold, less[i], less_or_eq[i], alpha_table, "comparison"));
//         }
//     };

//     void generate_r1cs_constraints();
//     void generate_r1cs_witness();
// };

// template <typename FieldT>
// void ReLU_vector_gadget_opt<FieldT>::generate_r1cs_constraints()
// {
//     for(size_t i = 0; i < n; i++){
//         comparison[i].generate_r1cs_constraints();
//         this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(in[i], FieldT::one() - less[i], out[i]),
//                                     FMT(this->annotation_prefix, "y[i] = x[i] * (! (x[i]<0) )"));
//     }
// }

// template <typename FieldT>
// void ReLU_vector_gadget_opt<FieldT>::generate_r1cs_witness()
// {
//     for(size_t i = 0; i < 2*n; i+=2){
//         pb.val(alpha_table[i]) = FieldT::zero();
//         pb.val(alpha_table[i+1]) = FieldT::one();
//     }
//     for(size_t i = 0; i < n; i++){
//         comparison[i].generate_r1cs_witness();
//     }
// }






/**
 * ReLU_vector_gadget with more constraints
 * namely, y = x * (1-(x<a)) + a * (x<a)
 *  if x<a, y == a
 *  else y==x
 * where in relu, a(threshold) is set to 0 
*/
template <typename FieldT>
class ReLU_vector_gadget_wmc : public gadget<FieldT>
{
private:
    pb_variable_array<FieldT> less_or_eq;
    std::vector<comparison_gadget<FieldT>> comparison;

public:
    size_t n;
    pb_variable_array<FieldT> out;
    pb_variable_array<FieldT> in;
    pb_variable<FieldT> threshold;
    pb_variable_array<FieldT> sym1;
    pb_variable_array<FieldT> sym2;
    pb_variable_array<FieldT> less;

    ReLU_vector_gadget_wmc(protoboard<FieldT> &pb,
                        pb_variable_array<FieldT> &out,
                        pb_variable_array<FieldT> &in,
                        const pb_variable<FieldT> &threshold,
                        size_t n,
                        const std::string &annotation_prefix = "") : gadget<FieldT>(pb, annotation_prefix), out(out), in(in), threshold(threshold), n(n)
    {
        // intermediate syms
        sym1.allocate(pb, n, FMT(this->annotation_prefix, "x*(!x<a)"));
        sym2.allocate(pb, n, FMT(this->annotation_prefix, "a*(x<a)"));
        less.allocate(pb, n, FMT(this->annotation_prefix, "less"));
        less_or_eq.allocate(pb, n, FMT(this->annotation_prefix, "less_or_eq"));
        for (size_t i=0; i < n; i++)
        {
            comparison.push_back(comparison_gadget<FieldT>(pb, 10, in[i], threshold, less[i], less_or_eq[i], "comparison"));
        }
    };

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

template <typename FieldT>
void ReLU_vector_gadget_wmc<FieldT>::generate_r1cs_constraints()
{

    for(size_t i = 0; i < n; i++){
        comparison[i].generate_r1cs_constraints();
        this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(in[i], FieldT::one() - less[i], sym1[i]),
                                    FMT(this->annotation_prefix, "x[i] * (! (x[i]<a) )"));
        this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(threshold, less[i], sym2[i]),
                                    FMT(this->annotation_prefix, "a * (x[i]<0) "));
        this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(sym1[i] + sym2[i], 1, out[i]),
                                    FMT(this->annotation_prefix, "y[i] = x[i] * (! (x[i]<0) ) + a * (x[i]<0)"));
        // this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(in[i], FieldT::one() - less[i], out[i]),
        //                             FMT(this->annotation_prefix, "y[i] = x[i] * (! (x[i]<0) )"));
    }
}

template <typename FieldT>
void ReLU_vector_gadget_wmc<FieldT>::generate_r1cs_witness()
{
    for(size_t i = 0; i < n; i++){
        comparison[i].generate_r1cs_witness();
        this->pb.val(sym1[i]) = this->pb.val(in[i]) * (FieldT::one() - this->pb.val(less[i]));
        this->pb.val(sym2[i]) = this->pb.val(threshold) * this->pb.val(less[i]);
    }
}



}// libsnark