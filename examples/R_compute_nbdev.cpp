#include "utils.h"
#include "glm.h"
#include "objects.h"

SEXP compute_nbdev (SEXP y, SEXP mu, SEXP phi, SEXP weights, SEXP dosum) {
    BEGIN_RCPP
    any_numeric_matrix counts(y);               //counts定义类型为any_numeric_matrix，取值为y
    const int num_tags=counts.get_nrow();      //nrow
    const int num_libs=counts.get_ncol();      //ncol
    std::vector<double> current(num_libs);      // std::vector<double> current = num_libs

    // Setting up means.
    Rcpp::NumericMatrix fitted(mu);             //fitted定义类型为NumericMatrix，取值为mu
    if (fitted.nrow()!=num_tags || fitted.ncol()!=num_libs) {
        throw std::runtime_error("dimensions of count and fitted value matrices are not equal");
    }

    // Setting up dispersions.
    compressed_matrix alld=check_CM_dims(phi, num_tags, num_libs, "dispersion", "count");

    // Seeing if we have to sum things together.
    bool sumtogether=check_logical_scalar(dosum, "summation specifier");

    if (sumtogether) {
        // Setting up weights.
        compressed_matrix allw(weights);      //等同于 compressed_matrix allw = weights;

        Rcpp::NumericVector output(num_tags);
        for (int tag=0; tag<num_tags; ++tag) {       //逐个读取基因
            counts.fill_row(tag, current.data());       //auto current=imat.row(index); std::copy(current.begin(), current.end(), ptr) 将counts对应tag行的内容写入current
            const double* dptr=alld.get_row(tag);     //提取dispersion
            const double* wptr=allw.get_row(tag);     //提取weight

            Rcpp::NumericMatrix::Row curmeans=fitted.row(tag);
            double& current_sumdev=output[tag];           //将current_sumdev指针指向output[tag]位置 vector长度可变
            auto cmIt=curmeans.begin();                   //cmIt为内存位置
            for (int lib=0; lib<num_libs; ++lib, ++cmIt) {      //逐个读取样本,current[lib]指各样本在此基因的count数
                current_sumdev += compute_unit_nb_deviance(current[lib], *cmIt, dptr[lib]) * wptr[lib];
            }
        }

        return output;
    } else {
        Rcpp::NumericMatrix output(num_tags, num_libs);
        for (int tag=0; tag<num_tags; ++tag) {
            counts.fill_row(tag, current.data());
            const double* dptr=alld.get_row(tag);

            auto curmeans=fitted.row(tag);
            auto cmIt=curmeans.begin();
            auto outvals=output.row(tag);
            auto ovIt=outvals.begin();

            for (int lib=0; lib<num_libs; ++lib, ++ovIt, ++cmIt) {
                (*ovIt) = compute_unit_nb_deviance(current[lib], *cmIt, dptr[lib]);
            }
       }

        return output;
    }
    END_RCPP
}

