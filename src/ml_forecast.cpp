#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector exact_rf_forecast(NumericMatrix embed_mat, int h, int lags,
                              Function predict_func, List model) {
    int n = embed_mat.nrow();
    NumericVector forecasts(h);
    
    // Create full working matrix (matches R's growing approach)
    NumericMatrix working_mat(h + n, lags + 1);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= lags; j++) {
            working_mat(h + i, j) = embed_mat(i, j);
        }
    }
    
    // Create lag names identical to R
    CharacterVector lag_names(lags);
    for (int j = 0; j < lags; j++) {
        lag_names[j] = "lag" + std::to_string(lags - j);
    }
    
    for (int i = 0; i < h; i++) {
        // Create prediction data frame (exact match to R)
        List newdata(lags);
        for (int j = 0; j < lags; j++) {
            newdata[j] = NumericVector::create(working_mat(h - i, j));
        }
        newdata.attr("names") = lag_names;
        newdata.attr("class") = "data.frame";
        newdata.attr("row.names") = 1;
        
        // Make prediction
        forecasts[i] = as<double>(predict_func(model, newdata));
        
        // Update working matrix (EXACT match to R's behavior)
        // 1. Create new row with shifted lags
        for (int j = lags; j > 1; j--) {
            working_mat(h - i - 1, j - 1) = working_mat(h - i, j - 2);
        }
        // 2. Insert forecast as first lag
        working_mat(h - i - 1, 0) = forecasts[i];
        // 3. Set y value
        working_mat(h - i - 1, lags) = forecasts[i];
    }
    
    // Return forecasts in reverse order to match R
    std::reverse(forecasts.begin(), forecasts.end());
    return forecasts;
}