data {
    int<lower=1> N;                   // # observations
    int<lower=2> P;                   // # populations
    vector[N] fst;                    // pairwise Fst
    vector[N] dist;                   // pairwise distances
    int<lower=1, upper=P> pop1[N];    // population id for first of pair
    int<lower=1, upper=P> pop2[N];    // population id for second of pair
}
parameters {
    real alpha[2];          // fixed intercept and slope
    vector[P] pop_int;      // population intercepts
    real<lower=0> sig_res;  // residual sd
    real<lower=0> sig_pop;  // population sd
}
transformed parameters {
    vector[N] fst_pred;     // predicted Fst values
    for (i in 1:N){
        fst_pred[i] = alpha[1] + pop_int[pop1[i]] + pop_int[pop2[i]] + alpha[2] * dist[i];
    }
}
model {
    alpha ~ normal(0, 1);
    sig_res ~ gamma(1.5, 3);
    sig_pop ~ gamma(1.5, 3);
    pop_int ~ normal(0, sig_pop);
    fst ~ normal(fst_pred, sig_res);
}
generated quantities {
    real log_lik[N];
    real log_lik_sum;
    for (i in 1:N){
        log_lik[i] = normal_lpdf(fst[i] | fst_pred[i], sig_res);
    }
    log_lik_sum = sum(log_lik);
}
