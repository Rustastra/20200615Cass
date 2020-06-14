library(data.table)
library(xlsx)
library(rstan)
library(dynpanel)
library(betareg)
library(plm)

#TODO: list of libraries; 
#details of PGMM
#Bayesian illustration#why we need hierarchical models
#quick intro into bayes: how it is different from freq. 
#hierarchical as a bonus


folder = "/home/labs/PROJECTS/Cass/git/"
dtlect = xlsx::read.xlsx(paste0(folder,"WEO_Data.xlsx"), sheetName = "WEO_Data", stringsAsFactors = F)
setDT(dtlect)
dtlect1 = dtlect[Subject.Descriptor %in% 
			c("Output gap in percent of potential GDP",
				"Gross domestic product, constant prices",
				"Current account balance",
				"General government primary net lending/borrowing"), ]

dtlect_melt = melt(
	dtlect1, 
	id.vars = c("Country", "Subject.Descriptor"), 
	measure.vars = c(paste0("X", 1980:2020)), 
	value.name="Value")

dtlect_use = dcast(dtlect_melt, "Country + variable ~ Subject.Descriptor", fun_aggregate = "sum")
colnames(dtlect_use)[2] = "Year"
colnames(dtlect_use)[3:6] = c("cur_ac", "deficit", "gdp","gap")

setDT(dtlect_use)
dtlect_use[, Year := substr(Year, 2, 5)]
dtlect_use[, Year := as.numeric(Year)]
dtlect_use[, gap:=as.numeric(gap)]
dtlect_use[, gdp:=as.numeric(gdp)]
dtlect_use[, deficit:=as.numeric(deficit)]
dtlect_use[, cur_ac:=as.numeric(cur_ac)]



setorder(dtlect_use, Country, Year)

dtlect_use[, gdp_lag1 := shift(gdp, type="lag", n = 1), by = Country]
dtlect_use[, curac_lag1 := shift(cur_ac, type="lag", n = 1), by = Country]
dtlect_use[ , gdp_lag1_sh := gdp_lag1/100]
dtlect_use[ , curac_lag1_sh := curac_lag1/100]
dtlect_use[ , curac_sh := cur_ac/100]


# =================================
# DATA DISTRIBUTION
# =================================
dtlect = dtlect_use[gdp < 0 & !is.na(gdp_lag1)]
dtlect[, gdp:= abs(gdp)]
dtlect[ , gdp_sh := gdp/100]

hist(dtlect[,gdp], breaks = 20)
plot(dtlect[,.(curac_lag1, gdp)])


# =================================
# NAIVE APPROACH
# =================================
model = lm("gdp_sh ~ gdp_lag1_sh + curac_lag1_sh", data = dtlect)
summary(model)


# =================================
# LOG-ODDS
# =================================
dtlect[,gdp_logodds := log(gdp_sh/(1-gdp_sh))]
model = lm("gdp_logodds ~ gdp_lag1_sh + curac_lag1_sh", 
	data = dtlect)
summary(model)



# =================================
# Beta regression
# =================================

# betareg(formula, data, subset, na.action, weights, offset,
# link = c("logit", "probit", "cloglog", "cauchit", "log", "loglog"),
# link.phi = NULL, type = c("ML", "BC", "BR"),
# control = betareg.control(...), model = TRUE,
# y = TRUE, x = FALSE, ...)

fla = "gdp_sh ~ gdp_lag1_sh + curac_lag1_sh"
model = betareg(fla, data = dtlect, link = "logit")
summary(model)

model = betareg(fla, data = dtlect, link = "probit")
summary(model)

model = betareg(fla, data = dtlect, link = "cloglog")
summary(model)

model = betareg(fla, data = dtlect, link = "cauchit")
summary(model)

model = betareg(fla, data = dtlect, link = "log")
summary(model)

model = betareg(fla, data = dtlect, link = "loglog")
summary(model)

# =================================
# German tank simulation
# =================================
# let's look at the distribution of sample max for N = 1000 and k = 100
N=1000;k=100
hist( unlist(lapply(1:50000, function(x) max(sample(1:N, k))) ), breaks = N)

N=2000;k=100
hist( unlist(lapply(1:50000, function(x) max(sample(1:N, k))) ), breaks = N)

#given k, how likely are certain N?
out = c()
for(n in k:(k*10)){
	out = c(out, dbinom(k, size=n, prob= 0.2)  )  #choose(n, k) * (1/n)^k )
}

plot(k:(k*10), out/sum(out))


# =================================
# BETA DISTRIBUTION CHARTS
# =================================


# Beta distributions
alpha = 0; beta = 0
hist(rbeta(1000, alpha, beta), col='blue', breaks = 100, main = paste0("Alpha: ", alpha, "; Beta: ", beta))

alpha = 1; beta = 1
hist(rbeta(50000, alpha, beta), col='blue', breaks = 50, main = paste0("Alpha: ", alpha, "; Beta: ", beta))

alpha = 10; beta = 1
hist(rbeta(50000, alpha, beta), col='blue', breaks = 50, main = paste0("Alpha: ", alpha, "; Beta: ", beta), xlim = c(0, 1))

alpha = 100; beta = 1
hist(rbeta(50000, alpha, beta), col='blue', breaks = 50, main = paste0("Alpha: ", alpha, "; Beta: ", beta), xlim = c(0, 1))

alpha = 1; beta = 10
hist(rbeta(50000, alpha, beta), col='blue', breaks = 50, main = paste0("Alpha: ", alpha, "; Beta: ", beta), xlim = c(0, 1))

alpha = 1; beta = 100
hist(rbeta(50000, alpha, beta), col='blue', breaks = 50, main = paste0("Alpha: ", alpha, "; Beta: ", beta), xlim = c(0, 1))

alpha = 100; beta = 100
hist(rbeta(50000, alpha, beta), col='blue', breaks = 50, main = paste0("Alpha: ", alpha, "; Beta: ", beta), xlim = c(0, 1))

alpha = 2; beta = 5
hist(rbeta(50000, alpha, beta), col='blue', breaks = 50, main = paste0("Alpha: ", alpha, "; Beta: ", beta), xlim = c(0, 1))

alpha = 0.1; beta = 0.1
hist(rbeta(50000, alpha, beta), col='blue', breaks = 50, main = paste0("Alpha: ", alpha, "; Beta: ", beta), xlim = c(0, 1))





# =================================
# Arellano-Bond GMM
# =================================


dta = dtlect_use
dta[ , gdp_sh := gdp/100]
dta[ , curac_sh := cur_ac/100]
dta[,gdp_logodds := log(gdp_sh/(1-gdp_sh))]

dta[,index := 1:.N,by=Country]
dta = dta[,.(Country, index, gdp_logodds, curac_sh, gdp)]

dta = dta[abs(gdp_logodds)!=Inf]

z1 <- pgmm(gdp_logodds ~ lag(gdp_logodds,1) + lag(curac_sh,1) | lag(gdp_logodds,2),
            data = dta, subset=(gdp<0), transformation = "d")
summary(z1)




# =================================
# Bayesian approach
# =================================
dta = dtlect[!is.na(gdp_sh + curac_lag1_sh + gdp_lag1_sh),
				.(Country,Year,gap, gdp_sh, curac_lag1_sh, gdp_lag1_sh)]
dta[,intercept:=1]
data_list = list(
	y = dta$gdp_sh,
	x = t(as.matrix(data.frame(
		"intercept" = dta$intercept,
		"curac" = dta$curac_lag1_sh,
		"y_lag" = dta$gdp_lag1_sh
		))),
	N = nrow(dta)
	)


model_string = "
	data {
		int N;
		real y[N];
		matrix[3, N] x;
	}

	parameters {
		real<lower = 0> phi;
		row_vector[3] coefs;
	}

	transformed parameters {
		vector[N] mu;
		vector[N] alpha_;
		vector[N] beta_;

		for(i in 1:N){
			mu[i] = 1/(1+exp(-coefs*x[,i]));
			alpha_[i] = mu[i]*phi;
			beta_[i] = (1-mu[i])*phi;
		}
	}

	model {
		for(i in 1:N){
			y[i] ~ beta(alpha_[i], beta_[i]);
		}
		coefs[1] ~ normal(0,10);
		coefs[2] ~ normal(0,10);
		coefs[3] ~ normal(0,10);
		phi ~ gamma(.1,.1);
	}
"

# takes some time to compile 
# about 3 minutes
model = stan_model(model_code = model_string)


samples = sampling(model, iter = 1500, data = data_list, chains = 5) #better to have 5k iterations and 16 chains
traceplot(samples, pars=c("phi","coefs[3]"))

#dim (df_1)
df_1 = as.data.frame(samples)
hist(df_1["coefs[1]"][[1]])









# =================================
# HIERARCHICAL Bayesian approach
# =================================
dta = dtlect[!is.na(gdp_sh + curac_lag1_sh + gdp_lag1_sh),
				.(Country,Year,gap, gdp_sh, curac_lag1_sh, gdp_lag1_sh)]
dta[,intercept:=1]
dta[, is_WEurope := 1]
dta[Country %in% 
		c("Belgium", 
			"Austria", 
			"Germany", 
			"France", 
			"Italy", 
			"Netherlands",
			"Finland",
			"Sweden",
			"Denmark",
			"Norway"), is_WEurope := 2]
data_list = list(
	y = dta$gdp_sh,
	curac = dta$curac_lag1_sh,
	y_lag = dta$gdp_lag1_sh,
	is_Europe = as.vector(dta$is_WEurope),
	K = 2, 
	Cats = c(1, 2),
	N = nrow(dta)
	)


model_string_hier = "
	data {
		int N;
		int K;
		int<lower=1,upper=2> is_Europe[N];
		real y[N];
		real y_lag[N];
		real curac[N];
	}

	parameters {
		real<lower = 0> phi;
		vector[K] b_intercept;
		real b_curac;
		vector[K] b_y_lag;
	}

	transformed parameters {
		vector[N] mu;
		vector[N] alpha_;
		vector[N] beta_;

		for(i in 1:N){
			mu[i] = 1/(1+exp(-( b_intercept[is_Europe[i]] + b_curac*curac[i] + b_y_lag[is_Europe[i]]*y_lag[i])));
			alpha_[i] = mu[i]*phi; 
			beta_[i] = (1-mu[i])*phi;
		}
	}

	model {
		for(i in 1:N){
			y[i] ~ beta(alpha_[i], beta_[i]);
		}
		b_curac ~ normal(0,100);
		for(i in 1:K){
			b_y_lag[i] ~ normal(0,100);
		}
		for(i in 1:K){
			b_intercept[i] ~ normal(0,100);
		}
		phi ~ gamma(.1,.1);
	}
"

model_hier = stan_model(model_code = model_string_hier)
samples_hier = sampling(model_hier, iter = 1000, data = data_list, chains = 5)
traceplot(samples_hier, pars=c("phi","b_y_lag[1]","b_y_lag[2]"))

traceplot(samples_hier, pars=c("phi","b_intercept[1]","b_intercept[2]"))

df_hier = as.data.frame(samples_hier)
hist(as.numeric(df_hier["b_y_lag[1]"][[1]]))
hist(as.numeric(df_hier["b_y_lag[2]"][[1]]))


pairs(samples_hier, pars=c("phi","b_y_lag[1]","b_y_lag[2]"))


